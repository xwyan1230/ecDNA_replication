# scATAC-seq data processing
# https://muon-tutorials.readthedocs.io/en/latest/single-cell-rna-atac/pbmc10k/2-Chromatin-Accessibility-Processing.html

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
# Import a module with ATAC-seq-related functions
from muon import atac as ac
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# This is the directory where those files are downloaded to
data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221102_scMultiome_learning/pbmc10k/03_analysis/'
output_dir = data_dir
os.chdir(data_dir)

mdata = mu.read("pbmc10k.h5mu")
print(mdata)

# #### ATAC
atac = mdata.mod['atac']

# ### PREPROCESSING
# ## QC
sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
print(atac)
sc.pl.violin(atac, ['total_counts', 'n_genes_by_counts'], jitter=0.4, multi_panel=True)

# Filter peaks which expression is not detected:
mu.pp.filter_var(atac, 'n_cells_by_counts', lambda x: x >= 10)
# This is analogous to
#   sc.pp.filter_genes(rna, min_cells=10)
# but does in-place filtering and avoids copying the object

# filter cells
mu.pp.filter_obs(atac, 'n_genes_by_counts', lambda x: (x >= 2000) & (x <= 15000))
# This is analogous to
#   sc.pp.filter_cells(atac, max_genes=15000)
#   sc.pp.filter_cells(atac, min_genes=2000)
# but does in-place filtering avoiding copying the object
mu.pp.filter_obs(atac, 'total_counts', lambda x: (x >= 4000) & (x <= 40000))
sc.pl.violin(atac, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True)
mu.pl.histogram(atac, ['n_genes_by_counts', 'total_counts'])

# ## ATAC specific QC
# There are a few expectations about how ATAC-seq data looks like as noted in the hitchhiker’s guide to ATAC-seq
# data analysis for instance.
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3#citeas

# # nucleosome signal
# Fragment size distribution typically reflects nucleosome binding pattern showing enrichment around values
# corresponding to fragments bound to a single nucleosome (between 147 bp and 294 bp) as well as nucleosome-free
# fragments (shorter than 147 bp).
atac.obs['NS'] = 1
print(atac.uns)
ac.pl.fragment_histogram(atac, region='chr1:1-2000000')
# The ratio of mono-nucleosome cut fragments to nucleosome-free fragments can be called nucleosome signal, and
# it can be estimated using a subset of fragments.
ac.tl.nucleosome_signal(atac, n=1e6)
mu.pl.histogram(atac, "nucleosome_signal", kde=False)

# # TSS enrichment
# We can expect chromatin accessibility enriched around transcription start sites (TSS) compared to accessibility
# of flanking regions. Thus this measure averaged across multiple genes can serve as one more quality control metric.
# The positions of transcription start sites can be obtained from the interval field of the gene annotation in the
# rna modality
# TSS enrichment function will return an AnnData object with cells x bases dimensions where bases correspond to
# positions around TSS and are defined by extend_upstream and extend_downstream parameters, each of them being 1000
# bp by default. It will also record tss_score in the original object.
tss = ac.tl.tss_enrichment(mdata, n_tss=1000)  # by default, features=ac.tl.get_gene_annotation_from_rna(mdata)
print(tss)
ac.pl.tss_enrichment(tss)

print(mdata)

# ## normalization
# Save original counts
atac.layers["counts"] = atac.X

# There can be multiple options for ATAC-seq data normalisation.

# One is latent semantic indexing that is frequently used for processing ATAC-seq datasets. First, it constructs
# term-document matrix from the original count matrix. Then the singular value decomposition (SVD) — the same
# technique that convential principal component analysis uses — is used to generate LSI components. Note that
# there are different flavours of computing TF-IDF, e.g. see this blog post about that.
# http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/
#
# TF-IDF normalisation is implemented in the muon’s ATAC module:
ac.pp.tfidf(atac, scale_factor=1e4)

# Here we will use the same log-normalisation and PCA that we are used to from scRNA-seq analysis. We notice
# on this data it yields PC & UMAP spaces similar to the one generated on scRNA-seq counts.
sc.pp.normalize_per_cell(atac, counts_per_cell_after=1e4)
sc.pp.log1p(atac)

# ## feature selection
# We will label highly variable peaks that we’ll use for downstream analysis.
sc.pp.highly_variable_genes(atac, min_mean=0.05, max_mean=1.5, min_disp=.5)
sc.pl.highly_variable_genes(atac)
print(np.sum(atac.var.highly_variable))

# ## scaling
# For uniformity, and for consequent visualisation, we’ll save log-transformed counts in a .raw slot:
atac.raw = atac

# ### ANALYSIS
# After filtering out low-quality cells, normalising the counts matrix, and selecting highly varianbe peaks,
# we can already use this data for multimodal integration.
#
# However, as in the case of gene expression, we will study this data individually first and will run PCA on
# the scaled matrix, compute cell neighbourhood graph, and perform clustering to define cell types. This might
# be useful later to compare cell type definition between modalities.

# ## LSI
# When working on TF-IDF counts, sc.tl.pca or ac.tl.lsi can be used to get latent components, e.g.:
ac.tl.lsi(atac)
# We find the first component is typically associated with number of peaks or counts per cell so it is reasonable
# to remove it:
atac.obsm['X_lsi'] = atac.obsm['X_lsi'][:, 1:]
atac.varm["LSI"] = atac.varm["LSI"][:, 1:]
atac.uns["lsi"]["stdev"] = atac.uns["lsi"]["stdev"][1:]
# The respective neighbourhood graph can be generated with sc.tl.neighbors:
sc.pp.neighbors(atac, use_rep="X_lsi", n_neighbors=10, n_pcs=30)

# ## PCA
# For this notebook, we are using PCA on the log-normalised counts in atac.X as described above.
sc.pp.scale(atac)
sc.tl.pca(atac)
# We can only colour our plots by cut counts in individual peaks with scanpy:
sc.pl.pca(atac, color=["n_genes_by_counts", "n_counts"])

# With muon’s ATAC module, we can plot average values for cut counts in peaks of different types (promoter/distal)
# that are assigned to respective genes — just by providing gene names.
#
# For that to work, we need the peak annotation table with gene -> peak correspondence. The peak_annotation.tsv
# file was detected and loaded automatically when we loaded the original data. Here is how the processed peak
# annotation table looks like:
print(atac.uns['atac']['peak_annotation'].tail())

# Now we can plot average cut values in peaks corresponding to genes just by providing a gene name.
# By default, values in atac.raw are used for plotting.
ac.pl.pca(atac, color=["BCL11B", "CCR6", "KLF4"], average="total")
# We can also average peaks of each type separately:
ac.pl.pca(atac, color="BCL11B", average="peak_type")
# We see how this component space here resembles the one based on gene expression from the previous notebook.
# Looking at top loadings of first two components, we see how peaks linked to BCL11B (ENSG00000127152) and
# KLF4 (ENSG00000136826) demarcate lympohoid / myeloid axis while peaks linked to CCR6 (ENSG00000112486)
# define B cell axis.
#
# Now we will compute a neighbourhood graph for cells that we’ll use for clustering later on.
sc.pp.neighbors(atac, n_neighbors=10, n_pcs=30)