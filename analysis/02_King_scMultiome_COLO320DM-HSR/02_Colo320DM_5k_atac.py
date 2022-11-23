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
data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = data_dir
os.chdir(data_dir)

rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

mdata = mu.read("%s.h5mu" % prefix)
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
print(atac.var_names)
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