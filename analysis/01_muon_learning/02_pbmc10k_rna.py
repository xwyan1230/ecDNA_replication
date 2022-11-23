# scRNA-seq data processing
# https://muon-tutorials.readthedocs.io/en/latest/single-cell-rna-atac/pbmc10k/1-Gene-Expression-Processing.html

# dataset required
# https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k?
# Filtered feature barcode matrix (HDF5)
# ATAC peak annotations based on proximal genes (TSV)
# ATAC Per fragment information file (TSV.GZ)
# ATAC Per fragment information index (TSV.GZ index)

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import muon as mu
import mudata as md
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# This is the directory where those files are downloaded to
data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221102_scMultiome_learning/pbmc10k/02_cell_ranger_output/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221102_scMultiome_learning/pbmc10k/03_analysis/'

# Remove file prefixes if any
prefix = "pbmc_granulocyte_sorted_10k_"
for file in os.listdir(data_dir):
    if file.startswith(prefix):
        new_filename = file[len(prefix):]
        os.rename(os.path.join(data_dir, file), os.path.join(data_dir, new_filename))

# #### LOAD DATA
#  muon will look for default files like atac_peak_annotation.tsv and atac_fragments.tsv.gz in the same folder
#  and will load peak annotation table and remember the path to the fragments file if they exist
mdata = mu.read_10x_h5(os.path.join(data_dir, "filtered_feature_bc_matrix.h5"))

mdata.var_names_make_unique()
print(mdata)

# #### RNA
rna = mdata.mod['rna']
print(rna)

# ### PREPROCESSING
# ## QC
rna.var['mt'] = rna.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
print(rna)
sc.pl.violin(rna, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

# filter genes which expression is not detected
mu.pp.filter_var(rna, 'n_cells_by_counts', lambda x: x >= 3)
# This is analogous to
#   sc.pp.filter_genes(rna, min_cells=3)
# but does in-place filtering and avoids copying the object

# filter cells
mu.pp.filter_obs(rna, 'n_genes_by_counts', lambda x: (x >= 200) & (x < 5000))
# This is analogous to
#   sc.pp.filter_cells(rna, min_genes=200)
#   rna = rna[rna.obs.n_genes_by_counts < 5000, :]
# but does in-place filtering avoiding copying the object
mu.pp.filter_obs(rna, 'total_counts', lambda x: x < 15000)
mu.pp.filter_obs(rna, 'pct_counts_mt', lambda x: x < 20)

sc.pl.violin(rna, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

# ## normalization
# We’ll normalise the data so that we get log-normalised counts to work with.
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)

# ## feature selection
# We will label highly variable genes that we’ll use for downstream analysis.
sc.pp.highly_variable_genes(rna, min_mean=0.02, max_mean=4, min_disp=0.5)
sc.pl.highly_variable_genes(rna)
np.sum(rna.var.highly_variable)

# ## scaling
# We’ll save log-normalised counts in a .raw slot:
rna.raw = rna
sc.pp.scale(rna, max_value=10)

# ### ANALYSIS
# ## PCA and neighbourhood graph
sc.tl.pca(rna, svd_solver='arpack')
sc.pl.pca(rna, color=['CD2', 'CD79A', 'KLF4', 'IRF8'])
sc.pl.pca_variance_ratio(rna, log=True)
# Now we can compute a neighbourhood graph for cells
sc.pp.neighbors(rna, n_neighbors=10, n_pcs=20)

# ## Non-linear dimensionality reduction and clustering
# With the neighbourhood graph computed, we can now perform clustering. We will use leiden clustering as an example
sc.tl.leiden(rna, resolution=.5)
# To visualise the results, we’ll first generate a 2D latent space with cells that we can colour according to their
# cluster assignment.
sc.tl.umap(rna, spread=1., min_dist=.5, random_state=11)
sc.pl.umap(rna, color="leiden", legend_loc="on data")

# ## Marker genes and celltypes
sc.tl.rank_genes_groups(rna, 'leiden', method='t-test')
result = rna.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.set_option('display.max_columns', 50)
print(pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(10))
sc.pl.rank_genes_groups(rna, n_genes=20, sharey=False)

# Exploring the data we notice clusters 9 and 15 seem to be composed of cells bearing markers for different cell
# lineages so likely to be noise (e.g. doublets). Cluster 12 has higher ribosomal gene expression when compared
# to other clusters. Cluster 16 seem to consist of proliferating cells.
# We will remove cells from these clusters before assigning cell types names to clusters.

mu.pp.filter_obs(rna, "leiden", lambda x: ~x.isin(["9", "15", "12", "16"]))
# Analogous to
#   rna = rna[~rna.obs.leiden.isin(["9", "15", "12", "16"])]
# but doesn't copy the object

new_cluster_names = {
    "0": "CD4+ memory T", "3": "CD8+ naïve T", "2": "CD4+ naïve T",
    "5": "CD8+ activated T", "7": "NK", "13": "MAIT",
    "6": "memory B", "10": "naïve B",
    "4": "CD14 mono", "1": "intermediate mono", "8": "CD16 mono",
    "11": "mDC", "14": "pDC",
}

rna.obs['celltype'] = rna.obs.leiden.astype("str").values
rna.obs.celltype = rna.obs.celltype.astype("category")
rna.obs.celltype = rna.obs.celltype.cat.rename_categories(new_cluster_names)

# We will also re-order categories for the next plots:
rna.obs.celltype.cat.reorder_categories([
    'CD4+ naïve T', 'CD4+ memory T', 'MAIT',
    'CD8+ naïve T', 'CD8+ activated T', 'NK',
    'naïve B', 'memory B',
    'CD14 mono', 'intermediate mono', 'CD16 mono',
    'mDC', 'pDC'], inplace=True)

# and take colours from a palette
cmap = plt.get_cmap('rainbow')
colors = cmap(np.linspace(0, 1, len(rna.obs.celltype.cat.categories)))

rna.uns["celltype_colors"] = list(map(matplotlib.colors.to_hex, colors))
sc.pl.umap(rna, color="celltype", legend_loc="on data")

# Finally, we’ll visualise some marker genes across cell types.
marker_genes = ['IL7R', 'TRAC',
                'ITGB1', # CD29
                'SLC4A10',
                'CD8A', 'CD8B', 'CCL5',
                'GNLY', 'NKG7',
                'CD79A', 'MS4A1', 'IGHM', 'IGHD',
                'IL4R', 'TCL1A',
                'KLF4', 'LYZ', 'S100A8', 'ITGAM', # CD11b
                'CD14', 'FCGR3A', 'MS4A7',
                'CST3', 'CLEC10A', 'IRF8', 'TCF4']

sc.pl.dotplot(rna, marker_genes, groupby='celltype');

# We will now write mdata object to an .h5mu file. It will contain all the changes we’ve done to the RNA modality
# (mdata.mod['rna']) inside it.
print(mdata.mod['atac'].uns)

mdata.mod['atac'].uns['files'] = dict(mdata.mod['atac'].uns['files'])
mdata.mod['atac'].uns['atac'] = dict(mdata.mod['atac'].uns['atac'])
mdata.write(os.path.join(output_dir, "pbmc10k.h5mu"))