import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/' \
           '02_cell_ranger_output_full/COLO320DM_5k_rep1_hg38/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

# LOAD DATA
mdata = mu.read_10x_h5(os.path.join(data_dir, "filtered_feature_bc_matrix.h5"))
mdata.var_names_make_unique()
print(mdata)

# RNA
rna = mdata.mod['rna']

# PREPROCESSING
# QC
rna.var['mt'] = rna.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(rna, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

# filter genes which expression is not detected
mu.pp.filter_var(rna, 'n_cells_by_counts', lambda x: x >= 3)
# filter cells
mu.pp.filter_obs(rna, 'n_genes_by_counts', lambda x: (x >= 200) & (x < 20000))
mu.pp.filter_obs(rna, 'total_counts', lambda x: x < 80000)
mu.pp.filter_obs(rna, 'pct_counts_mt', lambda x: x < 20)
sc.pl.violin(rna, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

# normalization
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)

# feature selection
sc.pp.highly_variable_genes(rna, min_mean=0.02, max_mean=4, min_disp=0.5)
sc.pl.highly_variable_genes(rna)
print(np.sum(rna.var.highly_variable))

# scaling
rna.raw = rna
sc.pp.scale(rna, max_value=10)
print(mdata)

# ANALYSIS
# following https://nbviewer.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
s_genes = [x.strip() for x in open('/Users/xwyan/Dropbox/LAB/Seqbasic/resources/'
                                   'regev_lab_cell_cycle_gene/regev_lab_cell_cycle_G1-S.txt')]
g2m_genes = [x.strip() for x in open('/Users/xwyan/Dropbox/LAB/Seqbasic/resources/'
                                     'regev_lab_cell_cycle_gene/regev_lab_cell_cycle_G2-M.txt')]
cell_cycle_genes = s_genes + g2m_genes
cell_cycle_genes = [x for x in cell_cycle_genes if x in rna.var_names]
print(len(cell_cycle_genes)) # this number can not be too small
# There are major differences in the way Scanpy and Seurat manage data, in particular we need to filter out cell
# cycle genes that are not present in our dataset to avoid errors.

sc.tl.score_genes_cell_cycle(rna, s_genes=s_genes, g2m_genes=g2m_genes)

# Here comes another difference from Seurat. The R package stores raw data, scaled data and variable genes
# information in separate slots, Scanpy instead keeps only one snapshot of the data. This implies that PCA is
# always calculated on the entire dataset. In order to calculate PCA reduction using only a subset of genes
# (like cell_cycle_genes), a trick should be used. Basically we create a dummy object to store information of
# PCA projection, which is then reincorporated into original dataset.
rna_cc_genes = rna[:, cell_cycle_genes]
sc.tl.pca(rna_cc_genes)
sc.pl.pca_scatter(rna_cc_genes, color='phase')
print(rna.obs['phase'])
print(mdata)

# In a typical use, you should regress out cell cycle effect following the below code. However, for current purpose,
# we want to separate cells based on their different cell cycle phase. Thus, do not regress out.
"""# As in the original vignette, cells can be easily separated by their cell cycle status when cell cycle genes
# are used. Now we can regress out both S score and G2M score.
sc.pp.regress_out(rna, ['S_score', 'G2M_score'])
sc.pp.scale(rna)

# Finally, we reproject dataset using cell cycle genes again. Since we regressed the scores, no effect of cell
# cycle is now evident.
rna_cc_genes = rna[:, cell_cycle_genes]
sc.tl.pca(rna_cc_genes)
sc.pl.pca_scatter(rna_cc_genes, color='phase')"""

# PCA and neighbourhood graph
sc.tl.pca(rna, svd_solver='arpack')
sc.pl.pca(rna, color='phase')
# Now we can compute a neighbourhood graph for cells
sc.pp.neighbors(rna, n_neighbors=10, n_pcs=20)

# ## Non-linear dimensionality reduction and clustering
# With the neighbourhood graph computed, we can now perform clustering. We will use leiden clustering as an example
sc.tl.leiden(rna, resolution=.5)
# To visualise the results, we’ll first generate a 2D latent space with cells that we can colour according to their
# cluster assignment.
sc.tl.umap(rna, spread=1., min_dist=.5, random_state=11)
sc.pl.umap(rna, color="leiden", legend_loc="on data")
sc.pl.umap(rna, color="phase", legend_loc="on data")

mdata.mod['atac'].uns['files'] = dict(mdata.mod['atac'].uns['files'])
mdata.mod['atac'].uns['atac'] = dict(mdata.mod['atac'].uns['atac'])
mdata.write(os.path.join(output_dir, "%s.h5mu" % prefix))