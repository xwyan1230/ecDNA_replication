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

data_dir1 = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320HSR_5k/' \
           '02_cellranger_output_King/COLO320HSR_5k_rep1_hg38/'
data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320HSR_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320HSR_5k/03_analysis/'
os.chdir(data_dir)
rep = 'rep1'
prefix = 'COLO320HSR_5K_%s' % rep

# LOAD DATA
mdata1 = mu.read_10x_h5(os.path.join(data_dir1, "filtered_feature_bc_matrix.h5"))
mdata1.var_names_make_unique()

# RNA
rna1 = mdata1.mod['rna']
rna1.var['mt'] = rna1.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.normalize_total(rna1, target_sum=1e4)
sc.pp.log1p(rna1)

df_rna = pd.DataFrame.sparse.from_spmatrix(rna1.X)
df_rna.columns = rna1.var.index
df_rna.index = rna1.obs.index
data_express = pd.DataFrame({'cell': df_rna.index, 'MYC': df_rna['MYC'], 'PVT1': df_rna['PVT1']})
print(data_express)

mdata = mu.read("%s.h5mu" % prefix)
rna = mdata.mod['rna']
print(rna)

rna.obs['gene_MYC'] = [data_express[data_express['cell'] == rna.obs.index[x]]['MYC'].tolist()[0]
                       if rna.obs.index[x] in data_express['cell'].tolist() else -1 for x in range(len(rna.obs.index))]
rna.obs['gene_PVT1'] = [data_express[data_express['cell'] == rna.obs.index[x]]['PVT1'].tolist()[0]
                        if rna.obs.index[x] in data_express['cell'].tolist() else -1 for x in range(len(rna.obs.index))]
print("1")

data_counts = pd.read_csv('%s_window1e6_sliding2e5_chr8_counts.txt' % prefix, na_values=['.'], sep='\t')
data_log2FC = pd.read_csv('%s_window1e6_sliding2e5_chr8_log2FC.txt' % prefix, na_values=['.'], sep='\t')
data_counts['cell'] = [x.split('#')[1] for x in data_counts.index]
data_log2FC['cell'] = [x.split('#')[1] for x in data_log2FC.index]

rna.obs['counts'] = [data_counts[data_counts['cell'] == rna.obs.index[x]]['X127600001'].tolist()[0]
                     if rna.obs.index[x] in data_counts['cell'].tolist() else 0 for x in range(len(rna.obs.index))]
rna.obs['log2FC'] = [data_log2FC[data_log2FC['cell'] == rna.obs.index[x]]['X127600001'].tolist()[0]
                     if rna.obs.index[x] in data_log2FC['cell'].tolist() else 0 for x in range(len(rna.obs.index))]
print("1")

rna_01267 = rna.copy()
mu.pp.filter_obs(rna_01267, "leiden", lambda x: x.isin(["0", "1", "2", "6", "7"]))
rna_01267 = rna_01267[~(rna_01267.obs.leiden.isin(["1", "2"]) & (rna_01267.obs.phase.isin(['G2M'])))].copy()

df = pd.DataFrame()
df['C1_01267'] = rna_01267.obs.index.tolist()
df.to_csv('%s%s_sorted_01267.txt' % (output_dir, prefix), index=False, header=False, sep='\t')

mu.pp.filter_obs(rna_01267, "log2FC", lambda x: x > 0)
sc.pl.umap(rna_01267, color='phase', legend_loc="on data", save='_01267_phase')
sc.pl.umap(rna_01267, color='log2FC', legend_loc="on data", save='_01267_chr8_X127600001_log2FC')
sc.pl.umap(rna_01267, color='MYC', legend_loc="on data", save='_01267_MYC')
sc.pl.umap(rna_01267, color='PVT1', legend_loc="on data", save='_01267_PVT1')

df = pd.DataFrame()
rna = rna_01267.copy()
mu.pp.filter_obs(rna, "gene_MYC", lambda x: x > 2.2)
df['MYC_high'] = rna.obs.index.tolist()
print(len(df))
df.to_csv('%s%s_sorted_01267_MYC_high.txt' % (output_dir, prefix), index=False, header=False, sep='\t')

df = pd.DataFrame()
rna = rna_01267.copy()
mu.pp.filter_obs(rna, "gene_MYC", lambda x: (x <= 2.2) & (x >= 0))
df['MYC_low'] = rna.obs.index.tolist()
print(len(df))
df.to_csv('%s%s_sorted_01267_MYC_low.txt' % (output_dir, prefix), index=False, header=False, sep='\t')

print("DONE!")



