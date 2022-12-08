import os
import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
from scipy.stats import ks_2samp
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

data_dir1 = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/' \
           '02_cellranger_output_King/COLO320DM_5k_rep1_hg38/'
data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

# LOAD DATA
mdata = mu.read_10x_h5(os.path.join(data_dir1, "filtered_feature_bc_matrix.h5"))
mdata.var_names_make_unique()

C2_high = pd.read_csv('%s%s_sorted_245_log2FC_3group1_high.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C2_low = pd.read_csv('%s%s_sorted_245_log2FC_3group1_low.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C2_mid = pd.read_csv('%s%s_sorted_245_log2FC_3group1_mid.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()

# RNA
rna = mdata.mod['rna']
rna.var['mt'] = rna.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.normalize_total(rna, target_sum=1e4)

cell = ['C2_high', 'C2_mid', 'C2_low']
celllist = [C2_high, C2_mid, C2_low]

rna_copy = rna.copy()
df_rna = pd.DataFrame.sparse.from_spmatrix(rna.X)
df_rna.columns = rna.var.index
df_rna.index = rna.obs.index
sc.pp.calculate_qc_metrics(rna_copy, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

df_cc = pd.DataFrame()
df_cc['gene'] = rna.var.index
df_cc['all'] = rna_copy.var['mean_counts'].tolist()

for i in range(len(cell)):
    rna_copy = rna.copy()
    rna_copy = rna_copy[rna_copy.obs.index.isin(celllist[i])]
    sc.pp.calculate_qc_metrics(rna_copy, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    df_cc[cell[i]] = rna_copy.var['mean_counts'].tolist()

df_cc.to_csv('%s%s_normalized_average_expression_C2_log2FC_3group.txt' % (output_dir, prefix), index=False, sep='\t')

print("Done!")