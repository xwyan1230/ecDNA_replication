import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
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
print(mdata)
C2 = pd.read_csv('%s%s_sorted_245.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()

# RNA
rna = mdata.mod['rna']


# PREPROCESSING
# QC
rna.var['mt'] = rna.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
rna_copy = rna.copy()
rna_copy = rna_copy[rna_copy.obs.index.isin(C2)]
# sc.pp.normalize_total(rna, target_sum=1e4)
# print(rna.var['mean_counts'])
# df_rna = pd.DataFrame.sparse.from_spmatrix(rna.X)
# print(df_rna[10:20])
sc.pp.calculate_qc_metrics(rna_copy, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
print(rna_copy.var['mean_counts'])
df_rna = pd.DataFrame.sparse.from_spmatrix(rna_copy.X)
print(df_rna)
# sc.pp.normalize_total(rna, target_sum=1e4)

print(rna)