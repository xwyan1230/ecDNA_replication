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

# RNA
rna = mdata.mod['rna']
rna.var['mt'] = rna.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.normalize_total(rna, target_sum=1e4)

rna_copy = rna.copy()
df_rna = pd.DataFrame.sparse.from_spmatrix(rna.X)
df_rna.columns = rna.var.index
df_rna.index = rna.obs.index
print(df_rna)

columns = df_rna.columns
stdev = []
for i in columns:
    stdev.append(np.std(df_rna[i].tolist()))

df_cc = pd.DataFrame()
df_cc['gene'] = rna.var.index
df_cc['all'] = stdev

df_cc.to_csv('%s%s_normalized_expression_stdev.txt' % (output_dir, prefix), index=False, sep='\t')

print("Done!")