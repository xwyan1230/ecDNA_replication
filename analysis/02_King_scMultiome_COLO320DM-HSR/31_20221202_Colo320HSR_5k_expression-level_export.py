import os
import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
from scipy.stats import ks_2samp
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

data_dir1 = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320HSR_5k/' \
           '02_cellranger_output_King/COLO320HSR_5k_rep1_hg38/'
data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320HSR_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320HSR_5k/03_analysis/'
rep = 'rep1'
prefix = 'COLO320HSR_5K_%s' % rep

# LOAD DATA
mdata = mu.read_10x_h5(os.path.join(data_dir1, "filtered_feature_bc_matrix.h5"))
mdata.var_names_make_unique()

C1 = pd.read_csv('%s%s_sorted_01267.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C1_MYC_high = pd.read_csv('%s%s_sorted_01267_MYC_high.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C1_MYC_low = pd.read_csv('%s%s_sorted_01267_MYC_low.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()

# RNA
rna = mdata.mod['rna']
rna.var['mt'] = rna.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.normalize_total(rna, target_sum=1e4)

rna_copy = rna.copy()
df_rna = pd.DataFrame.sparse.from_spmatrix(rna.X)
df_rna.columns = rna.var.index
df_rna.index = rna.obs.index
sc.pp.calculate_qc_metrics(rna_copy, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

df_cc = pd.DataFrame()
df_cc['gene'] = rna.var.index
df_cc['all'] = rna_copy.var['mean_counts'].tolist()

cell = ['C1', 'C1_MYC_high', 'C1_MYC_low']
celllist = [C1, C1_MYC_high, C1_MYC_low]

for i in range(len(cell)):
    rna_copy = rna.copy()
    rna_copy = rna_copy[rna_copy.obs.index.isin(celllist[i])]
    sc.pp.calculate_qc_metrics(rna_copy, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    df_cc[cell[i]] = rna_copy.var['mean_counts'].tolist()


def minuslogp_calculation(df1, df2):
    count = min([len(df1), len(df2)])
    out = []
    for i in df1.columns:
        p = -np.log(ks_2samp(df1[i].sample(n=count).tolist(), df2[i].sample(n=count).tolist())[1])
        out.append(p)
    return out


print("calculating p...")
df_cc['p_C1_high/low'] = minuslogp_calculation(df_rna[df_rna.index.isin(C1_MYC_high)], df_rna[df_rna.index.isin(C1_MYC_low)])

df_cc.to_csv('%s%s_normalized_average_expression.txt' % (output_dir, prefix), index=False, sep='\t')

print("Done!")