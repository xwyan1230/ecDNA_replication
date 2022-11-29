import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
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
df_rna = pd.DataFrame.sparse.from_spmatrix(rna.X)
df_rna.columns = rna.var.index
df_rna.index = rna.obs.index
# df_rna.to_csv('%s%s_rna.txt' % (output_dir, prefix), sep='\t')

C2 = pd.read_csv('%s%s_sorted_245.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C1 = pd.read_csv('%s%s_sorted_01.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C2_G1 = pd.read_csv('%s%s_sorted_245_G1.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C2_S = pd.read_csv('%s%s_sorted_245_S.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C2_G2M = pd.read_csv('%s%s_sorted_245_G2M.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C2_log2FC_high = pd.read_csv('%s%s_sorted_245_log2FC_high.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C2_log2FC_low = pd.read_csv('%s%s_sorted_245_log2FC_low.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()


def count_pct_non_zero(df: pd.DataFrame):
    out = []
    for column_name in df.columns:
        column = df[column_name]
        # Get the count of Zeros in column
        count = (column != 0).sum()
        pct_non_zero = count*100.0/len(df)
        out.append(pct_non_zero)
    return out


def minuslogp_calculation(df1, df2):
    count = min([len(df1), len(df2)])
    out = []
    for i in df1.columns:
        p = -np.log(ks_2samp(df1[i].sample(n=count).tolist(), df2[i].sample(n=count).tolist())[1])
        out.append(p)
    return out


df_cc = pd.DataFrame()
df_cc['gene'] = df_rna.columns
print("calculating pct...")
df_cc['C2'] = df_rna[df_rna.index.isin(C2)].mean(0).tolist()
df_cc['C2_pct'] = count_pct_non_zero(df_rna[df_rna.index.isin(C2)])
df_cc['C1'] = df_rna[df_rna.index.isin(C1)].mean(0).tolist()
df_cc['C1_pct'] = count_pct_non_zero(df_rna[df_rna.index.isin(C1)])
df_cc['C2_G1'] = df_rna[df_rna.index.isin(C2_G1)].mean(0).tolist()
df_cc['C2_G1_pct'] = count_pct_non_zero(df_rna[df_rna.index.isin(C2_G1)])
df_cc['C2_S'] = df_rna[df_rna.index.isin(C2_S)].mean(0).tolist()
df_cc['C2_S_pct'] = count_pct_non_zero(df_rna[df_rna.index.isin(C2_S)])
df_cc['C2_G2M'] = df_rna[df_rna.index.isin(C2_G2M)].mean(0).tolist()
df_cc['C2_G2M_pct'] = count_pct_non_zero(df_rna[df_rna.index.isin(C2_G2M)])
df_cc['C2_S-G2M'] = df_rna[df_rna.index.isin(C2_S+C2_G2M)].mean(0).tolist()
df_cc['C2_S-G2M_pct'] = count_pct_non_zero(df_rna[df_rna.index.isin(C2_S+C2_G2M)])
df_cc['C2_log2FC_high'] = df_rna[df_rna.index.isin(C2_log2FC_high)].mean(0).tolist()
df_cc['C2_log2FC_high_pct'] = count_pct_non_zero(df_rna[df_rna.index.isin(C2_log2FC_high)])
df_cc['C2_log2FC_low'] = df_rna[df_rna.index.isin(C2_log2FC_low)].mean(0).tolist()
df_cc['C2_log2FC_low_pct'] = count_pct_non_zero(df_rna[df_rna.index.isin(C2_log2FC_low)])
print("calculating p...")
df_cc['p_C2/C1'] = minuslogp_calculation(df_rna[df_rna.index.isin(C1)], df_rna[df_rna.index.isin(C2)])
df_cc['p_G1/S-G2M'] = minuslogp_calculation(df_rna[df_rna.index.isin(C2_G1)], df_rna[df_rna.index.isin(C2_S+C2_G2M)])
df_cc['p_log2FC_high/low'] = minuslogp_calculation(df_rna[df_rna.index.isin(C2_log2FC_high)], df_rna[df_rna.index.isin(C2_log2FC_low)])
df_cc.to_csv('%s%s_average-expression.txt' % (output_dir, prefix), index=False, sep='\t')

print("Done!")


