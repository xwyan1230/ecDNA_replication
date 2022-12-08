import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

df = pd.read_csv('%sCOLO320DM_scRNAseq_vs_RNAseq.txt' % data_dir, na_values=['.'], sep='\t')
df_background = df[(df['RNAseq_FPKM'] > 0.5) & (df['scRNAseq_mean-counts'] > 0)].copy()
# df_background = df[(df['RNAseq_FPKM'] > 10) | (df['RNAseq_TPM'] > 10)]
# gene_log2FC_high = pd.read_csv('%s%s_average-expression_enrich_in_C2_log2FC_high_gene-location.txt' % (data_dir, prefix), na_values=['.'], sep='\t')['gene'].tolist()
# df_log2FC_high = df_background[df_background['gene'].isin(gene_log2FC_high[:100])].copy()
print(len(df_background))
# print(len(df_log2FC_high))

fig, ax = plt.subplots()
ax.scatter(df_background['scRNAseq_mean-counts'], df_background['RNAseq_FPKM'], c='black', s=5, alpha=0.2)
# ax.scatter(df_log2FC_high['scRNAseq_mean-counts'], df_log2FC_high['RNAseq_FPKM'], c='r', s=5, alpha=0.5)
plt.xlim([0, 20])
plt.ylim([0, 1500])
plt.xlabel('scRNAseq_mean-counts')
plt.ylabel('RNAseq_FPKM')
plt.show()
# plt.savefig('%sfigures/%s_scRNAseq_vs_RNAseq_red-top100_log2FC_high.pdf' % (output_dir, prefix))

"""print(len(df_background[df_background['RNAseq_FPKM'] < 0.5]))  # cutoff
print(len(df_background[(df_background['RNAseq_FPKM'] < 10) & (df_background['RNAseq_FPKM'] >= 0.5)]))  # low
print(len(df_background[(df_background['RNAseq_FPKM'] < 1000) & (df_background['RNAseq_FPKM'] >= 10)]))  # medium
print(len(df_background[df_background['RNAseq_FPKM'] >= 1000]))  # high

print(len(df_log2FC_high[df_log2FC_high['RNAseq_FPKM'] < 0.5]))  # cutoff
print(len(df_log2FC_high[(df_log2FC_high['RNAseq_FPKM'] < 10) & (df_log2FC_high['RNAseq_FPKM'] >= 0.5)]))  # low
print(len(df_log2FC_high[(df_log2FC_high['RNAseq_FPKM'] < 1000) & (df_log2FC_high['RNAseq_FPKM'] >= 10)]))  # medium
print(len(df_log2FC_high[df_log2FC_high['RNAseq_FPKM'] >= 1000]))  # high"""