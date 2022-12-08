import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/01_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/01_analysis/'
os.chdir(data_dir)

df_HSR = pd.read_csv('%sCOLO320HSR_5K_rep1_normalized_expression_stdev.txt' % data_dir, na_values=['.'], sep='\t')
df_DM = pd.read_csv('%sCOLO320DM_5K_rep1_normalized_expression_stdev.txt' % data_dir, na_values=['.'], sep='\t')
data_low = pd.read_csv('COLO320DM_5K_rep1_normalized_average_expression_enrich_in_C2_log2FC_low_all.txt', na_values=['.'], sep='\t')
data_high_exclude = pd.read_csv('COLO320DM_5K_rep1_normalized_average_expression_enrich_in_C2_log2FC_high_all_exclude_ecDNA.txt', na_values=['.'], sep='\t')
data_high = pd.read_csv('COLO320DM_5K_rep1_normalized_average_expression_enrich_in_C2_log2FC_high_all.txt', na_values=['.'], sep='\t')
gene_wnt = pd.read_csv('genelist_Wnt_high_to_low.txt', na_values=['.'], sep='\t')['gene'].tolist() + \
           pd.read_csv('genelist_Wnt_low_to_high.txt', na_values=['.'], sep='\t')['gene'].tolist()
gene_wnt_full = pd.read_csv('genelist_Wnt_full.txt', na_values=['.'], sep='\t')['gene'].tolist()
# data_MYC_high = pd.read_csv('COLO320HSR_5K_rep1_normalized_average_expression_enrich_in_C1_MYC_high_all.txt', na_values=['.'], sep='\t')
# data_MYC_low = pd.read_csv('COLO320HSR_5K_rep1_normalized_average_expression_enrich_in_C1_MYC_low_all.txt', na_values=['.'], sep='\t')

df = pd.DataFrame({'gene': df_DM['gene'],
                   'DM': np.log1p(df_DM['all']),
                   'HSR': np.log1p(df_HSR['all']),
                   'DM/HSR': np.log1p(df_DM['all'])/np.log1p(df_HSR['all'])})

df_low = df[df['gene'].isin(data_low['gene'].tolist()[:50])].copy().reset_index(drop=True)
df_high_exclude = df[df['gene'].isin(data_high_exclude['gene'].tolist()[:50])].copy().reset_index(drop=True)
df_high = df[df['gene'].isin(data_high['gene'].tolist()[:50])].copy().reset_index(drop=True)
df_wnt = df[df['gene'].isin(gene_wnt)].copy().reset_index(drop=True)
df_wnt_full = df[df['gene'].isin(gene_wnt_full)].copy().reset_index(drop=True)
df_label = df[((df['DM/HSR']>2) & (df['DM'] > 1)) | ((df['DM/HSR']<0.5) & (df['HSR'] > 1))]
# df_label = df_label[df_label['gene'].isin(['NOTUM', 'GPC5', 'NKD1', 'RNF43'])]

df_sort = df[(df['DM/HSR'] > 2) & (df['DM'] > 1)]
# df_sort.to_csv('%sstdev_hit_DMvsHSR_2fold_DM_over1.txt' % output_dir, index=False, sep='\t')
print(df_sort)

df_sort1 = df[(df['DM/HSR'] <0.5) & (df['HSR'] > 1)]
# df_sort1.to_csv('%sstdev_hit_HSRvsDM_2fold_HSR_over1.txt' % output_dir, index=False, sep='\t')

limit = 5
fig, ax = plt.subplots(figsize=(9, 9))
ax.scatter(df['HSR'], df['DM'], c='black', alpha=0.1)
# ax.scatter(df_wnt_full['HSR'], df_wnt_full['DM'], c='orange')
# plt.scatter(df_wnt['HSR'], df_wnt['DM'], c='orange')
# plt.scatter(df_low['HSR'], df_low['DM'], c='b')
# plt.scatter(df_high['HSR'], df_high['DM'], c='orange')
# plt.scatter(df_high_exclude['HSR'], df_high_exclude['DM'], c='r')
for i in range(len(df_label)):
     ax.annotate(df_label['gene'].tolist()[i], (df_label['HSR'].tolist()[i], df_label['DM'].tolist()[i]))
plt.plot([0, limit], [0, limit], c='black', linestyle='--')
plt.xlim([0, limit])
plt.ylim([0, limit])
plt.xlabel('COLO320HSR expression stdev (log1p)')
plt.ylabel('COLO320DM expression stdev (log1p)')
plt.savefig('%sfigures/scatter_expression_stdev_label.pdf' % output_dir)
plt.show()