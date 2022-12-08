import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/01_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/01_analysis/'
os.chdir(data_dir)

df_HSR = pd.read_csv('%sCOLO320HSR_5K_rep1_normalized_average_expression.txt' % data_dir, na_values=['.'], sep='\t')
df_DM = pd.read_csv('%sCOLO320DM_5K_rep1_normalized_average_expression.txt' % data_dir, na_values=['.'], sep='\t')
data_low = pd.read_csv('COLO320DM_5K_rep1_normalized_average_expression_enrich_in_C2_log2FC_low_all.txt', na_values=['.'], sep='\t')
data_high_exclude = pd.read_csv('COLO320DM_5K_rep1_normalized_average_expression_enrich_in_C2_log2FC_high_all_exclude_ecDNA.txt', na_values=['.'], sep='\t')
data_high = pd.read_csv('COLO320DM_5K_rep1_normalized_average_expression_enrich_in_C2_log2FC_high_all.txt', na_values=['.'], sep='\t')
data_MYC_high = pd.read_csv('COLO320HSR_5K_rep1_normalized_average_expression_enrich_in_C1_MYC_high_all.txt', na_values=['.'], sep='\t')
data_MYC_low = pd.read_csv('COLO320HSR_5K_rep1_normalized_average_expression_enrich_in_C1_MYC_low_all.txt', na_values=['.'], sep='\t')
print(len(df_HSR))
print(len(df_DM))

df = pd.DataFrame({'gene': df_DM['gene'],
                   'log2FC_high_vs_low': df_DM['C2_log2FC_high']/(df_DM['C2_log2FC_high']+df_DM['C2_log2FC_low']),
                   'MYC_high_vs_low': df_HSR['C1_MYC_high']/(df_HSR['C1_MYC_high']+df_HSR['C1_MYC_low']),
                   'log2FC_p': df_DM['p_log2FC_high/low'],
                   'MYC_p': df_HSR['p_C1_high/low']})

df = df[df['log2FC_p'] > -np.log(0.05)].copy().reset_index(drop=True)

df_low = df[df['gene'].isin(data_low['gene'].tolist()[:50])].copy().reset_index(drop=True)
df_high_exclude = df[df['gene'].isin(data_high_exclude['gene'].tolist()[:50])].copy().reset_index(drop=True)
df_high = df[df['gene'].isin(data_high['gene'].tolist()[:50])].copy().reset_index(drop=True)
df_label1 = df_low[df_low['MYC_high_vs_low'] < 0.35]
gene_label = ['CDX2', 'LRATD2']
df_label2 = df_high[df_high['gene'].isin(gene_label)]
df_label = pd.concat([df_label1, df_label2], axis=0).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(9, 9))
ax.scatter(df['log2FC_high_vs_low'], df['MYC_high_vs_low'], c='black', alpha=0.5)
ax.scatter(df_low['log2FC_high_vs_low'], df_low['MYC_high_vs_low'], c='b')
ax.scatter(df_high['log2FC_high_vs_low'], df_high['MYC_high_vs_low'], c='orange')
ax.scatter(df_high_exclude['log2FC_high_vs_low'], df_high_exclude['MYC_high_vs_low'], c='r')
for i in range(len(df_label)):
    ax.annotate(df_label['gene'].tolist()[i], (df_label['log2FC_high_vs_low'].tolist()[i], df_label['MYC_high_vs_low'].tolist()[i]))
plt.plot([0, 1], [0, 1], c='black', linestyle='--')
plt.xlabel('COLO320DM_log2FC_high_vs_low')
plt.ylabel('COLO320HSR_MYC_high_vs_low')
plt.savefig('%sfigures/COLO320DM_log2FC_vs_COLO320HSR_MYC.pdf' % output_dir)
plt.show()


"""df = df[df['MYC_p'] > -np.log(0.05)].copy().reset_index(drop=True)
print(len(df))

df_low = df[df['gene'].isin(data_MYC_low['gene'].tolist()[:50])].copy().reset_index(drop=True)
df_high = df[df['gene'].isin(data_MYC_high['gene'].tolist()[:50])].copy().reset_index(drop=True)

fig, ax = plt.subplots(figsize=(9, 9))
ax.scatter(df['log2FC_high_vs_low'], df['MYC_high_vs_low'], c='black', alpha=0.5)
ax.scatter(df_low['log2FC_high_vs_low'], df_low['MYC_high_vs_low'], c='b')
ax.scatter(df_high['log2FC_high_vs_low'], df_high['MYC_high_vs_low'], c='r')
# for i in range(len(df_label)):
#     ax.annotate(df_label['gene'].tolist()[i], (df_label['log2FC_high_vs_low'].tolist()[i], df_label['MYC_high_vs_low'].tolist()[i]))
plt.plot([0, 1], [0, 1], c='black', linestyle='--')
plt.xlabel('COLO320DM_log2FC_high_vs_low')
plt.ylabel('COLO320HSR_MYC_high_vs_low')
plt.savefig('%sfigures/COLO320DM_log2FC_vs_COLO320HSR_MYC_plotting_MYC.pdf' % output_dir)
plt.show()"""
