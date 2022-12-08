import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
os.chdir(data_dir)
rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

df = pd.read_csv('%s%s_normalized_average_expression_C2_log2FC_3group.txt' % (data_dir, prefix), na_values=['.'], sep='\t')
df_gene = pd.read_csv('%s_normalized_average_expression_enrich_in_C2_log2FC_high_all.txt' % prefix, na_values=['.'], sep='\t')

df_heatmap = pd.DataFrame(columns=['gene', 'high', 'mid', 'low'])
for i in df_gene['gene'].tolist()[:50]:
    df_heatmap.loc[len(df_heatmap.index)] = [i, df[df['gene'] == i]['C2_high'].tolist()[0],
                                             df[df['gene'] == i]['C2_mid'].tolist()[0],
                                             df[df['gene'] == i]['C2_low'].tolist()[0]]

normalized_by = 'high'
df_heatmap.index = df_heatmap['gene']
df_heatmap = df_heatmap.sort_values(by=[normalized_by], ascending=False)
df_drop = df_heatmap.drop(['gene'], axis=1)

df_normalized = pd.DataFrame()
df_normalized['high'] = df_drop['high']/df_drop[normalized_by]
df_normalized['mid'] = df_drop['mid']/df_drop[normalized_by]
df_normalized['low'] = df_drop['low']/df_drop[normalized_by]
df_normalized = df_normalized.sort_values(by=['mid'], ascending=False)

plt.subplots(figsize=(9, 12))
sns.heatmap(df_drop, yticklabels=1)
plt.savefig('%sfigures/heatmap_log2FC_high_all_top50_3group1.pdf' % output_dir)
plt.show()

plt.subplots(figsize=(9, 12))
sns.heatmap(df_drop, yticklabels=1, vmax=20)
plt.savefig('%sfigures/heatmap_log2FC_high_all_top50_3group1_vmax20.pdf' % output_dir)
plt.show()

plt.subplots(figsize=(9, 12))
sns.heatmap(df_normalized, yticklabels=1)
plt.savefig('%sfigures/heatmap_log2FC_high_all_top50_3group1_normalized.pdf' % output_dir)
plt.show()