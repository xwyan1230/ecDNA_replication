import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
os.chdir(data_dir)
rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

df = pd.read_csv('%s%s_normalized_average_expression_C2_log2FC_3group.txt' % (data_dir, prefix), na_values=['.'], sep='\t')
df_gene = pd.read_csv('%s_normalized_average_expression_enrich_in_C2_log2FC_low_all.txt' % prefix, na_values=['.'], sep='\t')

df_heatmap = pd.DataFrame(columns=['gene', 'high', 'mid', 'low'])
for i in df_gene['gene'].tolist()[:30]:
    df_heatmap.loc[len(df_heatmap.index)] = [i, df[df['gene'] == i]['C2_high'].tolist()[0],
                                             df[df['gene'] == i]['C2_mid'].tolist()[0],
                                             df[df['gene'] == i]['C2_low'].tolist()[0]]
df_heatmap.index = df_heatmap['gene']
df_heatmap = df_heatmap.sort_values(by=['low'], ascending=False)
genelist = df_heatmap['gene'].tolist()

df_sns = pd.DataFrame()
df_sns['gene'] = genelist + genelist + genelist
df_sns['expression'] = df_heatmap['high'].tolist() + df_heatmap['mid'].tolist() + df_heatmap['low'].tolist()
df_sns['expression_group'] = ['high'] * len(df_heatmap) + ['mid'] * len(df_heatmap) + ['low'] * len(df_heatmap)

plt.subplots(figsize=(len(df_heatmap), 6))
sns.barplot(data=df_sns, x="gene", y="expression", hue='expression_group')
# plt.savefig('%sfigures/expression_level_%s.pdf' % (output_dir, genelist))
plt.show()