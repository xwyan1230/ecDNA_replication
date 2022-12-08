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

genelist = ['Notch_ligand_and_receptor']

df_express = pd.DataFrame()
for i in genelist:
    df = pd.read_csv('%s%s_expression_level.txt' % (data_dir, i), na_values=['.'], sep='\t')
    df_sub = df[df['DM_all'] > 0.2].copy()
    df_express = pd.concat([df_express, df_sub], axis=0)

df_express = df_express.reset_index()
df_express['high/low'] = df_express['DM_C2_log2FC_high']/(df_express['DM_C2_log2FC_high']+df_express['DM_C2_log2FC_low'])

plt.subplots(figsize=(9, 6))
sns.barplot(data=df_express, y='gene', x='high/low')
plt.axvline(x=0.5, linestyle='--', c='r')
plt.xlim([0.3, 0.7])
plt.savefig('%sfigures/log2FC_ratio_Notch_ligand_and_receptor.pdf' % output_dir)
plt.close()

"""df_sns = pd.DataFrame()
df_sns['gene'] = df['gene'].tolist() + df['gene'].tolist() + df['gene'].tolist() + df['gene'].tolist() + \
                 df['gene'].tolist() + df['gene'].tolist()
df_sns['expression'] = df['HSR_all'].tolist() + df['DM_all'].tolist() + df['DM_C1'].tolist() + df['DM_C2'].tolist() + \
                       df['DM_C2_log2FC_high'].tolist() + df['DM_C2_log2FC_low'].tolist()
df_sns['expression_group'] = ['HSR_all'] * len(df) + ['DM_all'] * len(df) + ['DM_C1'] * len(df) + ['DM_C2'] * len(df) + \
                             ['DM_C2_log2FC_high'] * len(df) + ['DM_C2_log2FC_low'] * len(df)

plt.subplots(figsize=(len(df), 6))
sns.barplot(data=df_sns, x="gene", y="expression", hue='expression_group')
plt.savefig('%sfigures/%s_expression_level.pdf' % (output_dir, genelist))
plt.close()"""