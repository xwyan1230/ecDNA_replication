import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# This is the directory where those files are downloaded to
data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = data_dir
os.chdir(data_dir)

rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

df = pd.read_csv('%s%s_normalized_average_expression.txt' % (data_dir, prefix), na_values=['.'], sep='\t')
df_hit = pd.read_csv('%s%s_normalized_average_expression_enrich_in_C2_log2FC_high_all_by_value.txt' % (data_dir, prefix), na_values=['.'], sep='\t')
print(len(df))

df['C2/C1+C2'] = df['C2']/(df['C1']+df['C2'])
df['gamma_C2/C1'] = np.abs(df['C2/C1+C2']-0.5)*df['p_C2/C1']
df_hit['C2/C1+C2'] = df_hit['C2']/(df_hit['C1']+df_hit['C2'])
df_hit['gamma_C2/C1'] = np.abs(df_hit['C2/C1+C2']-0.5)*df_hit['p_C2/C1']

df_sub = df[df['p_C2/C1'] > -np.log(0.05)].copy()

# 2d volcano plot
fig, ax = plt.subplots()
plt.scatter(df_sub['C2/C1+C2'], df_sub['p_C2/C1'], c='black', s=5, alpha=0.5)
plt.scatter(df_hit['C2/C1+C2'][:100], df_hit['p_C2/C1'][:100], c='r', s=5, alpha=1)
plt.ylim([0, 700])
plt.xlabel('C2/C1+C2')
plt.ylabel('-log(p)')
plt.savefig('%sfigures/%s_expression-level_log2FC_high_by_value_in_C2-vs-C1.pdf' % (output_dir, prefix))
plt.close()
