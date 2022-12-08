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

df = pd.read_csv('%s%s_normalized_average_expression_enrich_in_C2_log2FC_high_all.txt' % (data_dir, prefix), na_values=['.'], sep='\t')
df_bg = pd.read_csv('%s%s_normalized_average_expression.txt' % (data_dir, prefix), na_values=['.'], sep='\t')
print(len(df))

df['rank'] = df.index
df['C2/C1+C2'] = df['C2']/(df['C1']+df['C2'])
df['gamma_C2/C1'] = np.abs(df['C2/C1+C2']-0.5)*df['p_C2/C1']
df['high/low'] = df['C2_log2FC_high']/(df['C2_log2FC_high'] + df['C2_log2FC_low'])
df_bg['C2/C1+C2'] = df_bg['C2']/(df_bg['C1']+df_bg['C2'])
df_bg['gamma_C2/C1'] = np.abs(df_bg['C2/C1+C2']-0.5)*df_bg['p_C2/C1']
df_bg['high/low'] = df_bg['C2_log2FC_high']/(df_bg['C2_log2FC_high'] + df_bg['C2_log2FC_low'])

# global copy number
df_G1 = df[(df['C2/C1+C2'] > 0.4) & (df['C2/C1+C2'] < 0.6)].copy()
# df_G1.to_csv('%s%s_normalized_average_expression_enrich_in_C2_log2FC_low_global.txt' % (output_dir, prefix), index=False, sep='\t')
print(len(df_G1))
# C2_local
df_G2 = df[df['C2/C1+C2'] >= 0.7].copy()
# df_G2.to_csv('%s%s_normalized_average_expression_enrich_in_C2_log2FC_low_local_C2_cutoff7.txt' % (output_dir, prefix), index=False, sep='\t')
print(len(df_G2))
# C1_local
df_G3 = df[df['C2/C1+C2'] <= 0.3].copy()
# df_G3.to_csv('%s%s_normalized_average_expression_enrich_in_C2_log2FC_low_local_C1_cutoff3.txt' % (output_dir, prefix), index=False, sep='\t')
print(len(df_G3))

df_sub = df_bg[df_bg['p_C2/C1'] > -np.log(0.05)].copy()
df_sub1 = df_bg[df_bg['p_log2FC_high/low'] > -np.log(0.05)].copy()

# 2d volcano plot
fig, ax = plt.subplots()
plt.scatter(df_sub['C2/C1+C2'], df_sub['p_C2/C1'], c='black', s=5, alpha=0.5)
plt.scatter(df_G1['C2/C1+C2'][:50], df_G1['p_C2/C1'][:50], c='r', s=5, alpha=1)
plt.ylim([0, 700])
plt.xlabel('C2/C1+C2')
plt.ylabel('-log(p)')
plt.savefig('%sfigures/%s_expression-level_log2FC_high_global_in_C2-vs-C1.pdf' % (output_dir, prefix))
plt.close()

fig, ax = plt.subplots()
plt.scatter(df_sub1['high/low'], df_sub1['p_log2FC_high/low'], c='black', s=5, alpha=0.5)
plt.scatter(df_G1['high/low'][:50], df_G1['p_log2FC_high/low'][:50], c='r', s=5, alpha=1)
plt.xlabel('log2FC_high/low')
plt.ylabel('-log(p)')
plt.savefig('%sfigures/%s_expression-level_log2FC_high_global_in_log2FC_high_vs_low.pdf' % (output_dir, prefix))
plt.close()

fig, ax = plt.subplots()
ax.scatter(df_sub['C2/C1+C2'], df_sub['p_C2/C1'], c='black', s=5, alpha=0.5)
ax.scatter(df_G2['C2/C1+C2'][:20], df_G2['p_C2/C1'][:20], c='r', s=5, alpha=1)
for i in range(20):
    ax.annotate(df_G2['gene'].tolist()[i], (df_G2['C2/C1+C2'].tolist()[i], df_G2['p_C2/C1'].tolist()[i]))
plt.ylim([0, 700])
plt.xlabel('C2/C1+C2')
plt.ylabel('-log(p)')
plt.savefig('%sfigures/%s_expression-level_log2FC_high_local_C2_cutoff7_in_C2-vs-C1.pdf' % (output_dir, prefix))
plt.close()

fig, ax = plt.subplots()
ax.scatter(df_sub1['high/low'], df_sub1['p_log2FC_high/low'], c='black', s=5, alpha=0.5)
ax.scatter(df_G2['high/low'][:20], df_G2['p_log2FC_high/low'][:20], c='r', s=5, alpha=1)
for i in range(20):
    ax.annotate(df_G2['gene'].tolist()[i], (df_G2['high/low'].tolist()[i], df_G2['p_log2FC_high/low'].tolist()[i]))
plt.xlabel('log2FC_high/low')
plt.ylabel('-log(p)')
plt.savefig('%sfigures/%s_expression-level_log2FC_high_local_C2_in_log2FC_high_vs_low.pdf' % (output_dir, prefix))
plt.close()

fig, ax = plt.subplots()
plt.scatter(df_sub['C2/C1+C2'], df_sub['p_C2/C1'], c='black', s=5, alpha=0.5)
plt.scatter(df_G3['C2/C1+C2'][:20], df_G3['p_C2/C1'][:20], c='r', s=5, alpha=1)
for i in range(20):
    ax.annotate(df_G3['gene'].tolist()[i], (df_G3['C2/C1+C2'].tolist()[i], df_G3['p_C2/C1'].tolist()[i]))
plt.ylim([0, 700])
plt.xlabel('C2/C1+C2')
plt.ylabel('-log(p)')
plt.savefig('%sfigures/%s_expression-level_log2FC_high_local_C1_cutoff3_in_C2-vs-C1.pdf' % (output_dir, prefix))
plt.close()

fig, ax = plt.subplots()
ax.scatter(df_sub1['high/low'], df_sub1['p_log2FC_high/low'], c='black', s=5, alpha=0.5)
ax.scatter(df_G3['high/low'][:20], df_G3['p_log2FC_high/low'][:20], c='r', s=5, alpha=1)
for i in range(20):
    ax.annotate(df_G3['gene'].tolist()[i], (df_G3['high/low'].tolist()[i], df_G3['p_log2FC_high/low'].tolist()[i]))
plt.xlabel('log2FC_high/low')
plt.ylabel('-log(p)')
plt.savefig('%sfigures/%s_expression-level_log2FC_high_local_C1_in_log2FC_high_vs_low.pdf' % (output_dir, prefix))
plt.close()


