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

df = pd.read_csv('%s%s_average-expression_enrich_in_C2.txt' % (data_dir, prefix), na_values=['.'], sep='\t')
print(len(df))

df['S-G2M/G1+S-G2M'] = df['C2_S-G2M']/(df['C2_S-G2M']+df['C2_G1'])
df['gamma_S-G2M/G1'] = np.abs(df['S-G2M/G1+S-G2M']-0.5)*df['p_G1/S-G2M']

df_SG2Menrich = df[df['S-G2M/G1+S-G2M'] > 0.5].sort_values(by=['gamma_S-G2M/G1'], ascending=False)
print(len(df_SG2Menrich))
df_SG2Menrich.to_csv('%s%s_average-expression_enrich_in_C2_S-G2M.txt' % (output_dir, prefix), index=False, sep='\t')
df_G1enrich = df[df['S-G2M/G1+S-G2M'] < 0.5].sort_values(by=['gamma_S-G2M/G1'], ascending=False)
print(len(df_G1enrich))
df_G1enrich.to_csv('%s%s_average-expression_enrich_in_C2_G1.txt' % (output_dir, prefix), index=False, sep='\t')

# ishit = np.abs(df['C2/C1+C2']-0.5)*df['p_C2/C1'] >= 10
# print(len(df[ishit]))

# 2d volcano plot
fig, ax = plt.subplots()
plt.scatter(df['S-G2M/G1+S-G2M'], df['p_G1/S-G2M'], c='black', s=5, alpha=0.5)
# plt.scatter(df[ishit]['C2/C1+C2'], df[ishit]['p_C2/C1'], s=5, c='r')
# plt.plot(np.linspace(0, 1, 1000), np.abs(10 / np.linspace(-0.5, 0.5, 1000)), 'k--', lw=.5)
# for i in range(len(df_high)):
#      ax.annotate(df_high['gene'].tolist()[i], (df_high['S/G1'].tolist()[i], df_high['G2M/G1'].tolist()[i]))
# plt.axvline(x=cutoff, linestyle='--', c='r')
# plt.axhline(y=cutoff, linestyle='--', c='r')
# plt.ylim([0, 700])
plt.xlabel('S-G2M/G1+S-G2M')
plt.ylabel('-log(p)')
plt.show()
# plt.savefig('%sfigures/%s_expression-level_C2_S-G2M-vs-C2_G1.pdf' % (output_dir, prefix))
# plt.close()


# get genes that only express in S and G2M
# df_SG2M = df[(df['C2_G2M'] > 0) & (df['C2_S'] > 0) & (df['C2_G1'] == 0)].copy()
# df_SG2M.to_csv('%s%s_average-expression_S-G2M_only.txt' % (output_dir, prefix), index=False, sep='\t')
# print(len(df_SG2M))

# 3d plot
"""fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(df_subgenes['S/G1'], df_subgenes['G2M/G1'], df_subgenes['C2/C1'], s=5, alpha=0.2, c='black')
ax.scatter(df_high['S/G1'], df_high['G2M/G1'], df_high['C2/C1'], s=5, alpha=1, c='r')
ax.set_xlim([-50, 50])
ax.set_ylim([-50, 50])
ax.set_zlim([-50, 50])
ax.set_xlabel('S/G1')
ax.set_ylabel('G2M/G1')
ax.set_zlabel('C2/C1')
# plt.savefig('%sfigures/%s_sorted_245_S-G2M_cutoff-%s.pdf' % (output_dir, prefix, cutoff))
plt.show()"""

# sc.pl.heatmap(rna, genes[:50], groupby='phase')