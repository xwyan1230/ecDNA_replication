import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
os.chdir(data_dir)
rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

df1 = pd.read_csv('%s%s_average-expression_enrich_in_C2_log2FC_high_gene-location.txt' % (data_dir, prefix), na_values=['.'], sep='\t')
df2 = pd.read_csv('%s%s_average-expression_enrich_in_C2_log2FC_low_gene-location.txt' % (data_dir, prefix), na_values=['.'], sep='\t')
df_chr = pd.read_csv('%shg38_chromSizes.txt' % data_dir, na_values=['.'], sep='\t')

feature = ['chr_num', 'start', 'end']
for f in feature:
    df1[f] = [float(df1[f][i]) for i in range(len(df1))]
    df2[f] = [float(df2[f][i]) for i in range(len(df2))]

fig, ax = plt.subplots()
for i in range(len(df_chr)-1):
    ax.plot([0, float(df_chr['end'].tolist()[i])], [25-float(df_chr['chr_num'].tolist()[i]), 25-float(df_chr['chr_num'].tolist()[i])], c='black')

# ax.scatter(df1['start'][:100], 25-df1['chr_num'][:100], c='r', alpha=0.3)
# ax.scatter(df2['start'][:100], 25-df2['chr_num'][:100], c='b', alpha=0.3)
# df1_random = df1.sample(n=50).reset_index()
# df1_random = df1.copy()
# ax.scatter(df1_random['start'], 25-df1_random['chr_num'], c='r', alpha=0.03)
# df2_random = df2.sample(n=100).reset_index()
df2_random = df2.copy()
ax.scatter(df2_random['start'], 25-df2_random['chr_num'], c='b', alpha=0.03)
plt.savefig('%sfigures/%s_chr_location_C2_log2FC_low_background.pdf' % (output_dir, prefix))
plt.close()
