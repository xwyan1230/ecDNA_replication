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

genelist = ['Wnt_low_to_high']

df_list = pd.DataFrame()
for i in genelist:
    df = pd.read_csv('%sexpression_level_%s.txt' % (data_dir, i), na_values=['.'], sep='\t')
    df_list = pd.concat([df_list, df], axis=0)

df_list = df_list.reset_index()
# df_list = df_list[~df_list['gene'].isin(['TCF4'])].copy().reset_index(drop=True)
# df_list = df_list[(df_list['DM_all'] > 0.5) | (df_list['HSR_all'] > 0.5)].copy().reset_index()

"""df_sns = pd.DataFrame()
df_sns['HSR'] = df_list['HSR_all']
df_sns['DM'] = df_list['DM_all']

df_sns['C1'] = df_list['DM_C1']
df_sns['C2_high'] = df_list['DM_C2_high']
df_sns['C2_mid'] = df_list['DM_C2_mid']
df_sns['C2_low'] = df_list['DM_C2_low']
df_sns.index = df_list['gene']"""

df_sns = pd.DataFrame()
df_sns['HSR'] = np.log1p(df_list['HSR_all'])
df_sns['DM'] = np.log1p(df_list['DM_all'])

df_sns['C1'] = np.log1p(df_list['DM_C1'])
df_sns['C2_high'] = np.log1p(df_list['DM_C2_high'])
df_sns['C2_mid'] = np.log1p(df_list['DM_C2_mid'])
df_sns['C2_low'] = np.log1p(df_list['DM_C2_low'])
df_sns.index = df_list['gene']
df_sns = df_sns.sort_values(by='C2_low', ascending=False)

"""df_sns = pd.DataFrame()
df_sns['HSR'] = [df_list['HSR_all'][i]/max([df_list['DM_C2_high'][i], df_list['DM_C2_low'][i]]) for i in range(len(df_list))]
df_sns['DM'] = [df_list['DM_all'][i]/max([df_list['DM_C2_high'][i], df_list['DM_C2_low'][i]]) for i in range(len(df_list))]

df_sns['C1'] = [df_list['DM_C1'][i]/max([df_list['DM_C2_high'][i], df_list['DM_C2_low'][i]]) for i in range(len(df_list))]
df_sns['C2_high'] = [df_list['DM_C2_high'][i]/max([df_list['DM_C2_high'][i], df_list['DM_C2_low'][i]]) for i in range(len(df_list))]
df_sns['C2_mid'] = [df_list['DM_C2_mid'][i]/max([df_list['DM_C2_high'][i], df_list['DM_C2_low'][i]]) for i in range(len(df_list))]
df_sns['C2_low'] = [df_list['DM_C2_low'][i]/max([df_list['DM_C2_high'][i], df_list['DM_C2_low'][i]]) for i in range(len(df_list))]
df_sns.index = df_list['gene']"""

plt.subplots(figsize=(5, len(df_list)/2))
sns.heatmap(data=df_sns, square=True, vmin=0)
plt.savefig('%sfigures/heatmap_expression-level_HSR-vs-DM_Wnt_low_to_high_log1p.pdf' % output_dir)
plt.show()