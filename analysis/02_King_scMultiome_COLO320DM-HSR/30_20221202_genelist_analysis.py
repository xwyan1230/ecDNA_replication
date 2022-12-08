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

genelist = 'Wnt_low_to_high'

df = pd.read_csv('%sgenelist_%s.txt' % (data_dir, genelist), na_values=['.'], sep='\t')
df_DM = pd.read_csv('%sCOLO320DM_5K_rep1_normalized_average_expression.txt' % data_dir, na_values=['.'], sep='\t')
df_HSR = pd.read_csv('%sCOLO320HSR_5K_rep1_normalized_average_expression.txt' % data_dir, na_values=['.'], sep='\t')
df_3group = pd.read_csv('%sCOLO320DM_5K_rep1_normalized_average_expression_C2_log2FC_3group.txt' % data_dir, na_values=['.'], sep='\t')

for i in df['gene'].tolist():
    try:
        a = df_HSR[df_HSR['gene'] == i]['all'].tolist()[0]
    except:
        print(i)
"""df_sns = pd.DataFrame()
df_sns['gene'] = df['gene'].tolist() + df['gene'].tolist() + df['gene'].tolist() + df['gene'].tolist() + \
                 df['gene'].tolist() + df['gene'].tolist()
df_sns['expression'] = df['HSR_all'].tolist() + df['DM_all'].tolist() + df['DM_C1'].tolist() + df['DM_C2'].tolist() + \
                       df['DM_C2_log2FC_high'].tolist() + df['DM_C2_log2FC_low'].tolist()
df_sns['expression_group'] = ['HSR_all'] * len(df) + ['DM_all'] * len(df) + ['DM_C1'] * len(df) + ['DM_C2'] * len(df) + \
                             ['DM_C2_log2FC_high'] * len(df) + ['DM_C2_log2FC_low'] * len(df)

plt.subplots(figsize=(len(df), 6))
sns.barplot(data=df_sns[20:50], x="gene", y="expression", hue='expression_group')
# plt.savefig('%sfigures/expression_level_%s.pdf' % (output_dir, genelist))"""


df['HSR_all'] = [df_HSR[df_HSR['gene'] == i]['all'].tolist()[0] for i in df['gene'].tolist()]
df['DM_all'] = [df_DM[df_DM['gene'] == i]['all'].tolist()[0] for i in df['gene'].tolist()]
df['DM_C1'] = [df_DM[df_DM['gene'] == i]['C1'].tolist()[0] for i in df['gene'].tolist()]
df['DM_C2'] = [df_DM[df_DM['gene'] == i]['C2'].tolist()[0] for i in df['gene'].tolist()]
df['DM_C2_log2FC_high'] = [df_DM[df_DM['gene'] == i]['C2_log2FC_high'].tolist()[0] for i in df['gene'].tolist()]
df['DM_C2_log2FC_low'] = [df_DM[df_DM['gene'] == i]['C2_log2FC_low'].tolist()[0] for i in df['gene'].tolist()]
df['DM_C2_high'] = [df_3group[df_3group['gene'] == i]['C2_high'].tolist()[0] for i in df['gene'].tolist()]
df['DM_C2_mid'] = [df_3group[df_3group['gene'] == i]['C2_mid'].tolist()[0] for i in df['gene'].tolist()]
df['DM_C2_low'] = [df_3group[df_3group['gene'] == i]['C2_low'].tolist()[0] for i in df['gene'].tolist()]
df.to_csv('%sexpression_level_%s.txt' % (data_dir, genelist), index=False, sep='\t')
plt.show()