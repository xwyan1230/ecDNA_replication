import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
os.chdir(data_dir)
rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

df = pd.read_csv('%s%s_normalized_average_expression_enrich_in_C2_log2FC_high_all.txt' % (data_dir, prefix), na_values=['.'], sep='\t')
df_gene = pd.read_csv('%s%s_gene.txt' % (data_dir, prefix), na_values=['.'], sep='\t')
df['chr'] = [df_gene[df_gene['gene'] == i]['chr'].tolist()[0] for i in df['gene'].tolist()]
chr_num = []
for i in df['chr']:
    if type(i) == float:
        chr_num.append(np.nan)
    elif i.startswith('chr'):
        num_temp = i[3:]
        if num_temp == 'X':
            chr_num.append(23)
        elif num_temp == 'Y':
            chr_num.append(24)
        else:
            chr_num.append(float(num_temp))
    else:
        chr_num.append(np.nan)
df['chr_num'] = chr_num
df['start'] = [float(df_gene[df_gene['gene'] == i]['start'].tolist()[0]) for i in df['gene'].tolist()]
df['end'] = [float(df_gene[df_gene['gene'] == i]['end'].tolist()[0]) for i in df['gene'].tolist()]
df.to_csv('%s%s_normalized_average_expression_enrich_in_C2_log2FC_high_all.txt' % (data_dir, prefix), index=False, sep='\t')