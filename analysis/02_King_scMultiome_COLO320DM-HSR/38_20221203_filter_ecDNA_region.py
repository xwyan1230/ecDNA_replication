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

df = pd.read_csv('%secDNA_region.txt' % data_dir, na_values=['.'], sep='\t')
df_gene = pd.read_csv('%s_normalized_average_expression_enrich_in_C2_log2FC_high_all.txt' % prefix, na_values=['.'], sep='\t')

start = df['start'].tolist()[1]  # Hg38
end = df['end'].tolist()[1]  # Hg38

df_sort = df_gene[~((df_gene['chr_num'] == 8) & (((df_gene['start'] >= start) & (df_gene['start'] <= end)) | ((df_gene['end'] >= start) & (df_gene['end'] <= end))))].copy()
df_sort.to_csv('%s%s_normalized_average_expression_enrich_in_C2_log2FC_high_all_exclude_ecDNA.txt' % (output_dir, prefix), index=False, sep='\t')