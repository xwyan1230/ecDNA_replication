import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
os.chdir(data_dir)
rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

df = pd.read_csv('%s%s_normalized_average_expression_enrich_in_C2_log2FC_low_global.txt' % (data_dir, prefix), na_values=['.'], sep='\t')

summary = pd.DataFrame()
summary['gene'] = df['gene']
summary['rank'] = df['rank']
summary['chr_num'] = df['chr_num']
summary['start'] = df['start']
summary['end'] = df['end']
summary['scRNAseq_mean-counts'] = df['scRNAseq_mean-counts']
summary['RNAseq_FPKM'] = df['RNAseq_FPKM']
summary['log2FC_high/log2FC_high+low'] = df['C2_log2FC_high']/(df['C2_log2FC_high'] + df['C2_log2FC_low'])
summary['C2/C1+C2'] = df['C2/C1+C2']
summary.to_csv('%s%s_normalized_average_expression_enrich_in_C2_log2FC_low_global_summary.txt' % (data_dir, prefix), index=False, sep='\t')