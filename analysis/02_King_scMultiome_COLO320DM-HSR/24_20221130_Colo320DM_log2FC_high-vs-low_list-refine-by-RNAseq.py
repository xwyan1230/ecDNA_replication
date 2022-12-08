import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

df = pd.read_csv('%sCOLO320DM_scRNAseq_vs_RNAseq.txt' % data_dir, na_values=['.'], sep='\t')
print(len(df))
gene_log2FC_high = pd.read_csv('%s%s_normalized_average_expression_enrich_in_C2_log2FC_low_global.txt' % (data_dir, prefix), na_values=['.'], sep='\t')

gene_log2FC_high['scRNAseq_mean-counts'] = [df[df['gene'] == i]['scRNAseq_mean-counts'].tolist()[0] for i in gene_log2FC_high['gene'].tolist()]
gene_log2FC_high['RNAseq_FPKM'] = [df[df['gene'] == i]['RNAseq_FPKM'].tolist()[0] for i in gene_log2FC_high['gene'].tolist()]
gene_log2FC_high.to_csv('%s%s_normalized_average_expression_enrich_in_C2_log2FC_low_global.txt' % (output_dir, prefix), index=False, sep='\t')
