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

RNAseq = pd.read_csv('%sCOLO320DM_RNAseq_Rep1_FPKMandTPM.txt' % data_dir, na_values=['.'], sep='\t')
scRNAseq = pd.read_csv('%s%s_normalized_average_expression.txt' % (data_dir, prefix), na_values=['.'], sep='\t')

print(len(RNAseq))
print(len(scRNAseq))

df = pd.DataFrame()
df['gene'] = scRNAseq['gene']
df['scRNAseq_mean-counts'] = scRNAseq['all']
df['RNAseq_counts'] = [RNAseq[RNAseq['Geneid'] == i]['Counts'].tolist()[0] if i in RNAseq['Geneid'].tolist() else -1 for i in scRNAseq['gene'].tolist()]
df['RNAseq_FPM'] = [RNAseq[RNAseq['Geneid'] == i]['FPM'].tolist()[0] if i in RNAseq['Geneid'].tolist() else -1 for i in scRNAseq['gene'].tolist()]
df['RNAseq_FPKM'] = [RNAseq[RNAseq['Geneid'] == i]['FPKM'].tolist()[0] if i in RNAseq['Geneid'].tolist() else -1 for i in scRNAseq['gene'].tolist()]
df['RNAseq_FPK'] = [RNAseq[RNAseq['Geneid'] == i]['FPK'].tolist()[0] if i in RNAseq['Geneid'].tolist() else -1 for i in scRNAseq['gene'].tolist()]
df['RNAseq_TPM'] = [RNAseq[RNAseq['Geneid'] == i]['TPM'].tolist()[0] if i in RNAseq['Geneid'].tolist() else -1 for i in scRNAseq['gene'].tolist()]

df.to_csv('%sCOLO320DM_scRNAseq_vs_RNAseq.txt' % output_dir, index=False, sep='\t')

