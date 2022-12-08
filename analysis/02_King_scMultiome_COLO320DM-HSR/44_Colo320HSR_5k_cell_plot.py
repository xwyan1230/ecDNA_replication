import scanpy as sc
import pandas as pd
import muon as mu
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# This is the directory where those files are downloaded to
data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320HSR_5k/03_analysis/'
output_dir = data_dir
os.chdir(data_dir)

rep = 'rep1'
prefix = 'COLO320HSR_5K_%s' % rep

mdata = mu.read("%s.h5mu" % prefix)
rna = mdata.mod['rna']

data_counts = pd.read_csv('%s_window1e6_sliding2e5_chr8_counts.txt' % prefix, na_values=['.'], sep='\t')
data_log2FC = pd.read_csv('%s_window1e6_sliding2e5_chr8_log2FC.txt' % prefix, na_values=['.'], sep='\t')
data_counts['cell'] = [x.split('#')[1] for x in data_counts.index]
data_log2FC['cell'] = [x.split('#')[1] for x in data_log2FC.index]

rna.obs['counts'] = [data_counts[data_counts['cell'] == rna.obs.index[x]]['X127600001'].tolist()[0]
                     if rna.obs.index[x] in data_counts['cell'].tolist() else 0 for x in range(len(rna.obs.index))]
rna.obs['log2FC'] = [data_log2FC[data_log2FC['cell'] == rna.obs.index[x]]['X127600001'].tolist()[0]
                     if rna.obs.index[x] in data_log2FC['cell'].tolist() else 0 for x in range(len(rna.obs.index))]

mu.pp.filter_obs(rna, 'log2FC', lambda x: x > 0)
print(rna)

C1 = pd.read_csv('%s%s_sorted_01267.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C1_high = pd.read_csv('%s%s_sorted_01267_MYC_high.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C1_low = pd.read_csv('%s%s_sorted_01267_MYC_low.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()

rna = rna[rna.obs.index.isin(C1)]
celllist = [C1_high, C1_low]
cellname = ['C1_MYC_high', 'C1_MYC_low']

# sc.pl.umap(rna, color='log2FC', legend_loc="on data", save='_%s_C2_log2FC' % prefix)

for i in range(len(celllist)):
    rna.obs['cell'] = [1 if rna.obs.index[x] in celllist[i] else 0 for x in range(len(rna.obs.index))]
    sc.pl.umap(rna, color='cell', legend_loc="on data", save='_%s_cell_%s' % (prefix, cellname[i]))

print("DONE!")