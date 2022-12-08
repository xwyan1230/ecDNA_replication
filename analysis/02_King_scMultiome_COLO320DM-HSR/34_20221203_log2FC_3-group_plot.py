import scanpy as sc
import pandas as pd
import muon as mu
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# This is the directory where those files are downloaded to
data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = data_dir
os.chdir(data_dir)

rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

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

C2 = pd.read_csv('%s%s_sorted_245.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
# C2_log2FC_high = pd.read_csv('%s%s_sorted_245_log2FC_high.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
# C2_log2FC_low = pd.read_csv('%s%s_sorted_245_log2FC_low.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C2_log2FC_3group1_high = pd.read_csv('%s%s_sorted_245_log2FC_3group1_high.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C2_log2FC_3group1_low = pd.read_csv('%s%s_sorted_245_log2FC_3group1_low.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
C2_log2FC_3group1_mid = pd.read_csv('%s%s_sorted_245_log2FC_3group1_mid.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()

rna = rna[rna.obs.index.isin(C2)]
celllist = [C2_log2FC_3group1_low, C2_log2FC_3group1_high, C2_log2FC_3group1_mid]
cellname = ['C2_log2FC_3group1_low', 'C2_log2FC_3group1_high', 'C2_log2FC_3group1_mid']

i = 'MYC'
sc.pl.umap(rna, color=i, legend_loc="on data", save='_%s_C2_%s' % (prefix, i))

"""for i in range(len(celllist)):
    rna.obs['cell'] = [1 if rna.obs.index[x] in celllist[i] else 0 for x in range(len(rna.obs.index))]
    sc.pl.umap(rna, color='cell', legend_loc="on data", save='_%s_cell_%s' % (prefix, cellname[i]))"""

print("DONE!")