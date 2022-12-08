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

# data = pd.read_csv('COLO320DM_5K_rep1_normalized_average_expression_enrich_in_C2_log2FC_high_all.txt', na_values=['.'], sep='\t')
data = pd.read_csv('expression_level_Wnt_full.txt', na_values=['.'], sep='\t')
data = data[data['HSR_all'] > 0.5].copy().reset_index()

# i = 'KLK6'
sc.pl.umap(rna, color='phase', legend_loc="on data", save='_phase')
#for i in data['gene'].tolist():
#    sc.pl.umap(rna, color=i, legend_loc="on data", save='_%s_rna_%s' % (prefix, i))
     # sc.pl.umap(rna, color=i, legend_loc="on data")

print("DONE!")