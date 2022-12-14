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

mu.pp.filter_obs(rna, "leiden", lambda x: x.isin(["2", "4", "5"]))
rna = rna[~(rna.obs.leiden.isin(["4"]) & (rna.obs.phase.isin(['G2M'])))].copy()
# print(rna.var[rna.var.index == 'BX284613.2']['gene_ids'])

data = pd.read_csv('%s_normalized_average_expression_enrich_in_C2_log2FC_high_all.txt' % prefix, na_values=['.'], sep='\t')

# i = 'MYC'
# sc.pl.umap(rna, color=i, legend_loc="on data")

genelist = ['GLI3', 'PTCH1', 'CAMK2D', 'SMAD3', 'SMAD4', 'YAP1']
for i in genelist:
       sc.pl.umap(rna, color=i, legend_loc="on data", save='_%s_C2_%s' % (prefix, i))

# for i in data['gene'].tolist()[:50]:
#       sc.pl.umap(rna, color=i, legend_loc="on data", save='_%s_rna_%s' % (prefix, i))
     # sc.pl.umap(rna, color=i, legend_loc="on data")

print("DONE!")