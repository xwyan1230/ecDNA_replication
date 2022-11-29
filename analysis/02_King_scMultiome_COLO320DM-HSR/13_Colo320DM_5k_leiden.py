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
print(rna)

data = pd.read_csv('%s_average-expression_enrich_in_C2_log2FC_low.txt' % prefix, na_values=['.'], sep='\t')
# mu.pp.filter_obs(rna, "leiden", lambda x: x.isin(["2", "4", "5"]))
# sorted = '245'
# rna = rna[~(rna.obs.leiden.isin(["4"]) & (rna.obs.phase.isin(['G2M'])))].copy()
# print(rna.var[rna.var.index == 'BX284613.2']['gene_ids'])

# sc.pl.umap(rna, color='MYC', legend_loc="on data")
for i in data['gene'].tolist()[:50]:
     sc.pl.umap(rna, color=i, legend_loc="on data", save='_%s_rna_%s' % (prefix, i))
     # sc.pl.umap(rna, color=i, legend_loc="on data")

print("DONE!")