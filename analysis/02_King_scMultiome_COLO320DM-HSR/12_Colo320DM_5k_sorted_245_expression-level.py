import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
# Import a module with ATAC-seq-related functions
from muon import atac as ac
import matplotlib.pyplot as plt
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
print(mdata)

rna = mdata.mod['rna']
mu.pp.filter_obs(rna, "leiden", lambda x: x.isin(["2", "4", "5"]))
sorted = '245'
rna = rna[~(rna.obs.leiden.isin(["4"]) & (rna.obs.phase.isin(['G2M'])))].copy()
rna_sorted = rna.copy()
df_genes = pd.DataFrame()
df_genes['gene'] = rna.var.index.tolist()

for i in ['G1', 'S', 'G2M']:
    rna = rna_sorted.copy()
    mu.pp.filter_obs(rna, "phase", lambda x: x.isin([i]))
    sc.pp.calculate_qc_metrics(rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    df_genes[i] = rna.var.mean_counts.tolist()

df_genes['S/G1'] = df_genes['S']/df_genes['G1']
df_genes['G2M/G1'] = df_genes['G2M']/df_genes['G1']
df_genes.to_csv('%s%s_sorted_245_S-G2M.txt' % (output_dir, prefix), index=False, sep='\t')

df_subgenes = df_genes[(df_genes['S'] > 0.05) & (df_genes['G2M'] > 0.05) & (df_genes['G1'] > 0.05)].copy()
print(len(df_subgenes))

cutoff = 4

df_high = df_subgenes[(df_subgenes['S/G1'] > cutoff) & (df_subgenes['G2M/G1'] > cutoff)]
df_high.to_csv('%s%s_sorted_245_S-G2M_high_cutoff-%s.txt' % (output_dir, prefix, cutoff), index=False, sep='\t')

fig, ax = plt.subplots()
ax.scatter(df_subgenes['S/G1'], df_subgenes['G2M/G1'], s=5, alpha=0.5, c='black')
ax.scatter(df_high['S/G1'], df_high['G2M/G1'], s=5, alpha=1, c='r')
for i in range(len(df_high)):
    ax.annotate(df_high['gene'].tolist()[i], (df_high['S/G1'].tolist()[i], df_high['G2M/G1'].tolist()[i]))
plt.axvline(x=cutoff, linestyle='--', c='r')
plt.axhline(y=cutoff, linestyle='--', c='r')
plt.xlabel('S/G1')
plt.ylabel('G2M/G1')
plt.savefig('%sfigures/%s_sorted_245_S-G2M_cutoff-%s.pdf' % (output_dir, prefix, cutoff))
plt.close()


# sc.pl.heatmap(rna, genes[:50], groupby='phase')