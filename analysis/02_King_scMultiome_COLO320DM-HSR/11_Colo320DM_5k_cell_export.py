import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
# Import a module with ATAC-seq-related functions
from muon import atac as ac
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

data_counts = pd.read_csv('%s_window1e6_sliding2e5_chr8_counts.txt' % prefix, na_values=['.'], sep='\t')
data_log2FC = pd.read_csv('%s_window1e6_sliding2e5_chr8_log2FC.txt' % prefix, na_values=['.'], sep='\t')
data_counts['cell'] = [x.split('#')[1] for x in data_counts.index]
data_log2FC['cell'] = [x.split('#')[1] for x in data_log2FC.index]

counts = []
rna.obs['counts'] = [data_counts[data_counts['cell'] == rna.obs.index[x]]['X127600001'].tolist()[0]
                     if rna.obs.index[x] in data_counts['cell'].tolist() else 0 for x in range(len(rna.obs.index))]
rna.obs['log2FC'] = [data_log2FC[data_log2FC['cell'] == rna.obs.index[x]]['X127600001'].tolist()[0]
                     if rna.obs.index[x] in data_log2FC['cell'].tolist() else 0 for x in range(len(rna.obs.index))]

rna_245 = rna.copy()
mu.pp.filter_obs(rna_245, "leiden", lambda x: x.isin(["2", "4", "5"]))
rna_245 = rna_245[~(rna_245.obs.leiden.isin(["4"]) & (rna_245.obs.phase.isin(['G2M'])))].copy()
df = pd.DataFrame()
df['C2_245'] = rna_245.obs.index.tolist()
df.to_csv('%s%s_sorted_245.txt' % (output_dir, prefix), index=False, header=False, sep='\t')
rna_01 = rna.copy()
mu.pp.filter_obs(rna_01, "leiden", lambda x: x.isin(["0", "1"]))
rna_01 = rna_01[~(rna_01.obs.leiden.isin(["0"]) & (rna_01.obs.phase.isin(['G2M'])))].copy()
df = pd.DataFrame()
df['C1_01'] = rna_01.obs.index.tolist()
df.to_csv('%s%s_sorted_01.txt' % (output_dir, prefix), index=False, header=False, sep='\t')

for i in ['G1', 'S', 'G2M']:
    df = pd.DataFrame()
    rna = rna_245.copy()
    mu.pp.filter_obs(rna, "phase", lambda x: x.isin([i]))
    df[i] = rna.obs.index.tolist()
    df.to_csv('%s%s_sorted_245_%s.txt' % (output_dir, prefix, i), index=False, header=False, sep='\t')

df = pd.DataFrame()
rna = rna_245.copy()
mu.pp.filter_obs(rna, "log2FC", lambda x: x > 4.5)
df['log2FC_high'] = rna.obs.index.tolist()
df.to_csv('%s%s_sorted_245_log2FC_high.txt' % (output_dir, prefix), index=False, header=False, sep='\t')

df = pd.DataFrame()
rna = rna_245.copy()
mu.pp.filter_obs(rna, "log2FC", lambda x: (x <= 4.5) & (x > 0))
df['log2FC_low'] = rna.obs.index.tolist()
df.to_csv('%s%s_sorted_245_log2FC_low.txt' % (output_dir, prefix), index=False, header=False, sep='\t')

df = pd.DataFrame()
rna = rna_245.copy()
mu.pp.filter_obs(rna, "log2FC", lambda x: x > 4.8)
df['log2FC_high'] = rna.obs.index.tolist()
df.to_csv('%s%s_sorted_245_log2FC_3group1_high.txt' % (output_dir, prefix), index=False, header=False, sep='\t')

df = pd.DataFrame()
rna = rna_245.copy()
mu.pp.filter_obs(rna, "log2FC", lambda x: (x < 4) & (x > 0))
df['log2FC_low'] = rna.obs.index.tolist()
df.to_csv('%s%s_sorted_245_log2FC_3group1_low.txt' % (output_dir, prefix), index=False, header=False, sep='\t')

df = pd.DataFrame()
rna = rna_245.copy()
mu.pp.filter_obs(rna, "log2FC", lambda x: (x >= 4) & (x <= 4.8))
df['log2FC_mid'] = rna.obs.index.tolist()
df.to_csv('%s%s_sorted_245_log2FC_3group1_mid.txt' % (output_dir, prefix), index=False, header=False, sep='\t')

print("DONE!")



