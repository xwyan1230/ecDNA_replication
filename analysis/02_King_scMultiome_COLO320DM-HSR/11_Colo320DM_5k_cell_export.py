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
mu.pp.filter_obs(rna, "leiden", lambda x: x.isin(["2", "4", "5"]))
sorted = '245'
rna = rna[~(rna.obs.leiden.isin(["4"]) & (rna.obs.phase.isin(['G2M'])))].copy()
rna_sorted = rna.copy()

for i in ['G1', 'S', 'G2M']:
    df = pd.DataFrame()
    rna = rna_sorted.copy()
    mu.pp.filter_obs(rna, "phase", lambda x: x.isin([i]))
    df[i] = rna.obs.index.tolist()
    df.to_csv('%s%s_sorted_245_%s.txt' % (output_dir, prefix, i), index=False, header=False, sep='\t')

print("DONE!")



