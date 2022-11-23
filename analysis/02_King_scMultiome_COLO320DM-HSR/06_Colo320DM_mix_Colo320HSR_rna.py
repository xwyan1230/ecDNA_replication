import scanpy as sc
import pandas as pd
import anndata as ad
import muon as mu
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# This is the directory where those files are downloaded to
HSR_data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_scMultiome_ColoDM-ColoHSR/COLO320HSR_5k/03_analysis/'
DM_data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_scMultiome_ColoDM-ColoHSR/01_analysis/'
os.chdir(output_dir)

rep = 'rep1'
HSR_prefix = 'COLO320HSR_5K_%s' % rep
DM_prefix = 'COLO320DM_5K_%s' % rep

mdata_HSR = mu.read("%s%s.h5mu" % (HSR_data_dir, HSR_prefix))
mdata_DM = mu.read("%s%s.h5mu" % (DM_data_dir, DM_prefix))
rna_HSR = mdata_HSR.mod['rna']
rna_DM = mdata_DM.mod['rna']
rna_DM.obs['sample'] = ['DM'] * len(rna_DM.obs)
rna_DM.obs['leiden_ori'] = ['DM_%s' % i for i in rna_DM.obs['leiden'].tolist()]
rna_HSR.obs['sample'] = ['HSR'] * len(rna_HSR.obs)
rna_HSR.obs['leiden_ori'] = ['HSR_%s' % i for i in rna_HSR.obs['leiden'].tolist()]

rna = ad.concat([rna_HSR, rna_DM])
print(rna)

sc.pp.neighbors(rna, n_neighbors=10, n_pcs=20)
sc.tl.leiden(rna, resolution=.5)

sc.tl.umap(rna, spread=1., min_dist=.5, random_state=11)
sc.pl.umap(rna, color="leiden", legend_loc="on data", save='_%s_and_%s_rna_leiden.pdf' % (DM_prefix, HSR_prefix))
sc.pl.umap(rna, color="sample", legend_loc="on data", save='_%s_and_%s_rna_sample.pdf' % (DM_prefix, HSR_prefix))
sc.pl.umap(rna, color="leiden_ori", legend_loc="on data", save='_%s_and_%s_rna_leiden_ori.pdf' % (DM_prefix, HSR_prefix))
sc.pl.umap(rna, color="phase", legend_loc="on data", save='_%s_and_%s_rna_phase.pdf' % (DM_prefix, HSR_prefix))

print("DONE!")