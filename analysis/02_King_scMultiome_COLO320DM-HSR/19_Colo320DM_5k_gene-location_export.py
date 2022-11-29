import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
from scipy.stats import ks_2samp
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

data_dir1 = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/' \
           '02_cellranger_output_King/COLO320DM_5k_rep1_hg38/'
data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

# LOAD DATA
mdata = mu.read_10x_h5(os.path.join(data_dir1, "filtered_feature_bc_matrix.h5"))
mdata.var_names_make_unique()

# RNA
rna = mdata.mod['rna']
print(rna)
df_gene = pd.DataFrame(rna.var)
df_gene['gene'] = df_gene.index
df_gene['chr'] = [i.split(':')[0] for i in df_gene['interval'].tolist()]
print(df_gene)
df_gene['start'] = [i.split(':')[1].split('-')[0] if i != 'NA' else np.nan for i in df_gene['interval'].tolist()]
df_gene['end'] = [i.split(':')[1].split('-')[1] if i != 'NA' else np.nan for i in df_gene['interval'].tolist()]
print(df_gene)
df_gene.to_csv('%s%s_gene.txt' % (output_dir, prefix), index=False, sep='\t')

print("DONE!")
