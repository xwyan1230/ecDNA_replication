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
# sorted = '245'

data = pd.read_csv('%s_sorted_245_chr8_120000001-130000001_window1e5_sliding2e4_cellcycle_heatmap.txt' % prefix, na_values=['.'], sep='\t')

sns.heatmap(data)
plt.savefig('%sfigures/%s_rna_sorted_245_chr8_120000001-130000001_window1e5_sliding2e4_cellcycle_heatmap.pdf' % (output_dir, prefix))
plt.close()
# sns.heatmap(data, vmax=50)
# plt.savefig('%sfigures/%s_rna_chr8_120000001-130000001_window1e5_sliding2e4_cellcycle_heatmap.pdf' % (output_dir, prefix))
# plt.close()

print("DONE!")
