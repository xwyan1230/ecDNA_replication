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
# sorted = '01267'

data = pd.read_csv('%s_window10e6_sliding2e6_chr8_cellcycle_heatmap.txt' % prefix, na_values=['.'], sep='\t')

sns.heatmap(data)
plt.savefig('%sfigures/%s_rna_window10e6_sliding2e6_chr8_heatmap.pdf' % (output_dir, prefix))
plt.close()
sns.heatmap(data, vmax=300)
plt.savefig('%sfigures/%s_rna_window10e6_sliding2e6_chr8_heatmap_vmax300.pdf' % (output_dir, prefix))
plt.close()

print("DONE!")
