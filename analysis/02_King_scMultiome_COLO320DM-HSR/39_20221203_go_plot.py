import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import seaborn as sns
import textwrap
import os

data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
os.chdir(data_dir)
rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

df = pd.read_csv('%s_go_enrich_in_C2_log2FC_low_local_low_top50.txt' % prefix, na_values=['.'], sep='\t')
df['per'] = df.n_genes/df.n_go
df = df[0:10]

fig, ax = plt.subplots(figsize=(2, 9))
ax = fig.add_axes([0.15, 0.12, 0.4, 0.75])
cmap = mpl.cm.bwr_r
norm = mpl.colors.Normalize(vmin=df.p_corr.min(), vmax=df.p_corr.max())
mapper = cm.ScalarMappable(norm=norm, cmap=cm.bwr_r)
cbl = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
plt.savefig('%sfigures/go_plot_low_top50_palette.pdf' % output_dir)
plt.show()

fig, ax = plt.subplots(figsize=(9, 9))
fig.subplots_adjust(left=0.2)
ax = sns.barplot(data=df, x='per', y='term', palette=mapper.to_rgba(df.p_corr.values))
ax.set_yticklabels([textwrap.fill(e, 22) for e in df['term']])
plt.savefig('%sfigures/go_plot_low_top50.pdf' % output_dir)
plt.show()

