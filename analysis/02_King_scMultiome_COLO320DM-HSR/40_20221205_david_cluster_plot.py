import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import seaborn as sns
import textwrap
import os

data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/01_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/01_analysis/'
os.chdir(data_dir)
rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

df = pd.read_csv('david_stdev_HSR.txt', na_values=['.'], sep='\t')
df['per'] = df['count']/df['total']
df = df[0:10]

fig, ax = plt.subplots(figsize=(2, 9))
ax = fig.add_axes([0.15, 0.12, 0.4, 0.75])
cmap = mpl.cm.bwr_r
norm = mpl.colors.Normalize(vmin=0, vmax=2.5)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.bwr_r)
cbl = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
plt.savefig('%sfigures/david_plot_stdev_palette.pdf' % output_dir)
plt.show()

fig, ax = plt.subplots(figsize=(9, 9))
fig.subplots_adjust(left=0.2)
ax = sns.barplot(data=df, x='per', y='cluster', palette=mapper.to_rgba(df.group_ES.values))
ax.set_yticklabels([textwrap.fill(e, 22) for e in df['cluster']])
plt.savefig('%sfigures/david_plot_stdev_HSR.pdf' % output_dir)
plt.show()

