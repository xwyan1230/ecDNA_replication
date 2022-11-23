import scanpy as sc
import pandas as pd
import muon as mu
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# This is the directory where those files are downloaded to
data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_scMultiome_ColoDM-ColoHSR/COLO320HSR_5k/03_analysis/'
output_dir = data_dir
os.chdir(data_dir)

rep = 'rep1'
prefix = 'COLO320HSR_5K_%s' % rep

cellcycle_to_int = {
    "G1": 0,
    "S": 1,
    "G2M": 2,
}

data_counts = pd.read_csv('%s_window1e6_sliding2e5_chr8_counts.txt' % prefix, na_values=['.'], sep='\t')
data_log2FC = pd.read_csv('%s_window1e6_sliding2e5_chr8_log2FC.txt' % prefix, na_values=['.'], sep='\t')
data_counts['cell'] = [x.split('#')[1] for x in data_counts.index]
data_log2FC['cell'] = [x.split('#')[1] for x in data_log2FC.index]
mdata = mu.read("%s.h5mu" % prefix)
rna = mdata.mod['rna']
counts = []
rna.obs['counts'] = [data_counts[data_counts['cell'] == rna.obs.index[x]]['X127600001'].tolist()[0]
                     if rna.obs.index[x] in data_counts['cell'].tolist() else 0 for x in range(len(rna.obs.index))]
rna.obs['log2FC'] = [data_log2FC[data_log2FC['cell'] == rna.obs.index[x]]['X127600001'].tolist()[0]
                     if rna.obs.index[x] in data_log2FC['cell'].tolist() else 0 for x in range(len(rna.obs.index))]
rna.obs['phase_num'] = [cellcycle_to_int[i] for i in rna.obs['phase'].tolist()]
print(rna)

sc.pl.umap(rna, color="leiden", legend_loc="on data", save='_%s_rna_leiden.pdf' % prefix)
sc.pl.umap(rna, color="phase", legend_loc="on data", save='_%s_rna_phase.pdf' % prefix)
sc.pl.umap(rna, color="S_score", legend_loc="on data", save='_%s_rna_S_score.pdf' % prefix)
sc.pl.umap(rna, color="G2M_score", legend_loc="on data", save='_%s_rna_G2M_score.pdf' % prefix)
sc.pl.umap(rna, color="log2FC", legend_loc="on data", save='_%s_rna_window1e6_sliding2e5_chr8_X127600001_log2FC.pdf' % prefix)
sc.pl.umap(rna, color="counts", legend_loc="on data", vmin=10, vmax=500, save='_%s_rna_window1e6_sliding2e5_chr8_X127600001_counts.pdf' % prefix)
sc.pl.umap(rna, color="MKI67", legend_loc="on data", save='_%s_rna_MKI67.pdf' % prefix)
sc.pl.umap(rna, color="POLR3D", legend_loc="on data", save='_%s_rna_POLR3D.pdf' % prefix)
sc.pl.umap(rna, color="BRD4", legend_loc="on data", save='_%s_rna_BRD4.pdf' % prefix)

mu.pp.filter_obs(rna, "leiden", lambda x: x.isin(["0", "1", "2", "6", "7"]))
rna = rna[~(rna.obs.leiden.isin(["1", "2"]) & (rna.obs.phase.isin(['G2M'])))].copy()
mu.pp.filter_obs(rna, "log2FC", lambda x: x > 1)
sc.pl.umap(rna, color="leiden", legend_loc="on data", save='_%s_rna_sorted_leiden.pdf' % prefix)
sc.pl.umap(rna, color="phase", legend_loc="on data", save='_%s_rna_sorted_phase.pdf' % prefix)
sc.pl.umap(rna, color="S_score", legend_loc="on data", save='_%s_rna_sorted_S_score.pdf' % prefix)
sc.pl.umap(rna, color="G2M_score", legend_loc="on data", save='_%s_rna_sorted_G2M_score.pdf' % prefix)
sc.pl.umap(rna, color="log2FC", legend_loc="on data", save='_%s_rna_sorted_window1e6_sliding2e5_chr8_X127600001_log2FC.pdf' % prefix)
sc.pl.umap(rna, color="counts", legend_loc="on data", vmin=10, vmax=500, save='_%s_rna_sorted_window1e6_sliding2e5_chr8_X127600001_counts.pdf' % prefix)

plt.scatter(rna.obs['G2M_score'].tolist(), rna.obs['counts'].tolist(), c=rna.obs['phase_num'], alpha=0.5, s=3)
plt.ylim([0, 800])
plt.xlabel('G2M_score')
plt.ylabel('counts')
plt.savefig('%sfigures/%s_rna_sorted_G2M_score_vs_window1e6_sliding2e5_chr8_X127600001_counts.pdf' % (output_dir, prefix))
plt.close()

sns.boxplot(data=rna.obs, y='counts', x='phase_num', showfliers=False)
plt.savefig('%sfigures/%s_rna_sorted_window1e6_sliding2e5_chr8_X127600001_counts_by_phase.pdf' % (output_dir, prefix))

print("DONE!")



