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

cellcycle_range = [-0.2, -0.1, 0, 0.1, 0.2]
# [:,-0.2]: G1
# [-0.2, -0.1]: G1/S
# [-0.1, 0]: S (minor G1)
# [0, 0.1]: S (minor G2)
# [0.1, 0.2]: G2 (minor S)
# [0.2, :]: G2

data_counts = pd.read_csv('%s_window1e6_sliding2e5_chr8_counts.txt' % prefix, na_values=['.'], sep='\t')
data_log2FC = pd.read_csv('%s_window1e6_sliding2e5_chr8_log2FC.txt' % prefix, na_values=['.'], sep='\t')
data_counts['cell'] = [x.split('#')[1] for x in data_counts.index]
data_log2FC['cell'] = [x.split('#')[1] for x in data_log2FC.index]
mdata = mu.read("%s.h5mu" % prefix)
rna = mdata.mod['rna']

"""mu.pp.filter_obs(rna, "leiden", lambda x: x.isin(["0", "1", "2", "6", "7"]))
sorted = '01267'
rna = rna[~(rna.obs.leiden.isin(["1", "2"]) & (rna.obs.phase.isin(['G2M'])))].copy()"""
rna_sorted = rna.copy()

cellcycle = pd.DataFrame(columns=['G1', 'G1/S', 'S/G1', 'S/G2', 'G2/S', 'G2/M', 'region'])
cellcycle_heatmap = pd.DataFrame(columns=['G1', 'G1/S', 'S/G1', 'S/G2', 'G2/S', 'G2/M'])

regions = data_counts.columns[:-1]
print(regions)

for i in list(regions):
    print(i)
    rna = rna_sorted.copy()
    rna.obs['counts'] = [data_counts[data_counts['cell'] == rna.obs.index[x]][i].tolist()[0]
                         if rna.obs.index[x] in data_counts['cell'].tolist() else 0 for x in range(len(rna.obs.index))]
    rna.obs['log2FC'] = [data_log2FC[data_log2FC['cell'] == rna.obs.index[x]][i].tolist()[0]
                         if rna.obs.index[x] in data_log2FC['cell'].tolist() else 0 for x in range(len(rna.obs.index))]

    bin1 = rna.obs[rna.obs['G2M_score'] < -0.2]['counts'].median()
    bin2 = rna.obs[(rna.obs['G2M_score'] >= -0.2) & (rna.obs['G2M_score'] < -0.1)]['counts'].median()
    bin3 = rna.obs[(rna.obs['G2M_score'] >= -0.1) & (rna.obs['G2M_score'] < 0)]['counts'].median()
    bin4 = rna.obs[(rna.obs['G2M_score'] >= 0) & (rna.obs['G2M_score'] < 0.1)]['counts'].median()
    bin5 = rna.obs[(rna.obs['G2M_score'] >= 0.1) & (rna.obs['G2M_score'] < 0.2)]['counts'].median()
    bin6 = rna.obs[rna.obs['G2M_score'] >= 0.2]['counts'].median()

    cellcycle.loc[len(cellcycle.index)] = [bin1, bin2, bin3, bin4, bin5, bin6, i]
    cellcycle_heatmap.loc[len(cellcycle_heatmap.index)] = [bin1, bin2, bin3, bin4, bin5, bin6]

cellcycle.to_csv('%s%s_window1e6_sliding2e5_chr8_cellcycle.txt' % (output_dir, prefix), index=False, sep='\t')
cellcycle_heatmap.to_csv('%s%s_window1e6_sliding2e5_chr8_cellcycle_heatmap.txt' % (output_dir, prefix), index=False, sep='\t')

print("DONE!")



