import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import muon as mu
import os

# This is the directory where those files are downloaded to
data_dir1 = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/' \
           '02_cellranger_output_King/COLO320DM_5k_rep1_hg38/'
data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = data_dir
os.chdir(data_dir)

rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

mdata = mu.read_10x_h5(os.path.join(data_dir1, "filtered_feature_bc_matrix.h5"))
mdata.var_names_make_unique()
rna = mdata.mod['rna']
rna.var['mt'] = rna.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.normalize_total(rna, target_sum=1e4)

data_counts = pd.read_csv('%s_window1e6_sliding2e5_chr8_counts.txt' % prefix, na_values=['.'], sep='\t')
data_log2FC = pd.read_csv('%s_window1e6_sliding2e5_chr8_log2FC.txt' % prefix, na_values=['.'], sep='\t')
data_counts['cell'] = [x.split('#')[1] for x in data_counts.index]
data_log2FC['cell'] = [x.split('#')[1] for x in data_log2FC.index]

C2 = pd.read_csv('%s%s_sorted_245.txt' % (data_dir, prefix), header=None, na_values=['.'], sep='\t')[0].tolist()
data = pd.read_csv('%s_normalized_average_expression_enrich_in_C2_log2FC_low_all.txt' % prefix, na_values=['.'], sep='\t')
rna_C2 = rna[rna.obs.index.isin(C2)].copy()
df_rna = pd.DataFrame.sparse.from_spmatrix(rna_C2.X)
df_rna.columns = rna_C2.var.index
df_rna.index = rna_C2.obs.index
df_rna['log2FC'] = [data_log2FC[data_log2FC['cell'] == df_rna.index[x]]['X127600001'].tolist()[0]
                    if df_rna.index[x] in data_log2FC['cell'].tolist() else 0 for x in range(len(df_rna.index))]
df_filter = df_rna[df_rna['log2FC'] > 0].copy()

"""for i in data['gene'].tolist()[:50]:
    plt.scatter(df_filter[i], df_filter['log2FC'], s=5, alpha=0.5)
    plt.xlabel(i)
    plt.ylabel('log2FC')
    plt.savefig('%sfigures/log2FC_scatter_%s.pdf' % (output_dir, i))
    plt.show()"""

for i in data['gene'].tolist()[:50]:
    df_wo0 = df_filter[df_filter[i] > 0].copy()
    plt.scatter(df_wo0[i], df_wo0['log2FC'], s=5, alpha=0.5)
    plt.xlabel(i)
    plt.ylabel('log2FC')
    plt.savefig('%sfigures/log2FC_scatter_wo0_%s.pdf' % (output_dir, i))
    plt.show()

print("DONE!")