# follow directions from https://www.youtube.com/watch?v=ONiWugVEf2s
# ipythonNotebook is downloaded under Study-Go analysis

# create background list
# https://www.ncbi.nlm.nih.gov/gene/
# "9606" [Taxonomy ID] AND alive[property] AND genetype protein coding[Properties]
# for mouse, use 10090

# python /Users/xwyan/opt/anaconda3/envs/singlecell/bin/ncbi_gene_results_to_python.py -o genes_ncbi_human_proteincoding.py ./gene_result.txt

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
from scipy.stats import ks_2samp
from shared.genes_ncbi_human_proteincoding import GENEID2NT as GeneID2nt
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

data_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
output_dir = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221106_analysis_scMultiome_ColoDM-ColoHSR/COLO320DM_5k/03_analysis/'
os.chdir(data_dir)
rep = 'rep1'
prefix = 'COLO320DM_5K_%s' % rep

# obo_fname = download_go_basic_obo()
# fin_gene2go = download_ncbi_associations()
obodag = GODag("go-basic.obo")

mapper = {}
for key in GeneID2nt:
    mapper[GeneID2nt[key].Symbol] = GeneID2nt[key].GeneID

inv_map = {v: k for k, v in mapper.items()}

df = pd.read_csv('%s%s_average-expression_enrich_in_C2_log2FC_low.txt' % (data_dir, prefix), na_values=['.'], sep='\t')
genelist = df['gene'][:50]

for gene in genelist:
    print(gene)
    try:
        print(GeneID2nt[mapper[gene]].map_location)
    except:
        print("not protein coding")



