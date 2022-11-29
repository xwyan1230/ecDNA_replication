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
fin_gene2go = download_ncbi_associations()
obodag = GODag("go-basic.obo")

mapper = {}
for key in GeneID2nt:
    mapper[GeneID2nt[key].Symbol] = GeneID2nt[key].GeneID

inv_map = {v: k for k, v in mapper.items()}

# read ncbi's gene2go, store annotations in a list of namedtuples
objanno = Gene2GoReader(fin_gene2go, taxids=[9606])
# get namespace2association where:
#   namespace is:
#       BP: biological_process
#       MF: molecular_function
#       CC: cellular_component
#   association is a dict:
#       key: NCBI geneID
#       value: a set of GO IDs associated with that gene
ns2assoc = objanno.get_ns2assc()

goeaobj = GOEnrichmentStudyNS(
    GeneID2nt.keys(),  # list of human protein coding genes
    ns2assoc,  # geneID/GO associations
    obodag,  # ontologies
    propagate_counts=False,
    alpha=0.05,  # default significance cut-off
    methods=['fdr_bh'])  # default multipletest correction method

GO_items = []

temp = goeaobj.ns2objgoea['BP'].assoc
for item in temp:
    GO_items += temp[item]
temp = goeaobj.ns2objgoea['CC'].assoc
for item in temp:
    GO_items += temp[item]
temp = goeaobj.ns2objgoea['MF'].assoc
for item in temp:
    GO_items += temp[item]

# you can count number of genes associated with a certain go term by
# GO_items.count('GO:0001525')

df = pd.read_csv('%s%s_average-expression_enrich_in_C2_log2FC_low.txt' % (data_dir, prefix), na_values=['.'], sep='\t')
genelist_ori = df['gene'][:50]

# gene_remove = ['PVT1', 'PCAT1', 'LRATD2', 'MYC', 'CASC8', 'PRNCR1', 'AC084116.2']
# genelist = [e for e in genelist_ori if e not in gene_remove]
genelist = genelist_ori
print(len(genelist))


def go_it(test_genes):
    print(f'input genes: {len(test_genes)}')

    mapped_genes = []
    for gene in test_genes:
        try:
            mapped_genes.append(mapper[gene])
        except:
            pass
    print(f'mapped genes: {len(mapped_genes)}')

    goea_results_all = goeaobj.run_study(mapped_genes)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    GO = pd.DataFrame(list(map(lambda x: [x.GO, x.goterm.name, x.goterm.namespace, x.p_uncorrected, x.p_fdr_bh,
                                          x.ratio_in_study[0], x.ratio_in_study[1], GO_items.count(x.GO),
                                          list(map(lambda y: inv_map[y], x.study_items))], goea_results_sig)),
                      columns=['GO', 'term', 'class', 'p', 'p_corr', 'n_genes', 'n_study', 'n_go', 'study_genes'])

    GO = GO[GO.n_genes > 1]
    return GO


df_go = go_it(genelist)
df_go.to_csv('%s%s_go_enrich_in_C2_log2FC_low_top50.txt' % (output_dir, prefix), index=False, sep='\t')

