## inputs
f_loom_grn="../output/mc_mat.loom"
f_tfs="../cisTarget_db/hsa_hgnc_tfs.motifs-v10.txt"
f_motif_path="../cisTarget_db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
f_db_500bp="../cisTarget_db/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
f_db_10kb="../cisTarget_db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

## outputs
grn_output="../output/adj.tsv"
ctx_output="../output/reg.tsv"

#### 1. Build GRN 
time arboreto_with_multiprocessing.py $f_loom_grn $f_tfs --method grnboost2 --output $grn_output --num_workers 40 --seed 777

#### 2. Run cisTarget 
time pyscenic ctx $grn_output $f_db_500bp $f_db_10kb --annotations_fname $f_motif_path --expression_mtx_fname $f_loom_grn --output $ctx_output --num_workers 40 --mask_dropouts

#### 3. Run AUCell 
time pyscenic aucell $f_loom_grn $ctx_output --output ../output/SCENIC.loom --num_workers 40

import sys
import pandas as pd
from pyscenic.cli.utils import load_signatures

regulon_file = "../output/01_reg4.tsv"
project = "../output/micro"
min_regulon_size = 10

def get_motif_logo(regulon):
    base_url = "http://motifcollections.aertslab.org/v10nr_clust/logos/"
    for elem in regulon.context:
        if elem.endswith('.png'):
            return(base_url + elem)

# regulons to gmt file
print('''
##############################################
    1. Transform regulons to gmt file ...
##############################################
    ''')
sys.stdout.flush()
regulons = load_signatures(regulon_file)
select_cols = [i.name for i in regulons if len(i.genes) >= min_regulon_size]
gmt_file = project + ".regulons.gmt"
txt_file = project + ".regulons.txt"
fo1 = open(gmt_file, 'w')
fo2 = open(txt_file, 'w')
for i in regulons:
    if i.name in select_cols:
        motif = get_motif_logo(i)
        genes = "\t".join(i.genes)
        tf = "%s(%sg)" % (i.transcription_factor, len(i.genes))
        fo1.write("%s\t%s\t%s\n" % (tf, motif, genes))
        fo2.write("%s\t%s\t%s\n" % (tf, motif, genes.replace("\t",",")))
