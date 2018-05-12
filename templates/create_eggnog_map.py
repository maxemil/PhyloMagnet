#! ${params.python3}

from Bio import SeqIO
from collections import defaultdict
import pickle
from ete3 import ncbi_taxonomy
from ETE3_Utils import defaultdict_string, defaultdict_defaultdict


def get_clade_dict(tax_map):
    taxon_clade_mod = defaultdict(str)
    for k, v in tax_map.items():
        lineage = ncbi.get_lineage(v)
        lineage_string = ""
        for l in lineage:
            rank = ncbi.get_rank([l])[l]
            if rank != 'no rank':
                lineage_string = "{} {};".format(lineage_string, ncbi.get_taxid_translator([l])[l])
        taxon_clade_mod[k] = lineage_string
    return taxon_clade_mod


eggnog_id = "${fasta.simpleName}"
local_eggnog_map = {}
with open("${megan_eggnog_map}") as handle:
    for line in handle:
        if line.startswith('-'):
            continue
        else:
            line = line.split()
            local_eggnog_map[line[1]] = line[0]

ncbi = ncbi_taxonomy.NCBITaxa()

tax_map = defaultdict(defaultdict_string)
basename_files = "${fasta.baseName}".split('.')[1] if len("${fasta.baseName}".split('.')) > 1 else "${fasta.baseName}"
with open("{}_eggnog.map".format(basename_files), 'w') as eggnog_map, \
        open("{}.taxid.map".format(basename_files), 'w') as tax_map_pickel, \
        open("{}.class".format(basename_files), 'w') as class_map:
    class_map.write("%s\\t%s\\n" % (local_eggnog_map[eggnog_id], basename_files))
    for rec in SeqIO.parse("${fasta}", 'fasta'):
        s = rec.id.split('.')
        eggnog_map.write("%s\\t%s\\n" % (s[1], local_eggnog_map[eggnog_id]))
        tax_map[rec.id] = int(s[0])
    taxon_clade_mod = get_clade_dict(tax_map)
    for k,v in taxon_clade_mod.items():
        print("{}\t{}".format(k,v), file=tax_map_pickel)
