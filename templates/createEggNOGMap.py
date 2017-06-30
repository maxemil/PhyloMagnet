#! ${params.python3}

from Bio import SeqIO
from pathlib import Path
from ETE3_Utils import *
import pickle

eggnog_id = "${fasta.baseName}"
local_eggnog_map = {}
with open("${megan_eggnog_map}") as handle:
    for line in handle:
        if line.startswith('-'):
            continue
        else:
            line = line.split()
            local_eggnog_map[line[1]] = line[0]


def defaultdict_string():
    return ""

tax_map = defaultdict(defaultdict_string)
with open("${fasta.baseName}_eggnog.map", 'w') as eggnog_map, \
        open("${fasta.baseName}_taxid.map", 'wb') as tax_map_pickel, \
        open("${fasta.baseName}.class", 'w') as class_map:
    class_map.write("%s\\t%s\\n" % (local_eggnog_map[eggnog_id], eggnog_id))
    for rec in SeqIO.parse("${fasta}", 'fasta'):
        s = rec.id.split('.')
        eggnog_map.write("%s\\t%s\\n" % (s[1], local_eggnog_map[eggnog_id]))
        tax_map[rec.id] = int(s[0])
    clade_taxon_mod, taxon_clade_mod = get_clade_names(tax_map, "$params.taxonomy_level_trees")
    pickle.dump(taxon_clade_mod, tax_map_pickel, -1)
