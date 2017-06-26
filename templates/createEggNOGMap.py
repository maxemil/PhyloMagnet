#! ${params.python3}

from Bio import SeqIO
from pathlib import Path

eggnog_id = "${fasta.baseName}"
local_eggnog_map = {}
with open("${megan_eggnog_map}") as handle:
  for line in handle:
    if line.startswith('-'):
      continue
    else:
      line = line.split()
      local_eggnog_map[line[1]] = line[0]

with open("${fasta.baseName}_eggnog.map", 'w') as eggnog_map, \
      open("${fasta.baseName}_taxid.map", 'w') as tax_map, \
      open("${fasta.baseName}.class", 'w') as class_map:
  class_map.write("%s\\t%s\\n" % (local_eggnog_map[eggnog_id], eggnog_id))
  for rec in SeqIO.parse("${fasta}", 'fasta'):
    s = rec.id.split('.')
    eggnog_map.write("%s\\t%s\\n" % (s[1], local_eggnog_map[eggnog_id]))
    tax_map.write("%s\\t%s\\n" % (rec.id, s[0]))
