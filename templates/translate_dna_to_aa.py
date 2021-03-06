#! ${params.python3}

from Bio import SeqIO
import copy
import os

seqs = []
seqs_codon_table_6 = []
with open('%s.faa' % "${contig.baseName}", 'w') as out:
    for seq in SeqIO.parse("${contig}", 'fasta'):
        seq_prot = copy.deepcopy(seq)
        seq_prot.seq = seq.seq.translate(table=1)
        if '*' in seq_prot.seq:
            seq_prot.seq = seq.seq.translate(table=6)
            seqs_codon_table_6.append(seq_prot)
            continue
        seqs.append(seq_prot)
    SeqIO.write(seqs, out, 'fasta')

if os.stat('%s.faa' % "${contig.baseName}").st_size == 0:
    os.remove('%s.faa' % "${contig.baseName}")
        #SeqIO.write(seqs_codon_table_6, '%s_table6.faa' % "${contig.baseName}", 'fasta')
