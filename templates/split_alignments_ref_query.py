#! ${params.python3}
from Bio import SeqIO

def run(ref_aln, ref_query_aln, query_aln):
    refs = [rec.id for rec in SeqIO.parse(ref_aln, 'fasta')]
    queries = [rec for rec in SeqIO.parse(ref_query_aln, 'fasta') if not rec.id in refs]
    with open(query_aln, 'w') as out:
        for rec in queries:
            SeqIO.write(rec, out, 'fasta')

run("$refalignment", "$queryalignment", "${queryalignment.baseName}.queries.aln")
