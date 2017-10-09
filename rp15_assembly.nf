#!/usr/bin/env nextflow

params.references = ""
params.reads = ""
params.phylo_method = "fasttree"
params.outdir = "results"

references = Channel.fromPath(params.references)
Channel.fromPath(params.reads).into{reads; reads_extraction}


process makeDiamondDb {
  input:
  file '*' from references.collect()

  output:
  file "references.dmnd" into diamond_index

  script:
  """
  cat * > references.fasta
  diamond makedb --in references.fasta --db references.dmnd
  """
}


process alignReads {
  input:
  file index from diamond_index.first()
  file reads

  output:
  file "${reads.baseName}.sam" into aligned_reads

  script:
  """
  diamond blastx -q $reads --more-sensitive \
                 --unal 0 -e 1e-6 -o "${reads.baseName}.sam" \
                 -f sam --db $index
  """
}


process extractAlignedReads {
  input:
  file reads from reads_extraction
  file sam from aligned_reads

  output:
  file "${reads.baseName}.aligned.fastq" into selected_reads
  file "${reads.baseName}.aligned.ids" into selected_ids

  publishDir "${params.outdir}"

  script:
  """
  cut -f1 $sam | sort | uniq | sed "s/\\.[1-2]\$//g" | sort | uniq > ${reads.baseName}.aligned.ids_base
  sed "s/\$/.1/g" ${reads.baseName}.aligned.ids_base > ${reads.baseName}.aligned.ids
  sed "s/\$/.2/g" ${reads.baseName}.aligned.ids_base >> ${reads.baseName}.aligned.ids
  seqtk subseq $reads ${reads.baseName}.aligned.ids > ${reads.baseName}.aligned.fastq
  """
}

process assembleAlignedReads {
  input:
  file selected_reads

  output:
  file "assembly/final.contigs.fa" into assembled_reads

  publishDir "${params.outdir}", mode: 'copy'

  script:
  """
  /local/two/Software/megahit/megahit -r $selected_reads -t 10 --k-min 15 --k-max 151 --k-step 4 -o assembly
  """
}

process splitAssemblyFasta {
  input:
  file assembled_reads

  output:
  file '*.fasta' into contigs_single_fastas

  publishDir "${params.outdir}/assembly", mode: 'copy'

  script:
  """
  #! /usr/local/bin/anapy3
  from Bio import SeqIO

  for rec in SeqIO.parse('final.contigs.fa', 'fasta'):
    with open("%s.fasta" %rec.id, 'w') as out:
      SeqIO.write(rec, out, 'fasta')
  """
}

process predictGenes {
  input:
  file contig from contigs_single_fastas.flatten()

  output:
  file "PROKKA/${contig.baseName}.faa" into translated_contigs

  publishDir "${params.outdir}", mode: 'copy'

  script:
  """
  /local/two/Software/prokka-partial/bin/prokka --outdir PROKKA \
            --partialgenes --prefix ${contig.baseName} \
            --locustag ${contig.baseName} $contig
  """
}
