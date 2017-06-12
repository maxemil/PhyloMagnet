#!/usr/bin/env nextflow

/*
This is the main PhyloMagnet workflow file
author : max emil schön
email  : <max-emil.schon{at}icm.uu.se>
dependencies:
  python3
  Biopython
  diamond
  megan >= version 6 (gc-assembler and daa-meganizer tools)
  mafft >= version 7 (--addfragments option)
  trimal
  FastTree or iqtree-omp
*/

params.daa_meganizer = "/usr/local/bin/megan6-ce/tools/daa-meganizer"
params.gc_assembler = "/usr/local/bin/megan6-ce/tools/gc-assembler"
//params.python3 = "/usr/bin/env python3"
/* this is just because of my weird anaconda/python installation */
params.python3 = "/usr/local/bin/anapy3"
// can be either "fasttree" or "iqtree-omp"
params.phylo_method = "fasttree"
// should be 1 by default?
params.cpus = "40"
params.megan_eggnog_map = "eggnog.map"
params.reference_classes = "EGGNOG_List"
params.reference_dir = "references/"
params.queries_dir = "queries/"
params.project_list = "bioproject_result.txt"
params.taxonomy_level_trees = "class"

params.fastq = ""

// reads a list of Bioproject IDs, but testing only on one single ID
// IDs = Channel.from(file(params.project_list).readLines())
IDs = Channel.from('PRJNA104935')
EggNOGIDs = Channel.from(file(params.reference_classes).readLines())
megan_eggnog_map = Channel.from(file(params.megan_eggnog_map))

/*******************************************************************************
******************** Download and Prepare Section ******************************
*******************************************************************************/

/*
  download a description file from BioProject and filter all runs that seem to
  be shotgun metagenomic from the Illumina platorm
  Main process is in python3
*/
process filterRuns {
    input:
    val projectID from IDs

    output:
    set file('runs.txt'), val(projectID) into project_runs

    errorStrategy 'ignore'

    """
    #! ${params.python3}

    import re
    import sys
    import requests

    url = "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term=${projectID}"
    desc = requests.get(url).text

    regex_strings = ['WGS|shotgun', 'ILLUMINA', 'METAGENOMIC|metagenome']
    with open('runs.txt', 'w') as out:
      for line in desc.split():
        if all(re.search(pat, line) for pat in regex_strings):
          print(line.split(",")[0], file=out)
          #open(line.split(",")[0], 'w').close()
          #print("%s_%s" % ("${projectID}", line.split(",")[0]))
    """
}

/*
  for all valid runs specified in 'runs.txt', download fastQ files (in parallel).
  split reads if paired end data, but keep all in one file for diamond
  alignment and assembly
*/
process downloadFastQ {
    input:
    set file('runs.txt'), val(x) from project_runs

    output:
    file "*.fastq.gz" into fastq_files mode flatten

    publishDir params.queries_dir

    script:
    fastq_file = new File(params.fastq)
    if (fastq_file.exists())
      """
      ln -s ${fastq_file.getAbsolutePath()} .
      """
    else
      """
      #! /usr/bin/env bash
      while read run; do
        fastq-dump --gzip --readids --split-spot --skip-technical --clip \$run &
      done <  runs.txt
      wait
      """
}

/*
  Download all raw sequence files as well as the untrimmed alignment
  files for all provided EggNOG Ids.
*/
process downloadEggNOG {
    input:
    val id from EggNOGIDs

    output:
    file "${id}.fasta" into EggNOGFastas
    file "${id}.aln" into EggNOGAlignments

    publishDir params.reference_dir, mode: 'copy', overwrite: false

    """
    wget http://eggnogapi.embl.de/nog_data/text/fasta/${id} -O - | gunzip > ${id}.fasta
    wget http://eggnogapi.embl.de/nog_data/text/raw_alg/${id} -O - | gunzip > ${id}.aln
    """
}

/*
  redirect EggNOG fastA files to a concatenation process and to the map
  building process needed for meganization.
  Concatenate all fastA into one single references.fasta for diamond alignment
*/
EggNOGFastas.into {EggNOGFastas_concatenation; EggNOGFastas_mapping }
concatenated_references = EggNOGFastas_concatenation.collectFile(name: 'references.fasta', storeDir: params.reference_dir)

/*
  create a mapping file of reference sequence ID to EggNOG ID that can be
  used as synonyms file in megan
*/
process createEggNOGMap {
    input:
    file megan_eggnog_map from megan_eggnog_map.first()
    file fasta from EggNOGFastas_mapping

    // TODO dump eggnog_map or taxmap as python3 binary? to avoid parsing it again later
    output:
    file "${fasta.baseName}_eggnog.map" into eggnog_map
    file "${fasta.baseName}_taxid.map" into tax_map
    file "${fasta.baseName}.class" into class_map

    """
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
    """
}

/*
  concatenate single mapping files and keep track of long/short IDs for EggNOG
*/
eggnog_map_concat = eggnog_map.collectFile(name: 'eggnog.syn', storeDir: params.reference_dir)
tax_map_concat = tax_map.collectFile(name: 'tax.syn', storeDir: params.reference_dir)
class_map_concat = class_map.collectFile(name: 'class.map', storeDir: params.reference_dir)

/*
  prepare the diamond database from the concatenated fastA file
*/
process diamondMakeDB {
    input:
    file 'references.fasta' from concatenated_references

    output:
    file 'references.dmnd' into diamond_database

    publishDir params.reference_dir

    """
    diamond makedb --in references.fasta --db references.dmnd
    """
}

/*******************************************************************************
****************************** Assembly Section ********************************
*******************************************************************************/

/*
  Align each fastQ file against the diamond database, create daa files
*/
process alignFastQFiles {
    input:
    file fq from fastq_files
    file 'references.dmnd' from diamond_database.first()

    output:
    file "${fq.getName().minus(".fastq.gz")}.daa" into daa_files

    //publishDir 'queries', mode: 'copy'
    cpus = params.cpus

    """
    diamond blastx -q ${fq} --db references.dmnd -f 100 --unal 0 -e 1e-6 --out ${fq.getName().minus(".fastq.gz")}.daa
    """
}


/*
  Prepare daa files for gc-assembler, using the daa-meganizer tool and
  the synonyms mapping file created earlier, applies to the daa files in place,
  so copy them
  TODO may cause severe memory overhead in the future as we copy files a lot
*/
process meganizeDAAFiles {
    stageInMode 'copy'

    input:
    file daa from daa_files
    file eggnog_map from eggnog_map_concat.first()

    output:
    file daa into daa_files_meganized

    //publishDir 'queries', mode: 'copy', overwrite: true
    """
    mkdir references
    mv ${eggnog_map} references/
    ${params.daa_meganizer} --in ${daa} -fun EGGNOG -s2eggnog references/${eggnog_map}
    """
}


/*
  Perform Gene Centric Assembly for all classes of EggNOG. this will create
  fastA files of the form <runID>-<shortEggNOGID>.fasta
*/
process geneCentricAssembly {

    input:
    file daa from daa_files_meganized

    output:
    file '*.fasta' into assembled_contigs

    publishDir "${params.queries_dir}/${daa.baseName}", mode: 'copy'

    """
    ${params.gc_assembler} -i ${daa} -fun EGGNOG -id ALL -mic 99 -v --minAvCoverage 2
    """
}

/*******************************************************************************
************************* Alignment and Phylogeny Section **********************
*******************************************************************************/

/*
  Simple Biopython script to translate DNA to AA, first try table one, then table 6 if a stop codon is in sequence
  TODO include checks for this, but how can i make sure i translate correctly?
*/
process translateDNAtoAA {
    input:
    file contig from assembled_contigs.flatten()

    output:
    file '*.faa' into translated_contigs

    publishDir "${params.queries_dir}/${contig.baseName.minus(~/-.+/)}", mode: 'copy'

    """
    #! ${params.python3}

    from Bio import SeqIO
    import copy

    seqs = []
    seqs_codon_table_6 = []
    for seq in SeqIO.parse("${contig}", 'fasta'):
        seq_prot = copy.deepcopy(seq)
        seq_prot.seq = seq.seq.translate(table=1)
        if '*' in seq_prot.seq:
            seq_prot.seq = seq.seq.translate(table=6)
            seqs_codon_table_6.append(seq_prot)
            continue
        seqs.append(seq_prot)
        SeqIO.write(seqs, '%s.faa' % "${contig.baseName}", 'fasta')
        #SeqIO.write(seqs_codon_table_6, '%s_table6.faa' % "${contig.baseName}", 'fasta')
    """
}


/*
  All assembled contigs are aligned to the corresponding EggNOG alignment.
  This is relatively fast using the mafft --addfragments option.
*/
process alignContigs {
    input:
    file faa from translated_contigs
    file reference_alignment from EggNOGAlignments.toList()
    file class_map_concat from class_map_concat.first()

    output:
    file '*.aln' into aligned_contigs

    //publishDir 'queries', mode: 'copy'
    cpus params.cpus

    """
    id=\$(grep "${faa.baseName.minus(~/^.+-/)}" $class_map_concat | cut -f 2)".aln"
    mafft-fftnsi --adjustdirection --thread ${params.cpus} --addfragments ${faa} \$id > ${faa.baseName}.aln
    """
}


/*
  Trim all alignments of EggNOG references with contigs.
  Uses -gappyout for simplicity
*/
process trimContigAlignments {
    input:
    file aln from aligned_contigs.flatten()

    output:
    file "${aln.baseName}.trim.aln" into trimmed_contig_alignments

    publishDir "${params.queries_dir}/${aln.baseName.minus(~/-.+/)}", mode: 'copy'

    """
    trimal -in ${aln} -out ${aln.baseName}.trim.aln -gappyout
    """
}


/*
  Make trees for all alignments, either usig fasttree, iqtree-omp or RAxML
  try producing the same file name containing the tree with support values
  "<RUN-ID>-<EggNOG-ID>.trim.treefile"
*/
process buildTreeFromAlignment {
    input:
    file aln from trimmed_contig_alignments

    output:
    file "${aln.baseName}.treefile" into trees

    publishDir "${params.queries_dir}/${aln.baseName.minus(~/-.+/)}", mode: 'copy'
    cpus 20

    script:
    if (params.phylo_method == "iqtree")
      """
      iqtree-omp -s ${aln} -m LG -bb 1000 -nt 2 -pre ${aln.baseName}
      """
    else if (params.phylo_method == "fasttree")
      """
      FastTree ${aln} > ${aln.baseName}.treefile
      """
}

trees.into{ treesVisualize; treesMagnetize }

process makePDFsFromTrees {
    input:
    file tree from treesVisualize
    file tax_map_concat from tax_map_concat.first()

    output:
    file "${tree.baseName}.pdf" into pdfs

    publishDir "${params.queries_dir}/${tree.baseName.minus(~/-.+/)}", mode: 'copy'

    script:
    """
    #! ${params.python3}
    import sys
    sys.path.append('/local/two/Software/python_lib/')
    from ETE3_Utils import *
    from ETE3_styles import *
    from misc_utils import *
    import ete3

    def my_initiate_clades(tree, taxon_clade):
      for l in tree.traverse():
        l.add_feature(pr_name='clade', pr_value=taxon_clade[l.name])

    tree = parse_newick("$tree")
    set_node_style(tree, node_style_basic)
    taxon_clade, clade_taxon, taxon_prefix, prefix_taxon = read_prefix_map(tree, "$tax_map_concat")
    clade_taxon_mod, taxon_clade_mod = get_clade_names(taxon_clade, "$params.taxonomy_level_trees")
    my_initiate_clades(tree, taxon_clade_mod)

    leaves = set(tree.iter_leaves())
    while leaves:
      node = leaves.pop()
      mono_clade = get_mono_clade(node)
      ancestor = get_ancestor(list(mono_clade))
      if len(mono_clade) > 1 and node.clade:
        ancestor.name = node.clade
        for l in mono_clade:
          leaves.discard(l)
        for c in ancestor.get_children():
          ancestor.remove_child(c)

    tree.ladderize(direction=1)
    ts = tree_style_basic(layout_node_color, "${tree.baseName}")
    tree.render("${tree.baseName}.pdf", tree_style=ts)#w=1500, units="px",
    """
}

// process magnetizeTrees{
//     input:
//     file tree from treesMagnetize
//
//     output:
//
//     script:
//     """
//     #! ${params.python3}
//
//     tree = parse_newick("$tree")
//
//     """
//
//
// }
