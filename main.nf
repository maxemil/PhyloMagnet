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
def startup_message() {
    revision = grab_git_revision() ?: 'no idea... didnt find the git repo'
    log.info "=========================================================="
    log.info "                       PhyloMagnet"
    log.info "Author                            : Max Emil Schön"
    log.info "email                             : max-emil.schon@icm.uu.se"
    log.info "version                           : $revision"
    log.info "=========================================================="
    log.info "List of EggNOG classes            : $params.reference_classes"
    log.info "List of BioProject Ids            : $params.project_list"
    log.info "Output dir for queries            : $params.queries_dir"
    log.info "Output dir for references         : $params.reference_dir"
    log.info "Phylogenetic method               : $params.phylo_method"
    log.info "number of threads                 : $params.cpus"
    log.info "Map of EggNOG from MEGAN          : $params.megan_eggnog_map"
    log.info "taxonomic level to be analysed    : $params.taxonomy_level_trees"
    log.info "Use a local fastq file            : $params.fastq"
    log.info "Use run IDs instead of projects   : $params.is_runs"
    log.info "location of gc-assembler          : $params.gc_assembler"
    log.info "location of daa-meganizer         : $params.daa_meganizer"
    log.info "Python 3 binary used              : $params.python3"
    log.info ""
    log.info ""
}

startup_message()

// reads a list of Bioproject IDs, but testing only on one single ID
IDs = Channel.from(file(params.project_list).readLines())
// IDs = Channel.from('PRJNA104935')
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

    script:
    if (params.is_runs){
      """
      echo $projectID > runs.txt
      """
    }else{
      template 'filterRuns.py'
    }
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

    script:
    template 'createEggNOGMap.py'
}

/*
  concatenate single mapping files and keep track of long/short IDs for EggNOG
*/
eggnog_map_concat = eggnog_map.collectFile(name: 'eggnog.syn', storeDir: params.reference_dir)
// tax_map_concat = tax_map.collectFile(name: 'tax.syn', storeDir: params.reference_dir)
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

    script:
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

    script:
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
    script:
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

    script:
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

    script:
    template 'translateDNAtoAA.py'
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
    cpus 4
    // maxForks 1

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
      #iqtree-omp -s ${aln} -m LG -bb 1000 -nt 2 -pre ${aln.baseName}
      /local/two/Software/iqtree-1.6.beta3-Linux/bin/iqtree -fast -s ${aln} -m LG -nt 2 -pre ${aln.baseName}
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
    file '*' from tax_map.collect()

    output:
    file "${tree.baseName}.pdf" into pdfs
    file 'decision.txt' into decisions
    stdout x

    publishDir "${params.queries_dir}/${tree.baseName.minus(~/-.+/)}", mode: 'copy'

    beforeScript = {"ln -s \$(grep ${tree.baseName.minus(~/^.+-/).minus(~/.trim/)} $workflow.launchDir/${params.reference_dir}/class.map | cut -f 2)\"_taxid.map\" tax.map"}

    script:
    template 'makePDFfromTree.py'
}
x.subscribe{print it}
decisions_concat = decisions.collectFile(name: 'decisions.txt', storeDir: params.queries_dir)
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


def grab_git_revision() {
    if ( workflow.commitId ) { // it's run directly from github
        return workflow.commitId
    }

    // Try to find the revision directly from git
    head_pointer_file = file("${workflow.projectDir}/../.git/HEAD")
    if ( ! head_pointer_file.exists() ) {
        return ''
    }
    ref = head_pointer_file.newReader().readLine().tokenize()[1]

    ref_file = file("${workflow.projectDir}/../.git/$ref")
    if ( ! ref_file.exists() ) {
        return ''
    }
    revision = ref_file.newReader().readLine()

    return revision
}
