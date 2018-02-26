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
    log.info "Lineage(s) to look for            : $params.lineage"
    log.info ""
    log.info "Binaries location (use default if singularity image is used)"
    log.info "location of gc-assembler          : $params.gc_assembler"
    log.info "location of daa-meganizer         : $params.daa_meganizer"
    log.info "Python 3 binary used              : $params.python3"
    log.info ""
    log.info ""
}

startup_message()

// reads a list of Bioproject IDs, but testing only on one single ID
if (params.fastq) {
    ids = Channel.from(params.fastq)
}else {
    ids = Channel.from(file(params.project_list).readLines())
}
lineage_list = params.lineage.tokenize(',')
lineage = Channel.from(lineage_list)
// reads a list of Bioproject IDs, but testi

rank = Channel.from(params.taxonomic_rank)

eggNOGIDs = Channel.from(file(params.reference_classes).readLines())
eggNOGIDs_local = Channel.from(file(params.reference_classes))
Channel.from(file(params.megan_eggnog_map)).into { megan_eggnog_map; megan_eggnog_map_local }

if (params.local_ref) {
    local_ref = Channel.fromPath(params.local_ref)
}else {
    local_ref = Channel.empty()
}
// println(file(megan_eggnog_map.first()).isEmpty())

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
    val projectID from ids

    output:
    set file('runs.txt'), val(projectID) into project_runs

    errorStrategy 'ignore'
    tag "$projectID"

    script:
    if (params.is_runs){
      """
      echo $projectID > runs.txt
      """
    }else{
      template 'filter_runs.py'
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
    maxForks 2
    cpus 1
    tag "$x"

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

process includeLocalRef  {
  input:
  file "*" from local_ref.collect()
  file ids from eggNOGIDs_local.first()
  file megan_eggnog_map from megan_eggnog_map_local.first()

  output:
  file "*.fasta" into local_ref_eggnog mode flatten
  file "local_translation.txt" into local_ref_translation

  publishDir params.reference_dir, mode: 'copy'
  cache 'deep'

  script:
  """
  readlink -f $ids
  for f in *.fasta;
  do
    ID=\$(comm -3 <(grep -v '^-' $megan_eggnog_map | cut -f2 | cut -d ' ' -f 1 | sort ) <(cat $ids | sort) | tr -d '\\t' | head -n 1)
    cp -L \$f \$ID.fasta
    echo \$ID >> $ids
    echo \$(basename \$f)"\t"\$ID >> local_translation.txt
  done
  """
}

local_ref_eggnog.into {local_ref_eggnog_align; local_ref_eggnog_add}


process alignLocalRef {
  input:
  file fasta from local_ref_eggnog_align

  output:
  file "${fasta.baseName}.aln" into local_ref_align

  publishDir params.reference_dir, mode: 'copy'
  tag "${fasta.baseName}"

  script:
  """
  mafft-linsi --adjustdirection --thread ${task.cpus} $fasta > ${fasta.baseName}.aln
  """
}

/*
  Download all raw sequence files as well as the untrimmed alignment
  files for all provided EggNOG Ids.
*/
process downloadEggNOG {
    input:
    val id from eggNOGIDs

    output:
    file "${id}.fasta" into eggNOGFastas
    file "${id}.aln" into eggNOGAlignments

    publishDir params.reference_dir, mode: 'copy', overwrite: false
    tag "$id"

    """
    wget http://eggnogapi.embl.de/nog_data/text/fasta/${id} -O - | gunzip > ${id}.fasta || wget http://eggnogapi.embl.de/nog_data/text/fasta/${id} -O ${id}.fasta
    wget http://eggnogapi.embl.de/nog_data/text/raw_alg/${id} -O - | gunzip > ${id}.aln || wget http://eggnogapi.embl.de/nog_data/text/raw_alg/${id} -O ${id}.aln
    """
}

/*
  redirect EggNOG fastA files to a concatenation process and to the map
  building process needed for meganization.
  Concatenate all fastA into one single references.fasta for diamond alignment
*/
ref_alignments = eggNOGAlignments.mix(local_ref_align)
eggNOGFastas.mix(local_ref_eggnog_add).into {eggNOGFastas_concatenation; eggNOGFastas_mapping }
concatenated_references = eggNOGFastas_concatenation.collectFile(name: 'references.fasta', storeDir: params.reference_dir)

/*
  create a mapping file of reference sequence ID to EggNOG ID that can be
  used as synonyms file in megan
*/
process createEggNOGMap {
    input:
    file megan_eggnog_map from megan_eggnog_map.first()
    file fasta from eggNOGFastas_mapping


    tag "${fasta.simpleName}"
    // TODO dump eggnog_map or taxmap as python3 binary? to avoid parsing it again later
    // pickle!
    output:
    file "${fasta.baseName}_eggnog.map" into eggnog_map
    file "${fasta.baseName}_taxid.map" into tax_map
    file "${fasta.baseName}.class" into class_map

    script:
    template 'create_eggnog_map.py'
}

/*
  concatenate single mapping files and keep track of long/short IDs for EggNOG
*/
eggnog_map_concat = eggnog_map.collectFile(name: 'eggnog.syn', storeDir: params.reference_dir)
// tax_map_concat = tax_map.collectFile(name: 'tax.syn', storeDir: params.reference_dir)
class_map.collectFile(name: 'class.map', storeDir: params.reference_dir).into{class_map_concat; class_map_concat_copy}

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
    file "${fq.simpleName}.daa" into daa_files

    tag "${fq.simpleName}"
    //publishDir 'queries', mode: 'copy'

    script:
    """
    diamond blastx -q ${fq} --db references.dmnd -f 100 --unal 0 -e 1e-6 --out ${fq.simpleName}.daa --threads ${task.cpus}
    """
}


/*
  Prepare daa files for gc-assembler, using the daa-meganizer tool and
  the synonyms mapping file created earlier, applies to the daa files in place,
  so copy them
  TODO may cause severe memory overhead in the future as we copy files a lot
*/
process meganizeDAAFiles {

    input:
    file daa from daa_files
    file eggnog_map from eggnog_map_concat.first()

    output:
    file daa into daa_files_meganized

    stageInMode 'copy'
    tag "${daa.simpleName}"
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

    tag "${daa.simpleName}"
    publishDir "${params.queries_dir}/${daa.baseName}", mode: 'copy'

    script:
    """
    ${params.gc_assembler} -i ${daa} -fun EGGNOG -id ALL -mic 99 -v --minAvCoverage 2 -t ${task.cpus}
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
    file '*.faa' optional true into translated_contigs

    tag "${contig.simpleName}"
    publishDir "${params.queries_dir}/${contig.simpleName.minus(~/-.+/)}", mode: 'copy'

    script:
    template 'translate_dna_to_aa.py'
}


/*
  All assembled contigs are aligned to the corresponding EggNOG alignment.
  This is relatively fast using the mafft --addfragments option.
*/
process alignContigs {
    input:
    file faa from translated_contigs
    file reference_alignment from ref_alignments.toList()
    file class_map_concat from class_map_concat.first()

    output:
    file '*.aln' into aligned_contigs

    tag "${faa.baseName}"
    //publishDir 'queries', mode: 'copy'

    script:
    """
    id=\$(grep "^${faa.baseName.minus(~/^.+-/)}\\b" $class_map_concat | cut -f 2)".aln"
    mafft-fftnsi --adjustdirection --thread ${task.cpus} --addfragments ${faa} \$id > ${faa.baseName}.aln
    if [ ! -s ${faa.baseName}.aln ]
    then
      echo "the alignment file is empty, presumably mafft crashed and we retry...."
      rm ${faa.baseName}.aln
    fi
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
    // stageInMode 'copy'
    tag "${aln.baseName}"

    """
    #sed -i -E "/^>/! s/X/-/g" ${aln}
    sed -E "/^>/! s/X/-/g" ${aln} > ${aln.baseName}.clean
    trimal -in ${aln.baseName}.clean -out ${aln.baseName}.trim.aln -gappyout
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
    file "${aln.simpleName}.treefile" into trees

    publishDir "${params.queries_dir}/${aln.simpleName.minus(~/-.+/)}", mode: 'copy'
    tag "${aln.simpleName}"

    script:
    if (params.phylo_method == "iqtree")
      """
      #iqtree-omp -s ${aln} -m LG -bb 1000 -nt 2 -pre ${aln.simpleName}
      iqtree-1.6.beta4 -fast -s ${aln} -m LG -nt ${task.cpus} -pre ${aln.simpleName}
      """
    else if (params.phylo_method == "fasttree")
      """
      FastTree ${aln} > ${aln.simpleName}.treefile
      """
}

process magnetizeTrees {
    input:
    file tree from trees
    file '*' from tax_map.collect()
    file '*' from class_map_concat_copy.first()
    each lineage from lineage
    val rank from rank.first()

    output:
    file "${tree.baseName}.pdf" optional true into pdfs
    file 'decision.txt' into decisions
    stdout x
    tag "${tree.baseName} - $lineage"

    publishDir "${params.queries_dir}/${tree.simpleName.minus(~/-.+/)}", mode: 'copy'

    // beforeScript = {"ln -s \$(grep '^${tree.baseName.minus(~/^.+-/).minus(~/.trim/)}\\b' $workflow.launchDir/${params.reference_dir}/class.map | cut -f 2)\"_taxid.map\" tax.map"}

    script:
    template 'magnetize_tree.py'
}
x.subscribe{print it}

decisions_concat = decisions.collectFile(name: 'tree_decisions.txt', storeDir: params.queries_dir)

process decideSamples {
    input:
    file tree_decisions from decisions_concat

    output:
    file "sample_decisions.txt" into sample_decisions

    publishDir "${params.queries_dir}", mode: 'copy'

    script:
    template "decide_samples.py"
}

// code from J. Viklund of SciLifeLab Uppsala
def grab_git_revision() {
    if ( workflow.commitId ) { // it's run directly from github
        return workflow.commitId
    }

    // Try to find the revision directly from git
    head_pointer_file = file("${workflow.projectDir}/.git/HEAD")
    if ( ! head_pointer_file.exists() ) {
        return ''
    }
    ref = head_pointer_file.newReader().readLine().tokenize()[1]

    ref_file = file("${workflow.projectDir}/.git/$ref")
    if ( ! ref_file.exists() ) {
        return ''
    }
    revision = ref_file.newReader().readLine()

    return revision
}
