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
    log.info "taxonomic level to be analysed    : $params.taxonomic_rank"
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

local_ref = optional_channel(params.local_ref)

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
  cp \$(readlink -f $ids) ${ids.simpleName}_local.txt
  for f in *.fasta;
  do
    ID=\$(comm -3 <(grep -v '^-' $megan_eggnog_map | cut -f2 | cut -d ' ' -f 1 | sort ) <(cat ${ids.simpleName}_local.txt | sort) | tr -d '\\t' | head -n 1)
    cp -L \$f \$ID.\${f%%.*}.fasta
    echo \$ID >> ${ids.simpleName}_local.txt
    echo \${f%%.*}"\t"\$ID >> local_translation.txt
  done
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
    // file "${id}.aln" into eggNOGAlignments

    publishDir params.reference_dir, mode: 'copy', overwrite: false
    tag "$id"

    """
    wget http://eggnogapi.embl.de/nog_data/text/fasta/${id} -O - | gunzip > ${id}.fasta || wget http://eggnogapi.embl.de/nog_data/text/fasta/${id} -O ${id}.fasta
    """
}

references_fastas = eggNOGFastas.mix(local_ref_eggnog)

/*
  create a mapping file of reference sequence ID to EggNOG ID that can be
  used as synonyms file in megan
*/
process createMappingFiles {
    input:
    file megan_eggnog_map from megan_eggnog_map.first()
    file fasta from references_fastas
    file local_ref_translation from local_ref_translation

    output:
    file "*.eggnog.map" into eggnog_map
    file "*.taxid.map" into tax_map
    file "*.class" into class_map
    file "*.unique.fasta" into references_unique_fastas

    tag "${fasta.baseName}"

    script:
    template 'create_eggnog_map.py'
}
references_unique_fastas.into{references_unique_fastas_align; references_unique_fastas_concat}

process alignReferences {
  input:
  file fasta from references_unique_fastas_align

  output:
  file "${fasta.baseName}.aln" into ref_alignments

  publishDir params.reference_dir, mode: 'copy'
  tag "${fasta.baseName}"
  stageInMode 'copy'

  script:
  if (params.align_method == "prank")
    """
    sed -i 's/\\*//g' $fasta
    prank -protein -d=$fasta -o=${fasta.baseName} -f=fasta
    mv ${fasta.baseName}.best.fas ${fasta.baseName}.aln
    """
  else if (params.align_method.startsWith("mafft"))
    """
    sed -i 's/\\*//g' $fasta
    $params.align_method --adjustdirection --thread ${task.cpus}  $fasta > ${fasta.baseName}.aln
    """
}


/*
  Concatenate all fastA into one single references.fasta for diamond alignment
*/
concatenated_references = references_unique_fastas_concat.collectFile(name: 'references.fasta', storeDir: params.reference_dir)

/*
  concatenate single mapping files and keep track of long/short IDs for EggNOG
*/
eggnog_map_concat = eggnog_map.collectFile(name: 'eggnog.syn', storeDir: params.reference_dir)
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
    ${params.daa_meganizer} --in ${daa} -s2eggnog ${eggnog_map}
    """
}


/*
  Perform Gene Centric Assembly for all classes of EggNOG. this will create
  fastA files of the form <runID>-<shortEggNOGID>.fasta
*/
process geneCentricAssembly {
    input:
    file daa from daa_files_meganized
    file class_map from class_map_concat.first()

    output:
    file '*.fasta' into assembled_contigs mode flatten

    tag "${daa.simpleName}"
    publishDir "${params.queries_dir}/${daa.baseName}", mode: 'copy'

    script:
    """
    ${params.gc_assembler} -i ${daa} -fun EGGNOG -id ALL -mic 99 -vv --minAvCoverage 2 -t ${task.cpus}
    while read id name ; do
      if [ -e ${daa.simpleName}-\$id.fasta ];
      then
        mv ${daa.simpleName}-\$id.fasta ${daa.simpleName}-\$name.fasta
      fi
    done < $class_map
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
    file contig from assembled_contigs

    output:
    file '*.faa' optional true into translated_contigs

    tag "${contig.simpleName}"
    publishDir "${params.queries_dir}/${contig.simpleName.tokenize('-')[0]}", mode: 'copy'

    script:
    template 'translate_dna_to_aa.py'
}

process buildTreefromReferences {
  input:
  file reference_alignment from ref_alignments

  output:
  set file("${reference_alignment.simpleName}.treefile"), file("${reference_alignment.simpleName}.modelfile"), file("$reference_alignment") into reference_trees

  publishDir "${params.reference_dir}", mode: 'copy'
  tag "${reference_alignment.simpleName}"

  script:
  if (params.phylo_method.startsWith("iqtree"))
    """
    ${params.phylo_method} -s ${reference_alignment} -m LG+G+F -nt AUTO -ntmax ${task.cpus} -pre ${reference_alignment.simpleName}
    raxml -f e -s ${reference_alignment} -t ${reference_alignment.simpleName}.treefile -T ${task.cpus} -n file -m PROTGAMMALGF
    mv RAxML_info.file ${reference_alignment.simpleName}.modelfile
    """
  else if (params.phylo_method == "fasttree")
    """
    FastTree -log ${reference_alignment.simpleName}.log -lg ${reference_alignment} > ${reference_alignment.simpleName}.treefile
    raxml -f e -s ${reference_alignment} -t ${reference_alignment.simpleName}.treefile -T ${task.cpus} -n file -m PROTGAMMALGF
    mv RAxML_info.file ${reference_alignment.simpleName}.modelfile
    """
  else if (params.phylo_method == "raxml")
    """
    raxml -f e -s ${reference_alignment} -t ${reference_alignment.simpleName}.treefile -T ${task.cpus} -n file -m PROTGAMMALGF
    """
}

process alignQueriestoRefMSA {
  input:
  set file(reftree), file(modelinfo), file(refalignment) from reference_trees
  each contigs from translated_contigs

  output:
  set file("$reftree"), file("$modelinfo"), file("$refalignment"), file("${contigs.simpleName}.ref.aln") into aligned_queries

  when:
  "${refalignment.simpleName}" == "${contigs.simpleName.tokenize('-')[1]}"

  tag "${contigs.simpleName} - ${refalignment.simpleName}"

  script:
  """
  trimal -in $refalignment -out ${refalignment.simpleName}.phy -phylip
  papara -t $reftree -s ${refalignment.simpleName}.phy -q $contigs -a -n ${contigs.simpleName} -r
  trimal -in papara_alignment.${contigs.simpleName} -out ${contigs.simpleName}.ref.aln -fasta
  """
}

process placeContigsOnRefTree {
  input:
  set file(reftree), file(modelinfo), file(refalignment), file(queryalignment) from aligned_queries

  output:
  file "${queryalignment.simpleName}.jplace" into placed_contigs

  tag "${queryalignment.simpleName} - ${refalignment.simpleName}"

  script:
  """
  epa-ng --ref-msa $refalignment --tree $reftree --query $queryalignment --model $modelinfo --no-heur
  mv epa_result.jplace ${queryalignment.simpleName}.jplace
  """
}

process assignContigs {
  input:
  file placed_contigs from placed_contigs
  each file(tax_map) from tax_map

  output:
  file "${placed_contigs.simpleName}.csv" into profiles
  file "${placed_contigs.simpleName}.assign" into assignments
  file "${placed_contigs.simpleName}.svg" into colored_tree_svg
  file "${placed_contigs.simpleName}.newick" into placement_tree

  publishDir "${params.queries_dir}/${placed_contigs.baseName.tokenize('-')[0]}", mode: 'copy'

  when:
  "${tax_map.simpleName}" == "${placed_contigs.baseName.tokenize('-')[1]}"

  tag "${placed_contigs.simpleName}"

  script:
  """
  gappa analyze assign --jplace-path $placed_contigs --taxon-file $tax_map
  mv profile.csv ${placed_contigs.simpleName}.csv
  mv per_pquery_assign ${placed_contigs.simpleName}.assign
  gappa analyze visualize-color --jplace-path $placed_contigs --write-svg-tree
  mv tree.svg ${placed_contigs.simpleName}.svg
  gappa analyze graft --name-prefix 'Q_' --jplace-path $placed_contigs
  """
}


process magnetizeTrees {
    input:
    file profile from profiles
    each lineage from lineage
    val rank from rank.first()

    output:
    file 'decision.txt' into decisions
    stdout x

    tag "${profile.baseName} - $lineage"
    publishDir "${params.queries_dir}/${profile.simpleName.tokenize('-')[0]}", mode: 'copy'

    script:
    template 'magnetize_tree.py'
}
x.subscribe{print it}

decisions_concat = decisions.collectFile(name: 'tree_decisions.txt', storeDir: params.queries_dir)

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


def optional_channel(argument) {
  if(argument != ""){
    handle = file( argument )
    if (argument && handle.exists()) {
      return Channel.fromPath(argument)
    }else if (argument) {
      return Channel.from(argument)
    }else{
      return Channel.empty()
    }
  }else{
    return Channel.empty()
  }
}
