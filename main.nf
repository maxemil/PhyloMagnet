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
    log.info "Use a local fastq file            : $params.fastq"
    log.info "Use run IDs instead of projects   : $params.is_runs"
    log.info "Lineage(s) to look for            : $params.lineage"
    log.info "Sequence DB to download from      : $params.database"
    log.info ""
    log.info "Binaries location (use default if singularity image is used)"
    log.info "location of MEGAN6                : $params.megan_dir"
    log.info "Python 3 binary used              : $params.python3"
    log.info ""
    log.info ""
}

startup_message()

// reads a list of Bioproject IDs, but testing only on one single ID

ids = optional_channel_from_file(params.project_list)

lineage_list = params.lineage.tokenize(',')
lineage = Channel.from(lineage_list)
// reads a list of Bioproject IDs, but testi


eggNOGIDs = optional_channel_from_file(params.reference_classes)
eggNOGIDs_local = optional_channel_from_path(params.reference_classes)

optional_channel_from_path(params.local_ref).into {local_ref_include_raw; local_ref_map; local_ref_save}
fastq_files = optional_channel_from_path(params.fastq)
local_ref_save.subscribe{ it.copyTo("${params.reference_dir}/${it.simpleName}/${it.simpleName}.fasta") }
local_ref_include = local_ref_include_raw.ifEmpty("execute anyway")

Channel.from(file(params.megan_vmoptions)).into { megan_vmoptions_meganizer; megan_vmoptions_assembler }

ref_packages = optional_channel_from_path(params.reference_packages)

warning_no_queries = "WARNING: no query sequences or IDs given, but will continue to prepare references"
if (!params.fastq && !params.project_list) {
  print(warning_no_queries)
}


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
    file 'runs.txt' into project_runs

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

all_runs = project_runs.flatMap {it -> it.readLines() }


/*
  for all valid runs specified in 'runs.txt', download fastQ files (in parallel).
  split reads if paired end data, but keep all in one file for diamond
  alignment and assembly
*/
process downloadFastQ {
    input:
    val run_id from all_runs

    output:
    file "*.fastq.gz" into fastq_files_SRA mode flatten

    publishDir params.queries_dir
    tag "$run_id"

    script:
    template 'download_fastq.py'
}

fastq_files_all = fastq_files.mix(fastq_files_SRA)

process unpackRefPackage {
  input:
  file rpkg from ref_packages

  output:
  file "${rpkg.simpleName}.fasta" into local_ref_rpkg
  file "${rpkg.simpleName}.unique.fasta" into local_ref_unique_rpkg
  set file("${rpkg.simpleName}.treefile"), file("${rpkg.simpleName}.modelfile"), file("${rpkg.simpleName}.unique.aln") into reference_trees_rpkg
  file "${rpkg.simpleName}.taxid.map" into tax_map_rpkg
  file "${rpkg.simpleName}.eggnog.map" into eggnog_map_rpkg
  file "${rpkg.simpleName}.class" into class_rpkg

  publishDir "${params.reference_dir}/${rpkg.simpleName}", mode: 'copy'
  tag "${rpkg.simpleName}"

  script:
  """
  tar xzfv $rpkg
  """
}

local_ref_include_w_rpkg = local_ref_include.mix(local_ref_rpkg)
class_rpkg.into{class_rpkg_map_raw; class_rpkg_concat}

class_rpkg_map = class_rpkg_map_raw.ifEmpty("execute anyway")

process includeLocalRef {
  input:
  file "*" from local_ref_include_w_rpkg.collect()
  file "*" from class_rpkg_map.collect()

  output:
  file "eggnog.map" into extended_eggnog_map
  file "data.jar" into data_jar

  publishDir params.reference_dir, pattern: '*.jar'

  script:
  """
  cp -L $params.megan_dir/jars/data.jar .
  jar -xf data.jar
  mv resources/files/eggnog.map .

  if [ ! -n "\$(echo *.class)" ];
  then
    for c in *.class;
    do
      if [ ! \$(grep "\${c%%.*}\\s" eggnog.map) ];
      then
        cat \$c >> eggnog.map
      fi
    done
  fi

  for f in *.fasta;
  do
    if [ ! \$(grep "\${f%%.*}\\s" eggnog.map) ];
    then
      ID="\$((\$(tail -n 1 eggnog.map | cut -f1) + 1))"
      printf "%i\\t%s\\n" "\$ID" "\${f%%.*}" >> eggnog.map
    fi
  done

  cp eggnog.map resources/files/
  jar -cvf data.jar resources/*
  rm -rf resources
  """
}

data_jar.into{data_jar_meganizer; data_jar_assembler}

/*
  Download all raw sequence files as well as the untrimmed alignment
  files for all provided EggNOG Ids.
*/
process downloadEggNOG {
    input:
    val id from eggNOGIDs

    output:
    file "${id}.fasta" into eggNOGFastas

    publishDir "${params.reference_dir}/$id", mode: 'copy'
    tag "$id"

    """
    wget http://eggnogapi.embl.de/nog_data/text/fasta/${id} -O - | gunzip > ${id}.fasta || wget http://eggnogapi.embl.de/nog_data/text/fasta/${id} -O ${id}.fasta
    """
}

references_fastas = eggNOGFastas.mix(local_ref_map)

/*
  create a mapping file of reference sequence ID to EggNOG ID that can be
  used as synonyms file in megan
*/
process createMappingFiles {
    input:
    file megan_eggnog_map from extended_eggnog_map
    file fasta from references_fastas

    output:
    file "*.eggnog.map" into eggnog_map
    file "*.taxid.map" into tax_map
    file "*.class" into class_map
    file "*.unique.fasta" into references_unique_fastas

    publishDir "${params.reference_dir}/${fasta.simpleName}", mode: 'copy'
    tag "${fasta.baseName}"

    script:
    template 'create_eggnog_map.py'
}

tax_map_combined = tax_map.mix(tax_map_rpkg)

references_unique_fastas.into{references_unique_fastas_align; references_unique_fastas_concat_raw}

references_unique_fastas_concat = references_unique_fastas_concat_raw.mix(local_ref_unique_rpkg)

/*
  Concatenate all fastA into one single references.fasta for diamond alignment
*/
concatenated_references = references_unique_fastas_concat.collectFile(name: 'references.fasta', storeDir: params.reference_dir)

/*
  concatenate single mapping files and keep track of long/short IDs for EggNOG
*/
eggnog_map_concat = eggnog_map.mix(eggnog_map_rpkg).collectFile(name: 'eggnog.syn', storeDir: params.reference_dir)
class_map_concat = class_map.mix(class_rpkg_concat).collectFile(name: 'class.map', storeDir: params.reference_dir)

/*
  prepare the diamond database from the concatenated fastA file
*/
process diamondMakeDB {
    input:
    file 'references.fasta' from concatenated_references

    output:
    file 'references.dmnd' into diamond_database

    publishDir params.reference_dir, mode: 'copy'

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
    file fq from fastq_files_all
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
    file "data.jar" from data_jar_meganizer
    file "MEGAN.vmoptions" from megan_vmoptions_meganizer

    output:
    file daa into daa_files_meganized

    stageInMode 'copy'
    tag "${daa.simpleName}"
    //publishDir 'queries', mode: 'copy', overwrite: true

    script:
    """
    vmOptions=\$(grep '^-' MEGAN.vmoptions | tr "\\n" " ")
    java \$vmOptions -cp "$params.megan_dir/jars/MEGAN.jar:data.jar" \
          megan.tools.DAAMeganizer --in ${daa} -s2eggnog ${eggnog_map}
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
    file "data.jar" from data_jar_assembler
    file "MEGAN.vmoptions" from megan_vmoptions_assembler

    output:
    file '*.fasta' into assembled_contigs mode flatten

    tag "${daa.simpleName}"
    publishDir "${params.queries_dir}/${daa.baseName}", mode: 'copy'

    script:
    """
    vmOptions=\$(grep '^-' MEGAN.vmoptions | tr "\\n" " ")
    java \$vmOptions -cp "$params.megan_dir/jars/MEGAN.jar:data.jar" \
          megan.tools.GCAssembler -i ${daa} -fun EGGNOG -id ALL -mic 99 -vv --minAvCoverage 2 -t ${task.cpus}

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


process alignReferences {
  input:
  file fasta from references_unique_fastas_align

  output:
  file "${fasta.baseName}.aln" into ref_alignments

  publishDir "${params.reference_dir}/${fasta.simpleName}", mode: 'copy'
  tag "${fasta.simpleName}"
  stageInMode 'copy'

  script:
  if (params.align_method == "prank")
    """
    prank -protein -d=$fasta -o=${fasta.baseName} -f=fasta
    sed -i '/^>/! s/[U|*|X]/-/g' ${fasta.baseName}.best.fas
    trimal -in ${fasta.baseName}.best.fas -out ${fasta.baseName}.aln -gt 0 -fasta
    """
  else if (params.align_method.startsWith("mafft"))
    """
    $params.align_method --adjustdirection --anysymbol --thread ${task.cpus}  $fasta > ${fasta.baseName}.mafft.aln
    sed -i '/^>/! s/[U|*|X]/-/g' ${fasta.baseName}.mafft.aln
    trimal -in ${fasta.baseName}.mafft.aln -out ${fasta.baseName}.aln -gt 0 -fasta
    """
}


process buildTreefromReferences {
  input:
  file reference_alignment from ref_alignments

  output:
  set file("${reference_alignment.simpleName}.treefile"), file("${reference_alignment.simpleName}.modelfile"), file("$reference_alignment") into reference_trees

  publishDir "${params.reference_dir}/${reference_alignment.simpleName}", mode: 'copy'
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

reference_trees_combined = reference_trees.mix(reference_trees_rpkg)

process alignQueriestoRefMSA {
  input:
  set file(reftree), file(modelinfo), file(refalignment) from reference_trees_combined
  each contigs from translated_contigs

  output:
  set file("$reftree"), file("$modelinfo"), file("$refalignment"), file("${contigs.simpleName}.ref.aln") into aligned_queries

  tag "${contigs.simpleName} - ${refalignment.simpleName}"
  publishDir "${params.queries_dir}/${contigs.simpleName.tokenize('-')[0]}", mode: 'copy'

  when:
  "${refalignment.simpleName}" == "${contigs.simpleName.tokenize('-')[1]}"

  script:
  """
  trimal -in $refalignment -out ${refalignment.simpleName}.phy -phylip
  papara -t $reftree -s ${refalignment.simpleName}.phy -q $contigs -a -n ${contigs.simpleName} -r
  trimal -in papara_alignment.${contigs.simpleName} -out ${contigs.simpleName}.ref.aln -fasta
  """
}


process splitAlignmentsRefQuery {
  input:
  set file(reftree), file(modelinfo), file(refalignment), file(queryalignment) from aligned_queries

  output:
  set file("$reftree"), file("$modelinfo"), file("$refalignment"), file("${queryalignment.simpleName}.queries.aln") into aligned_queries_split

  tag "${queryalignment.simpleName} - ${refalignment.simpleName}"

  script:
  template 'split_alignments_ref_query.py'
}


process placeContigsOnRefTree {
  input:
  set file(reftree), file(modelinfo), file(refalignment), file(queryalignment) from aligned_queries_split

  output:
  file "${queryalignment.simpleName}.jplace" into placed_contigs

  tag "${queryalignment.simpleName} - ${refalignment.simpleName}"
  publishDir "${params.queries_dir}/${queryalignment.simpleName.tokenize('-')[0]}", mode: 'copy'

  script:
  """
  epa-ng --ref-msa $refalignment --tree $reftree --query $queryalignment --model $modelinfo --no-heur
  mv epa_result.jplace ${queryalignment.simpleName}.jplace
  """
}


process assignContigs {
  input:
  file placed_contigs from placed_contigs
  each file(tax_map) from tax_map_combined

  output:
  set file("${placed_contigs.simpleName}.csv"), file("$tax_map") into profiles
  file "${placed_contigs.simpleName}.assign" into assignments
  file "${placed_contigs.simpleName}.svg" into colored_tree_svg
  file "${placed_contigs.simpleName}.newick" into placement_tree

  publishDir "${params.queries_dir}/${placed_contigs.baseName.tokenize('-')[0]}", mode: 'copy'
  tag "${placed_contigs.simpleName}"

  when:
  "${tax_map.simpleName}" == "${placed_contigs.baseName.tokenize('-')[1]}"


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
    set file(profile), file(tax_map) from profiles
    each lineage from lineage

    output:
    file "decision_${profile.baseName.tokenize('-')[1]}.txt" into decisions
    stdout x

    tag "${profile.baseName} - $lineage"
    publishDir "${params.queries_dir}/${profile.simpleName.tokenize('-')[0]}", mode: 'copy', overwrite: 'true'

    script:
    template 'magnetize_tree.py'
}
x.subscribe{print it}

boolean fileSuccessfullyDeleted = new File("${params.queries_dir}/tree_decisions.txt").delete()
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


def optional_channel_from_path(argument) {
  if(argument != ""){
    return Channel.fromPath(argument)
  }else{
    return Channel.empty()
  }
}


def optional_channel_from_file(argument) {
  if(argument != ""){
    handle = file( argument )
    if ( handle.exists() ){
      return  Channel.from(handle.readLines())
    }else{
      return Channel.empty()
    }
  }else{
    return Channel.empty()
  }
}
