[![Docs Status](https://readthedocs.org/projects/phylomagnet/badge/?version=latest)](http://phylomagnet.readthedocs.io/en/latest/)
[![Build Status](https://travis-ci.org/maxemil/PhyloMagnet.svg?branch=master)](https://travis-ci.org/maxemil/PhyloMagnet)
[![Hosted](https://img.shields.io/badge/hosted-singularity--hub-blue.svg)](https://www.singularity-hub.org/collections/978)

# PhyloMagnet
## Pipeline for screening metagenomes, looking for arbitrary lineages, using gene-centric assembly methods and phylogenetics

## Abstract
Given a list of Bioproject IDs, PhyloMagnet downloads WGS Metagenome Illumina runs, aligns them against EggNOG references using diamond, then perform gene-centric-assembly and compute phylogenies for the gene clusters present in the references and query.

## Quick installation & usage
```bash
# download the image with all tools installed
singularity pull --name PhyloMagnet.simg shub://maxemil/PhyloMagnet:latest

# get versions of tools used in the pipeline:
./PhyloMagnet.simg

# execute the pipeline with nextflow
nextflow run main.nf --reference_classes eggnog.txt \
                     --project_list bioprojects.txt \
                     --database 'ena'
                     --phylo_method 'fasttree' \
                     --align_method 'mafft-einsi'
                     --queries_dir queries_output \
                     --reference_dir ref_output \
                     --megan_vmoptions 'MEGAN.vmoptions' \
                     --lineage 'Rickettsiales','Pelagibacterales' \
                     --cpus 20
```
## Author

* [Max Emil Schön](https://github.com/maxemil), PhD student @ [Ettemalab](https://www.ettemalab.org), Uppsala University
