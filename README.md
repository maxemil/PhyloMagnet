[![Docs Status](https://readthedocs.org/projects/phylomagnet/badge/?version=latest)](http://phylomagnet.readthedocs.io/en/latest/)
[![Build Status](https://travis-ci.org/maxemil/PhyloMagnet.svg?branch=master)](https://travis-ci.org/maxemil/PhyloMagnet)
[![Hosted](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://www.singularity-hub.org/collections/482)

# PhyloMagnet
## Pipeline for screening metagenomes, looking for arbitrary lineages, using gene-centric assembly methods and phylogenetics

## Abstract
Given a list of Bioproject IDs, PhyloMagnet downloads WGS Metagenome Illumina runs, aligns them against EggNOG references using diamond, then perform gene-centric-assembly and compute phylogenies for the gene clusters present in the references and query.

## Quick installation & usage
```bash
# download the image with all tools installed
singularity pull --name PhyloMagnet.simg shub://maxemil/PhyloMagnet:latest

# execute the pipeline with nextflow
nextflow run main.nf --reference_classes eggnog.txt \
                     --project_list bioprojects.txt \
                     --phylo_method fasttree \
                     --queries_dir queries_output \
                     --reference_dir ref_output \
                     --megan_eggnog_map eggnog_map.txt \
                     --lineage Rickettsiales \
                     --taxonomy_level_trees order \
                     --cpus 20
```

## Author

* [Max Emil Schön](https://github.com/maxemil), PhD student @ [Ettemalab](https://www.ettemalab.org), Uppsala University
