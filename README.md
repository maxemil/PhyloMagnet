[![Docs Status](https://readthedocs.org/projects/phylomagnet/badge/?version=latest)](http://phylomagnet.readthedocs.io/en/latest/)
[![Build Status](https://travis-ci.org/maxemil/PhyloMagnet.svg?branch=master)](https://travis-ci.org/maxemil/PhyloMagnet)

# PhyloMagnet
## Pipeline for screening metagenomes, looking for arbitrary lineages, using gene-centric assembly methods and phylogenetics

## Abstract
Given a list of Bioproject IDs, the pipeline should get WGS Metagenome Illumina runs, download them and align them against a reference using diamond, then perform gc-assembly and compute phylogenies for the gene clusters present in the references

## Installation & Usage
The documentation and usage examples are availabel at http://phylomagnet.readthedocs.io/en/latest/

### Future Plans
* try including RP15 pipeline, so producing linked contigs for genes (paired end information?)
  * Intergenic space length?
  * pull out reads, then mini-meta-assembly?
