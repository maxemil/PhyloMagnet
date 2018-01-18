[![Docs Status](https://readthedocs.org/projects/phylomagnet/badge/?version=latest)](phylomagnet.readthedocs.io)
# Pipeline for combing through metagenomes (PhyloMagnet)
## looking for arbitrary lineages, using gene-centric assembly methods

## (Abstract)
Given a list of Bioproject IDs, the pipeline should get WGS Metagenome Illumina runs, download them and align them against a reference using diamond, then perform gc-assembly and compute phylogenies for the gene clusters present in the references

## Installation & Usage
The documentation and usage examples are availabel at phylomagnet.readthedocs.io

## Future Plans
* try including RP15 pipeline, so producing linked contigs for genes (paired end information?)
  * Intergenic space length?
  * pull out reads, then mini-meta-assembly?
