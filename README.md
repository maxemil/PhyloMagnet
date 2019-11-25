[![Docs Status](https://readthedocs.org/projects/phylomagnet/badge/?version=latest)](http://phylomagnet.readthedocs.io/en/latest/)
[![Build Status](https://travis-ci.org/maxemil/PhyloMagnet.svg?branch=master)](https://travis-ci.org/maxemil/PhyloMagnet)
[![Hosted](https://img.shields.io/badge/hosted-singularity--hub-blue.svg)](https://www.singularity-hub.org/collections/978)

# PhyloMagnet
## Pipeline for screening metagenomes, looking for arbitrary lineages, using gene-centric assembly methods and phylogenetics

## Abstract
Motivation: Metagenomic and metatranscriptomic sequencing analyses have become increasingly popular tools for producing massive amounts of short-read data, often used for the reconstruction of draft genomes or the detection of (active) genes in microbial communities. Unfortunately, sequence assemblies of such datasets generally remain a computationally challenging task. Frequently, researchers are only interested in a specific group of organisms or genes; yet, the assembly of multiple datasets only to identify candidate sequences for a specific question is sometimes prohibitively slow, forcing researchers to select a subset of available datasets to address their question. Here we present PhyloMagnet, a workflow to screen meta-omics datasets for taxa and genes of interest using gene-centric assembly and phylogenetic placement of sequences.
Results: Using PhyloMagnet, we could identify up to 87% of the genera in an in vitro mock community with variable abundances, while the false positive predictions per single gene tree ranged from 0% to 23%. When applied to a group of metagenomes for which a set of MAGs have been published, we could detect the majority of the taxonomic labels that the MAGs had been annotated with. In a metatranscriptomic setting the phylogenetic placement of assembled contigs corresponds to that of transcripts obtained from transcriptome assembly. See https://github.com/maxemil/PhyloMagnet-benchmarks for benchmark experiments.

## Quick installation & usage
For detailed documentation, please visit http://phylomagnet.readthedocs.io/en/latest/
```bash
# download the image with all tools installed using singularity 3x
singularity pull --name PhyloMagnet.sif shub://maxemil/PhyloMagnet:latest

# get versions of tools used in the pipeline:
singularity exec PhyloMagnet.sif conda list -n PhyloMagnet-<version>

# execute the test pipeline with nextflow
nextflow run main.nf \
          -with-singularity PhyloMagnet.sif \
          --is_runs true \
          --fastq "test/*rpoB.fastq.gz" \
          --reference_packages "test/rpkgs/*" \
          --lineage "order" \
          --megan_eggnog_map eggnog.map \
          --cpus 2 \
          --is_runs true \
          --queries_dir test/queries \
          --reference_dir test/references \
          --phylo_method 'fasttree' \
          --align_method 'mafft-fftnsi' \
          -w test/work -resume
```
## Citing PhyloMagnet

PhyloMagnet is published in Bioinformatics:
Max E Schön, Laura Eme, Thijs J G Ettema, PhyloMagnet: fast and accurate screening of short-read meta-omics data using gene-centric phylogenetics, Bioinformatics, btz799, https://doi.org/10.1093/bioinformatics/btz799

Please make sure to also cite all tools that are used in the pipeline if you use it for your research! Visit http://phylomagnet.readthedocs.io/en/latest/ or see the startup message for details.
