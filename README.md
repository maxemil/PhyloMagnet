# Pipeline for combing through metagenomes (PhyloMagnet)
## looking for arbitrary lineages, using gene-centric assembly methods


## (Abstract)
Given a list of Bioproject IDs, the pipeline should get WGS Metagenome Illumina runs, download them and align them against a reference using diamond, then perform gc-assembly and compute phylogenies for the gene clusters present in the references

## Installation & Usage
* Install NextFlow to execute the workflow (www.nextflow.io)
```
curl -fsSL get.nextflow.io | bash
```
* If your third-party binaries are not in your path, specify their location in the command line (or change the script locally)
* execute the workflow like this:
```
nextflow run PhyloComb.nf --cpus 30 --reference-classes EGGNOG_List --project-list bioprojects.txt
```
* I do most testing using Bioproject PRJNA104935 and a small list of EggNOG IDs, but you can use any BioProject that contains WGS Illumina (Meta)genome reads and any EggNOG Ids for which there are fasta and aln files available.
* this will create something like the following directory tree:
```
-rw-r--r-- 1 MaxEmil users      40 Apr 11 11:47 EGGNOG_List
-rw-r--r-- 1 MaxEmil users 1712353 Apr 11 11:57 eggnog.map
drwxr-xr-x 1 MaxEmil users     354 Apr 28 13:10 references
drwxr-xr-x 1 MaxEmil users      52 Apr 28 13:10 work

```
## Development PLans
* try including RP15 pipeline, so producing linked contigs for genes (paired end information?)
* Intergenic space?
* Refine interface to be able to use own fastq files
* Try using pplacer or similar instead of aligning and redoing phylogenies for each marker multiple times....

## comparison to PhyloSift
* HMM search (might be more sensitive...)
* FastTree for phylogenetics
* No internal assembly, so either assembled metagnomes or raw reads  
* Makes a whole community analyses to see what is there, so much slower but a little bit more general

## Log

### 17-04-03
* Rewrote the start of the current pipeline using sciluigi as a workflow manager in python3

### 17-04-07
* Check PhyloSift, but using assembled genes should be better than raw reads for taxonomy
* Maybe also generate Krona plots?
* provide synonyms file to daa-meganizer in the form `<sequenceID>\t<EGGNOGID>`
* where the EGGNOG ID is found in the .map file provided by megan (COG0085 becomes 85)
* and the sequenceID is actually a part of the sequence id in the fasta file, so 1000565.METUNv1_02263 is `<TaxID>.<sequenceID>`
* It should be possible to simply translate output of gc-assembly into proteins, but check for right codon table. (i.e. no stop codons in the middle of the sequence)


### 17-04-13
* Put code on bibucket
* download reference alignments and trees from EGGNOG, then use --addfragments in mafft to add contigs.

### 17-04-18
* Download aligned references for all EggNOG classes, then only add assembled contigs

### 17-04-20
* changed language (again...) to NextFlow

### 17-04-24
* NextFlow implementation ready up to alignment and trimming of assembled contigs.

### 17-04-28
* Put third-party binary paths in the params, so you can change them easily on the CL
