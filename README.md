# Pipeline for combing through metagenomes (PhyloMagnet)
## looking for arbitrary lineages, using gene-centric assembly methods


## (Abstract)
Given a list of Bioproject IDs, the pipeline should get WGS Metagenome Illumina runs, download them and align them against a reference using diamond, then perform gc-assembly and compute phylogenies for the gene clusters present in the references

## Installation & Usage
* Install NextFlow to execute the workflow (www.nextflow.io)
```
curl -fsSL get.nextflow.io | bash
```

* [RECOMMENDED] Install Singularity to easily use all third-party software
```
#see further inforamtion on http://singularity.lbl.gov/install-linux
git clone https://github.com/singularityware/singularity.git
cd singularity
./autogen.sh
./configure --prefix=/usr/local
make
sudo make install
```
* If your third-party binaries are not in your path, specify their location in the command line (or change the configuration file locally)

* execute the workflow like this:
```
nextflow run main.nf --cpus 30 \
              --reference-classes EGGNOG_List \
              --project-list bioprojects.txt \
              --fastq local_file.fq \
              --phylo_method iqtree \
              --is_runs true \
              --local_ref local_references.fasta \
              --lineage 'Rickettsiales' \
              --taxonomy_level_trees 'order' \
              -with-singularity PhyloMagnet.img \
              -with-report \
              -resume
```
* I do most testing using Bioproject PRJNA104935 and a small list of EggNOG IDs, but you can use any BioProject that contains WGS Illumina (Meta)genome reads and any EggNOG Ids for which there are fasta and aln files available.
* this will create something like the following directory tree:
```
-rw-r--r-- 1 MaxEmil users      40 Apr 11 11:47 EGGNOG_List
-rw-r--r-- 1 MaxEmil users 1712353 Apr 11 11:57 eggnog.map
drwxr-xr-x 1 MaxEmil users     354 Apr 28 13:10 references
drwxr-xr-x 1 MaxEmil users      52 Apr 28 13:10 work
```

* change memory usage for daa-meganizer and gc-assembler with
```
sudo singularity exec --writable PhyloMagnet.img vim /usr/local/megan/MEGAN.vmoptions
```
## Future Plans
* try including RP15 pipeline, so producing linked contigs for genes (paired end information?)
  * Intergenic space length?
  * pull out reads, then mini-meta-assembly?
