Example Usage
=============
You can use PhyloMagnet in a number of different scenarios:

Using Bioproject IDs and eggNOG references
------------------------------------------

.. code-block:: bash

   nextflow run maxemil/PhyloMagnet --reference_classes eggnog.txt \
                         --project_list bioprojects.txt \
                         --phylo_method fasttree \
                         --queries_dir queries_output \
                         --reference_dir ref_output \
                         --megan_eggnog_map eggnog_map.txt \
                         --lineage Rickettsiales \
                         --taxonomy_level_trees order \
                         --cpus 20



Using SRA run IDs and local references
--------------------------------------

.. code-block:: bash

   nextflow run maxemil/PhyloMagnet --local_ref customOG.fasta \
                         --project_list sra_runs.txt \
                         --is_runs true \
                         --phylo_method fasttree \
                         --queries_dir queries_output \
                         --reference_dir ref_output \
                         --megan_eggnog_map eggnog_map.txt \
                         --lineage Rickettsiales \
                         --taxonomy_level_trees order \
                         --cpus 20



Using local FastQ file and local + eggNOG references
----------------------------------------------------

.. code-block:: bash

   nextflow run maxemil/PhyloMagnet --reference_classes eggnog.txt \
                         --local_ref customOG.fasta \
                         --phylo_method fasttree \
                         --fastq local_metagenome.fastq.gz
                         --queries_dir queries_output \
                         --reference_dir ref_output \
                         --megan_eggnog_map eggnog_map.txt \
                         --lineage Rickettsiales \
                         --taxonomy_level_trees order \
                         --cpus 20
