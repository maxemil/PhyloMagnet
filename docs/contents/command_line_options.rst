Command line options
=============
To get a quick overview of all available command line options, run

.. code-block:: bash

  nextflow run maxemil/PhyloMagnet --help


Query option
-----
Fastq input (``--fastq``)
``````
Provide one or several short read samples (space separated list) in fastq format to PhyloMagnet. The files can be in gzip compressed form.

BioProject list (``--project_list``)
``````
Provide a list of BioProject identifiers (e.g. PRJNA324704) in a single file, one ID per line.

BioProject run IDs (``--is_runs``)
``````
Boolean value that modifies the ``--project_list`` option such that it expects run IDs instead of Project IDs (e.g. SRR3656745). Default is ``false``

Database (``--database``)
``````
Choose if you would like to download the sra data from ncbi's (``ncbi``) or ena's (``ena``) servers (depending on where you are, this might make a huge speed difference). Default is ``ena``

References otions
------
EggNOG reference Ids (``--reference_classes``)
``````
A file with identifiers from the EggNOG database, e.g. ``COG0051`` or ``ENOG410XPEN``. One ID per line.

Archived reference packages (``--reference_packages``)
``````
A single path to one or several compressed references packages (see the utility script ``make_reference_packages.sh`` in the ``utils`` folder) from a previous PhyloMagnet run. Can contain wildcards if put in double quotes. e.g. ``"my_rpkg/*.tgz"``

Local reference sequences (``--local_ref``)
``````
A single path to one or several local fasta files containing orthologous groups of proteins. the referece sequences should be annotated with their taxonomy ID in the NCBI taxonomy (e.g. ``562.NC_011750`` as the sequence record's header, ``562`` being the taxID).

Output options
--------
Output queries (``--queries_dir``)
``````
Output directory for queries; assembled contigs, placement results and summary tables/figures get saved here. Defaults to ``queries``.

Output references (``--reference_dir``)
``````
Output directory for references; fasta files, alignments, model files and trees get saved here. Defaults to ``references``.

Run parameters
----------
Phylogenetic method (``--phylo_method``)
``````
Phylogentic tool used for the reconstruction of the reference tree. Only used for references from EggNOG and local files, not packages. Accepted values: ``iqtree``, ``fasttree``, ``raxml``.

Alignment method (``--align_method``)
``````
Alignment tool used to compute reference alignments. Only used for references from EggNOG and local files, not packages. Accepted values: ``maftt-*``, ``prank``.

Taxonomic lineage (``--lineage``)
``````
the lineage(s) to report occurences for. Can be either a list of labels provided by the user (e.g. ``Rickettsiales,Holosporales``), a taxonmic rank (e.g. ``family``), or both (e.g. ``Rickettsia,family``)

No. of CPUs (``--cpus``)
``````
No. of CPUs to use. For more fine-grained uage of resources per process change the file ``nextflow.config``.
