Command line options
=============
To get a quick overview of all available command line options, run

.. code-block:: bash

  nextflow run maxemil/PhyloMagnet --help


Query option
-----
Fastq input (``--fastq``)
``````
Provide one or several short read samples (with wildcards in double quotes) in
fastq format to PhyloMagnet. The files can be in gzip compressed form.
E.g. ``"fastq/*.fastq.gz"``.


BioProject list (``--project_list``)
``````
Provide a list of BioProject identifiers (e.g. ``PRJNA324704``) in a single
file, one ID per line.


BioProject run IDs (``--is_runs``)
``````
Boolean value that modifies the ``--project_list`` option such that it expects
run IDs instead of Project IDs (e.g. ``SRR3656745``). Default is ``false``


Database (``--database``)
``````
Choose if you would like to download the sra data from ncbi's (``ncbi``) or
ena's (``ena``) servers (depending on where you are, this might make a huge
speed difference). Default is ``ena``


References otions
------
EggNOG reference Ids (``--reference_classes``)
``````
A file with identifiers from the EggNOG database, e.g. ``COG0051`` or
``ENOG410XPEN``. One ID per line.


Archived reference packages (``--reference_packages``)
``````
A single path to one or several compressed references packages (see the utility
script ``make_reference_packages.sh`` in the ``utils`` folder) from a previous
PhyloMagnet run. Can contain wildcards if put in double quotes. e.g.
``"my_rpkg/*.tgz"``


Local reference sequences (``--local_ref``)
``````
A single path to one or several local fasta files containing orthologous groups
of proteins. the referece sequences should be annotated with their taxonomy ID
in the NCBI taxonomy (e.g. ``562.NC_011750`` as the sequence record's header,
``562`` being the taxID).


Output options
--------
Output queries (``--queries_dir``)
``````
Output directory for queries; assembled contigs, placement results and summary
tables/figures get saved here. Defaults to ``queries``.


Output references (``--reference_dir``)
``````
Output directory for references; fasta files, alignments, model files and trees
get saved here. Defaults to ``references``.


Run parameters
----------
Phylogenetic method (``--phylo_method``)
``````
Phylogentic tool used for the reconstruction of the reference tree. Only used
for references from EggNOG and local files, not packages. Accepted values:
``iqtree``, ``fasttree``, ``raxml``.


Alignment method (``--align_method``)
``````
Alignment tool used to compute reference alignments. Only used for references
from EggNOG and local files, not packages.
Accepted values: ``maftt-*``, ``prank``.


Taxonomic lineage (``--lineage``)
``````
the lineage(s) to report occurences for. Can be either a list of labels provided
by the user (e.g. ``Rickettsiales,Holosporales``), a taxonmic rank (e.g.
``family``), or both (e.g. ``Rickettsia,family``)


No. of CPUs (``--cpus``)
``````
No. of CPUs to use. For more fine-grained uage of resources per process change
the file ``nextflow.config``.


Threshold for plotting taxonomic labels (``--plot_threshold``)
``````
Threshold for filtering low frequent taxon labels from summary plots. e.g. label
would not get plotted when present in only 1 out of 4 trees for default value.
Threshold is checked per sample. Default: 0.25, accepted values between o and 1


Threshold for including taxonomic labels (``--aLWR_threshold``)
``````
Threshold of accumulated likelihood weight ratio (aLWR, see ``gappa``'s
documentation) to include labels in the summary table. Default 0.8, accepted
values between o and 1


MEGAN VM options file (``--megan_vmoptions``)
``````
File with options that are passed on to the Java virtual machine running MEGAN.
Most importantly, state here the amount of memory that is available, e.g.
`-Xmx16G` for 16GB. By default the file is expected to be in the execution
directory. Example file `MEGAN.vmoptions` is included in the repository.


Location of MEGAN (``--megan_dir``)
``````
The directory MEGAN's source files are located. When Manually installing MEGAN
it could be something like `/usr/local/megan`. Leave as default if the
singularity image is used.


Location of Python3 (``--python3``)
``````
Location of the python3 executable that has all needed packages available.
Should usually be `/usr/bin/env python3`, leave as default is using the
singularity image.
