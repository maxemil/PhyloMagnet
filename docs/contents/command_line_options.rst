Command line options
=============

Query option
-----
Fastq input (``--fastq``)
``````
Provide one or several short read samples (space separated list) in fastq format to PhyloMagnet.

BioProject list (``--project_list``)
``````
Provide a list of BioProject identifiers (e.g. PRJNA324704) in a single file, one ID per line.

BioProject run IDs (``--is_runs``)
``````
Boolean value that modifies the ``--project_list`` option such that it expects run IDs instead of Project IDs (e.g. SRR3656745). Default is ``false``

Database (``--database``)
``````
Choose if you would like to download the sra data from ncbi's (``ncbi``) or ena's (``ena``) servers (depending on where you are, this might make a huge difference). Default is ``ena``

References otions
------

Run parameters
----------
