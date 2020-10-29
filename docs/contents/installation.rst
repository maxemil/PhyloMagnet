Installation Instructions for PhyloMagnet
=========================================
Nextflow
--------
To run PhyloMagnet, you will need to install `Nextflow <https://www.nextflow.io/>`_, a pipeline execution framework. This is however quite simple:

.. code-block:: bash

  curl -s https://get.nextflow.io | bash


Singularity
-----------
We provide a singularity container with all necessary tools installed and configured. To use it, you first need to install `Singularity` (either `version 2 <http://singularity.lbl.gov/install-linux>`_ or `version 3 <https://www.sylabs.io/guides/3.2/user-guide/quick_start.html#quick-installation-steps>`_) itself:

.. code-block:: bash

  VERSION=2.6
  wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz
  tar xvf singularity-$VERSION.tar.gz
  cd singularity-$VERSION
  ./configure --prefix=/usr/local
  make
  sudo make install

  # for installation of version 3 see the Singularity website


Singularity container
----------------
Then, download the container from singularity-hub or build it locally with the singularity recipe:

.. code-block:: bash

  # singularity 2

  singularity pull --name PhyloMagnet.simg shub://maxemil/PhyloMagnet:latest
  # or
  sudo singularity build PhyloMagnet.simg Singularity

  # singularity 3

  singularity pull --name PhyloMagnet.sif shub://maxemil/PhyloMagnet:latest
  # or
  sudo singularity build PhyloMagnet.sif Singularity

.. note::

  We highly recommend using the provided Singularity container to install all needed Software. Installing everything directly on the machine can be achieved using the conda environment.yml file.

If you want to see the version of the tools installed in the container, simply use conda to list all installed packages:

.. code-block:: bash

  singularity exec PhyloMagnet.sif conda list -n PhyloMagnet-<version>


PhyloMagnet
--------

now you can either get PhyloMagnet from github (clone or download from github.com/maxemil/PhyloMagnet) or let Nextflow handle that as well:

.. code-block:: bash

  nextflow run maxemil/PhyloMagnet --help

  # or
  git clone https://github.com/maxemil/PhyloMagnet

  nextflow run PhyloMagnet/main.nf --help
