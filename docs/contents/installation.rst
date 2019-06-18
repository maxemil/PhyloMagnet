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

Then, download the container from singularity-hub or build it locally with the singularity recipe:

.. code-block:: bash
  # singularity 2

  singularity pull --name PhyloMagnet.simg shub://maxemil/PhyloMagnet:latest

  sudo singularity build PhyloMagnet.simg Singularity

  # singularity 3

  singularity pull library://maxemil/default/phylomagnet:0.6

  sudo singularity build PhyloMagnet.sif Singularity

.. note::

  We highly recommend using the provided Singularity container to install all needed Software! Installing everything manually can be very cumbersome. If you want to do that, look at the %post section of the singularity file to get an idea what tools are needed and how to install them

PhyloMagnet
--------

now you can either get PhyloMagnet from github (clone or download from github.com/maxemil/PhyloMagnet) or let Nextflow handle that as well:

.. code-block:: bash

  nextflow run maxemil/PhyloMagnet --help

  # or
  git clone github.com/maxemil/PhyloMagnet

  nextflow run PhyloMagnet/main.nf --help
