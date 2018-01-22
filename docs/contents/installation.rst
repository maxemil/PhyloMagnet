Installation Instructions for PhyloMagnet
=========================================
Nextflow
--------
To run PhyloMagnet, you will need to install `Nextflow <https://www.nextflow.io/>`_, a pipeline execution framwork. This is however quite simple:

.. code-block:: bash

  curl -s https://get.nextflow.io | bash


Singularity
-----------
We provide a singularity container with all necessary tools installed and configured. To use it, you first need to install `Singularity <http://singularity.lbl.gov/install-linux>`_ itself:

.. code-block:: bash

  VERSION=2.4
  wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz
  tar xvf singularity-$VERSION.tar.gz
  cd singularity-$VERSION
  ./configure --prefix=/usr/local
  make
  sudo make install


Then, download the container from singularity-hub or build it locally with the singularity recipe:

.. code-block:: bash

  singularity pull --name PhyloMagnet.simg shub://maxemil/PhyloMagnet:latest

  sudo singularity build PhyloMagnet.simg Singularity


When building the image from the singularity recipe yourself, you can adjust the JVM options for MEGAN by changing the file MEGAN.vmoptions. This is necessary if you are planning to use large FastQ files or need other specific options (especially memory requirements) for the daa-meganizer and the gc-assembler programs.

.. note::

  We highly recommend using the provided Singularity container to install all needed Software! Installing everything manually can be very cumbersome. If you want to do that, look at the %post section of the singularity file for the needed tools and how to install them
