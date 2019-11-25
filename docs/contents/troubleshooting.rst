Troubleshooting
=============

Out of Memory
-------
A common source of error is the memory allocation to MEGAN6. If the MEGAN process crashes with a 'Out of Memory' Error, try to increase the parameter in the 'MEGAN.vmoptions' file, e.g. set it to 16GB:

.. code-block:: bash

    -Xmx16G

The memory should be around the same as the size of the query (fastq) file.

data.jar not found
--------
Another problem might be that the 'includeLocalRef' process exits saying it cannot find the 'data.jar'. Make sure you are using the location of MEGAN6, which is where PhyloMagnet copies the jar from.
