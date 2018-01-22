sudo apt-get -qq update
sudo apt-get install singularity-container

wget -qO- get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

singularity pull --name PhyloMagnet.simg shub://maxemil/PhyloMagnet:latest
