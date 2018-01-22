wget -O- http://neuro.debian.net/lists/trusty.de-md.libre | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
sudo apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9

sudo apt-get -qq update
sudo apt-get install singularity-container

wget -qO- get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

singularity pull --name PhyloMagnet.simg shub://maxemil/PhyloMagnet:latest
