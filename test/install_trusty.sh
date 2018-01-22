sudo apt-get -qq update

VERSION=2.4
wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz
tar xvf singularity-$VERSION.tar.gz
cd singularity-$VERSION
./configure --prefix=/usr/local
make
sudo make install

wget -qO- get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

singularity pull --name PhyloMagnet.simg shub://maxemil/PhyloMagnet:latest
