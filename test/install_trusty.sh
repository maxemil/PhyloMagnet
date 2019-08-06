sudo apt-get -qq update
sudo apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    golang-go


echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
  echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
  source ~/.bashrc

export VERSION=3.3.0 && # adjust this as necessary \
  wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz && \
  tar -xzf singularity-${VERSION}.tar.gz && \
  cd singularity

./mconfig && \
  make -C builddir && \
  sudo make -C builddir install
#
# VERSION=2.4
# wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz
# tar xvf singularity-$VERSION.tar.gz
# cd singularity-$VERSION
# ./configure --prefix=/usr/local
# make
# sudo make install

wget -qO- get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
