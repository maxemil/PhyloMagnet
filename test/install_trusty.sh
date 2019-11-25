SINGULARITY_BASE="${GOPATH}/src/github.com/sylabs/singularity"
export PATH="${GOPATH}/bin:${PATH}"

mkdir -p "${GOPATH}/src/github.com/sylabs"
cd "${GOPATH}/src/github.com/sylabs"

git clone -b vault/release-3.4 https://github.com/sylabs/singularity
cd singularity

./mconfig -v -p /usr/local
make -j `nproc 2>/dev/null || echo 1` -C ./builddir all
sudo make -C ./builddir install
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
