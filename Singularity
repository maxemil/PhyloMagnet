Bootstrap: docker
From: debian:stretch

%labels
Maintainer	max-emil.schon@icm.uu.se

%environment
PYTHONPATH='/usr/local/custom_python3_lib/'
export PYTHONPATH
VERSION_MEGAN="6_14_2"
export VERSION_MEGAN
VERSION_PYTHON3=$(python3 --version | cut -d' ' -f2)
export VERSION_PYTHON3
VERSION_TRIMAL=$(trimal --version | cut -d ' ' -f2)
export VERSION_TRIMAL
VERSION_MAFFT="7.407"
export VERSION_MAFFT
VERSION_DIAMOND="0.9.24"
export VERSION_DIAMOND
VERSION_FASTQ_DUMP=$(fastq-dump --version | cut -d':' -f2)
export VERSION_FASTQ_DUMP
VERSION_FASTTREE=$(FastTree 2>&1 >/dev/null | grep Usage | cut -d' ' -f 5)
export VERSION_FASTTREE
VERSION_IQTREE="1.6.9"
export VERSION_IQTREE
VERSION_PRANK=$(prank -version | grep PRANK | cut -d' ' -f4)
export VERSION_PRANK
VERSION_EPA_NG=$(epa-ng --version | cut -d' ' -f 2)
export VERSION_EPA_NG
VERSION_PAPARA=$(papara 2>&1 >/dev/null | grep 'papara_core' | cut -d' ' -f 5)
export VERSION_PAPARA
VERSION_RAXML=$(raxml -version | grep 'RAxML version' | cut -d' ' -f 5)
export VERSION_RAXML
VERSION_GAPPA='v0.0.0'
export VERSION_GAPPA
VERSION_ETE3=$(python3 -c 'import ete3; print(ete3.__version__)')
export VERSION_ETE3
VERSION_BIOPYTHON=$(python3 -c 'import Bio; print(Bio.__version__)')
export VERSION_BIOPYTHON
VERSION_PANDAS=$(python3 -c 'import pandas; print(pandas.__version__)')
export VERSION_PANDAS
VERSION_NUMPY=$(python3 -c 'import numpy; print(numpy.__version__)')
export VERSION_NUMPY

%post
######## base system ########
apt-get update
apt-get clean
apt-get install --no-install-recommends -qy \
                  default-jdk \
                  git \
                  wget \
                  vim \
                  tk \
                  zlib1g-dev \
                  libxml-simple-perl \
                  libtime-piece-perl \
                  libdigest-md5-file-perl \
                  python3 \
                  python3-pip \
                  python3-setuptools \
                  python3-dev \
                  python3-tk \
                  qt5-default \
                  xvfb \
                  autotools-dev \
                  libtool \
                  flex \
                  bison \
                  cmake \
                  automake \
                  autoconf \
                  build-essential \
                  ca-certificates

######## MEGAN6 ########
cd /usr/local/
wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/MEGAN_Community_unix_"$VERSION_MEGAN".sh
chmod +x MEGAN_Community_unix_"$VERSION_MEGAN".sh
./MEGAN_Community_unix_"$VERSION_MEGAN".sh -q

######## python ########
pip3 install wheel
pip3 install biopython ete3 scipy pandas seaborn xvfbwrapper pyqt5 requests
rm /usr/local/bin/ete3

######## PRANK #########
cd /usr/local/
wget http://wasabiapp.org/download/prank/development/prank.source.170703.tgz
tar -xvf prank.source.170703.tgz
cd /usr/local/development/src
make
ln -s /usr/local/development/src/prank /usr/local/bin/prank

######## trimal #########
cd /usr/local
git clone https://github.com/scapella/trimal.git
cd trimal/source
make
ln -s /usr/local/trimal/source/trimal /usr/local/bin/

######## MAFFT #########
cd /usr/local/
wget https://mafft.cbrc.jp/alignment/software/mafft-"$VERSION_MAFFT"-without-extensions-src.tgz
tar xfvz mafft-"$VERSION_MAFFT"-without-extensions-src.tgz
cd  /usr/local/mafft-"$VERSION_MAFFT"-without-extensions/core
make clean
make
make install

######## diamond #########
cd /usr/local
mkdir diamond
cd diamond
wget https://github.com/bbuchfink/diamond/releases/download/v"$VERSION_DIAMOND"/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
ln -s /usr/local/diamond/diamond /usr/local/bin/

######## SRA toolkit #########
apt-get install --no-install-recommends -qy sra-toolkit

######## fasttree #########
cd /usr/local/bin
wget http://www.microbesonline.org/fasttree/FastTreeMP
chmod +x FastTreeMP
cp FastTreeMP FastTree

######## IQtree #########
cd /usr/local
wget https://github.com/Cibiv/IQ-TREE/releases/download/v"$VERSION_IQTREE"/iqtree-"$VERSION_IQTREE"-Linux.tar.gz
tar -xvf iqtree-"$VERSION_IQTREE"-Linux.tar.gz
ln -s /usr/local/iqtree-"$VERSION_IQTREE"-Linux/bin/iqtree /usr/local/bin/iqtree

####### EPA-ng #########
cd /usr/local
git clone https://github.com/Pbdas/epa-ng.git
cd /usr/local/epa-ng
make
ln -s /usr/local/epa-ng/bin/epa-ng /usr/local/bin/epa-ng

######## standard-raxml ########
cd /usr/local/
git clone https://github.com/stamatak/standard-RAxML.git standard-raxml
cd /usr/local/standard-raxml
make -f Makefile.PTHREADS.gcc
rm *.o
ln -s /usr/local/standard-raxml/raxmlHPC-PTHREADS /usr/local/bin/raxml

####### parapa #########
cd /usr/local/bin
wget https://sco.h-its.org/exelixis/resource/download/software/papara_nt-2.5-static_x86_64.tar.gz
tar -xvf papara_nt-2.5-static_x86_64.tar.gz
rm papara_nt-2.5-static_x86_64.tar.gz
mv papara_static_x86_64 papara

####### gappa #########
cd /usr/local/
git clone --recursive https://github.com/lczech/gappa.git
cd /usr/local/gappa
make
ln -s /usr/local/gappa/bin/gappa /usr/local/bin/gappa

%test
/usr/local/megan/tools/daa-meganizer -h
/usr/local/megan/tools/gc-assembler -h
trimal --version
mafft --version
diamond version
fastq-dump --version
which FastTree
iqtree -h
python3 --version
epa-ng --version
raxml -h
papara -h
gappa --help
prank -version
python3 -c "import ete3"

%runscript
echo "MEGAN6: "$VERSION_MEGAN
echo "trimal: "$VERSION_TRIMAL
echo "MAFFT: "$VERSION_MAFFT
echo "diamond: "$VERSION_DIAMOND
echo "fastq-dump: "$VERSION_FASTQ_DUMP
echo "FastTree: "$VERSION_FASTTREE
echo "IQ-TREE: "$VERSION_IQTREE
echo "Python3: "$VERSION_PYTHON3
echo "Prank: "$VERSION_PRANK
echo "RAxML: "$VERSION_RAXML
echo "epa-ng: "$VERSION_EPA_NG
echo "PaPaRa: "$VERSION_PAPARA
echo "gappa: "$VERSION_GAPPA
echo "Biopython: "$VERSION_BIOPYTHON
echo "Pandas: "$VERSION_PANDAS
echo "Numpy: "$VERSION_NUMPY
echo "ETE3: "$VERSION_ETE3
