Bootstrap: docker
From: debian:stretch

%files
lib/*.py /usr/local/
MEGAN.vmoptions /usr/local/

%labels
Maintainer	max-emil.schon@icm.uu.se

%environment
PYTHONPATH='/usr/local/custom_python3_lib/'
export PYTHONPATH

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
                  python-qt4 \
                  libxml-simple-perl \
                  libtime-piece-perl \
                  libdigest-md5-file-perl \
                  cpanminus \
                  python3 \
                  python3-pip \
                  python3-setuptools \
                  python3-dev \
                  xvfb \
                  build-essential

######## MEGAN6 ########
cd /usr/local/
wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/MEGAN_Community_unix_6_11_1.sh
chmod +x MEGAN_Community_unix_6_11_1.sh
./MEGAN_Community_unix_6_11_1.sh -q
mv /usr/local/MEGAN.vmoptions /usr/local/megan/MEGAN.vmoptions

######## python ########
pip3 install biopython ete3 scipy pandas seaborn xvfbwrapper
mkdir -p /usr/local/custom_python3_lib/
mv /usr/local/*.py /usr/local/custom_python3_lib/

######## MAFFT #########
cd /usr/local/
wget https://mafft.cbrc.jp/alignment/software/mafft-7.312-without-extensions-src.tgz
tar xfvz mafft-7.312-without-extensions-src.tgz
cd  /usr/local/mafft-7.312-without-extensions/core
make clean
make
make install

######## trimal #########
cd /usr/local
git clone https://github.com/scapella/trimal.git
cd trimal/source
make
ln -s /usr/local/trimal/source/trimal /usr/local/bin/

######## diamond #########
cd /usr/local
mkdir diamond
cd diamond
wget https://github.com/bbuchfink/diamond/releases/download/v0.9.21/diamond-linux64.tar.gz
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
wget https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.3/iqtree-1.6.3-Linux.tar.gz
tar -xvf iqtree-1.6.3-Linux.tar.gz
ln -s /usr/local/iqtree-1.6.3-Linux/bin/iqtree /usr/local/bin/iqtree

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
