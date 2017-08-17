Bootstrap: docker
From: finalduty/archlinux:daily

%post

######## base system ########
echo "Server = http://mirror.de.leaseweb.net/archlinux/\$repo/os/\$arch" >> /etc/pacman.d/mirrorlist
echo "[lambdait]" >> /etc/pacman.conf
echo "SigLevel = Never" >> /etc/pacman.conf
echo "Server = https://lambda.informatik.uni-tuebingen.de/repo/mypkgs/" >> /etc/pacman.conf

pacman -Syu --noconfirm
pacman -S --noconfirm base-devel jdk git wget expect tk python-pyqt4

######## MEGAN6 ########
cd /usr/local/

wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/MEGAN_Community_unix_6_8_12.sh

chmod +x MEGAN_Community_unix_6_8_12.sh
./MEGAN_Community_unix_6_8_12.sh -q

######## python ########
pacman -S --noconfirm python3 python-pip

pip install biopython ete3 scipy pandas seaborn xvfbwrapper

mkdir /usr/local/custom_python3_lib/

######## MAFFT #########
pacman -S --noconfirm mafft

######## trimal #########
git clone https://github.com/scapella/trimal.git

cd trimal/source
make
ln -s /usr/local/trimal/source/trimal /usr/local/bin/

cd ../..

######## diamond #########
mkdir diamond
cd diamond
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.9/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
ln -s /usr/local/diamond/diamond /usr/local/bin/

cd ..

######## SRA toolkit #########
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -xvf sratoolkit.current-centos_linux64.tar.gz
ln -s /usr/local/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump /usr/local/bin/

######## fasttree #########
cd /usr/local/bin
wget http://www.microbesonline.org/fasttree/FastTreeMP
chmod +x FastTreeMP
ln -s /usr/local/bin/FastTreeMP /usr/local/bin/FastTree

######## IQtree #########
cd /usr/local
wget https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.beta4/iqtree-1.6.beta4-Linux.tar.gz
tar -xvf iqtree-1.6.beta4-Linux.tar.gz
ln -s /usr/local/iqtree-1.6.beta4-Linux/bin/iqtree /usr/local/bin/iqtree-1.6.beta4


%files
/local/two/Software/python_lib/*.py /usr/local/custom_python3_lib/

%environment
PYTHONPATH='/usr/local/custom_python3_lib/'
export PYTHONPATH

%test
/usr/local/megan/tools/daa-meganizer -h
/usr/local/megan/tools/gc-assembler -h
trimal --version
mafft --version
diamond version
fastq-dump --version
FastTree
iqtree-1.6.beta4 -h

%labels
Maintainer	max-emil.schon@icm.uu.se
