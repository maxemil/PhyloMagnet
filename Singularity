Bootstrap: docker
From: finalduty/archlinux:daily

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
echo "Server = http://mirror.de.leaseweb.net/archlinux/\$repo/os/\$arch" >> /etc/pacman.d/mirrorlist
echo "[lambdait]" >> /etc/pacman.conf
echo "SigLevel = Never" >> /etc/pacman.conf
echo "Server = https://lambda.informatik.uni-tuebingen.de/repo/mypkgs/" >> /etc/pacman.conf

pacman -Syu --noconfirm
pacman -S --noconfirm base-devel \
                      jdk \
                      git \
                      wget \
                      expect \
                      tk \
                      python-pyqt4 \
                      vim \
                      perl-time-piece \
                      perl-xml-simple \
                      perl-digest-md5 \
                      cpanminus

######## MEGAN6 ########
cd /usr/local/
wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/MEGAN_Community_unix_6_8_12.sh
chmod +x MEGAN_Community_unix_6_8_12.sh
./MEGAN_Community_unix_6_8_12.sh -q
mv /usr/local/MEGAN.vmoptions /usr/local/megan/MEGAN.vmoptions

######## python ########
pacman -S --noconfirm python3 python-pip
pip install biopython ete3 scipy pandas seaborn xvfbwrapper
mkdir -p /usr/local/custom_python3_lib/
mv /usr/local/*.py /usr/local/custom_python3_lib/

######## MAFFT #########
pacman -S --noconfirm mafft

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
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.9/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
ln -s /usr/local/diamond/diamond /usr/local/bin/

######## SRA toolkit #########
cd /usr/local
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

######## SEQTK #########
cd /usr/local
git clone https://github.com/lh3/seqtk.git
cd seqtk
make
ln -s /usr/local/seqtk/seqtk /usr/local/bin/seqtk

######## MEGAHIT #########
cd /usr/local
git clone https://github.com/voutcn/megahit.git
cd megahit
make
ln -s /usr/local/megahit/megahit /usr/local/bin/megahit

######## BioPerl #########
cd /usr/local/
/usr/bin/vendor_perl/cpanm Bio::Perl

######## PROKKA-partial #########
cd /usr/local
git clone https://github.com/jennahd/prokka.git
for binary in /usr/local/prokka/bin/*; do ln -s $binary /usr/local/bin; done
prokka --setupdb

%test
/usr/local/megan/tools/daa-meganizer -h
/usr/local/megan/tools/gc-assembler -h
trimal --version
mafft --version
diamond version
fastq-dump --version
which FastTree
iqtree-1.6.beta4 -h
which seqtk
megahit --version
prokka --version
python3 --version
