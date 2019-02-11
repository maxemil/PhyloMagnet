From:continuumio/miniconda3:4.5.12
Bootstrap:docker

%labels
    MAINTAINER Max Emil Schön <max-emil.schon@icm.uu.se>
    DESCRIPTION Singularity image containing all requirements for the PhyloMagnet pipeline
    VERSION 0.6

%environment
    PATH=/opt/conda/envs/PhyloMagnet-0.6/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    apt-get update
    apt-get install -y procps libxtst6
    apt-get clean -y

    ####### parapa #########
    cd /usr/local/bin
    wget https://sco.h-its.org/exelixis/resource/download/software/papara_nt-2.5-static_x86_64.tar.gz
    tar -xvf papara_nt-2.5-static_x86_64.tar.gz
    rm papara_nt-2.5-static_x86_64.tar.gz
    mv papara_static_x86_64 papara

    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a


%test
/opt/conda/envs/PhyloMagnet-0.6/opt/megan-6.12.3/tools/daa-meganizer -h
/opt/conda/envs/PhyloMagnet-0.6/opt/megan-6.12.3/tools/gc-assembler -h
trimal --version
mafft --version
diamond version
fastq-dump --version
which FastTree
iqtree -h
python3 --version
raxmlHPC-PTHREADS -version
epa-ng --version
gappa --help
prank -version
python3 -c "import ete3"
