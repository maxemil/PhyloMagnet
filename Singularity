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
/opt/conda/envs/PhyloMagnet-0.6/bin/trimal --version
/opt/conda/envs/PhyloMagnet-0.6/bin/mafft --version
/opt/conda/envs/PhyloMagnet-0.6/bin/diamond version
/opt/conda/envs/PhyloMagnet-0.6/bin/fastq-dump --version
which /opt/conda/envs/PhyloMagnet-0.6/bin/FastTree
/opt/conda/envs/PhyloMagnet-0.6/bin/iqtree -h
/opt/conda/envs/PhyloMagnet-0.6/bin/python3 --version
/opt/conda/envs/PhyloMagnet-0.6/bin/raxmlHPC-PTHREADS -version
/opt/conda/envs/PhyloMagnet-0.6/bin/epa-ng --version
/opt/conda/envs/PhyloMagnet-0.6/bin/gappa --help
/opt/conda/envs/PhyloMagnet-0.6/bin/prank -version
/opt/conda/envs/PhyloMagnet-0.6/bin/python3 -c "import ete3"
/usr/local/bin/papara
