from ubuntu:latest

RUN apt-get update -y
RUN apt-get upgrade -y
RUN DEBIAN_FRONTEND="noninteractive" apt-get install -y python python3 python-pip python3-pip emacs curl wget zip unzip tar zlib1g-dev liblzo2-dev tzdata

RUN curl -O https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
RUN bash Miniconda2-latest-Linux-x86_64.sh -b -p /root/miniconda2
RUN rm Miniconda2-latest-Linux-x86_64.sh
ENV PATH /root/miniconda2/bin:$PATH
RUN conda update conda
RUN conda update --all
RUN conda install -y -c anaconda pandas=0.24.0 biopython=1.73
RUN conda install -y -c conda-forge python-lzo=1.12 python-xxhash=1.4.3
RUN conda install -y -c bioconda pysam=0.15.2 star=2.6.1b samtools=1.9.0 mapsplice=2.2.1
RUN conda install -y -c bioconda star-fusion=1.5.0 bowtie=1.1.2
RUN conda install -y -c bioconda bowtie2=2.3.4.3
RUN conda install -y -c bioconda bx-python=0.8.2 gffutils=0.10.1
RUN conda install -y -c bioconda fastqc=0.11.9
RUN conda install -y -c bioconda r-optparse
RUN conda install -y -c r r=3.6.0 r-xml=3.98 r-tidyverse=1.2.1 r-randomforest=4.6_14
RUN conda install -y -c bioconda perl-parallel-forkmanager
RUN conda clean -ay

COPY ./code/ /code/
ENV PERL5LIB "$PERL5LIB:/code/SOAPfuse/1.27/source/bin/perl_module/"

#ENV PATH "$PATH:/code/STAR-2.6.1d/bin/Linux_x86_64"
#ENV PATH "$PATH:/code/ncbi-blast-2.8.1+/bin"
