from ubuntu:20.04

RUN apt-get update -y --fix-missing
RUN DEBIAN_FRONTEND="noninteractive" apt-get install -y software-properties-common git python python3 python3-pip emacs curl wget zip unzip tar zlib1g-dev liblzo2-dev libcurl4-openssl-dev libssl-dev libncurses5 tzdata

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
RUN conda install -y -c bioconda perl-parallel-forkmanager=2.02
RUN conda clean -ay

# These commands should install R-4.2.0
# Maybe try to specify version to make installation more stable

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-get install -y r-base


# Install R packages

RUN R --slave -e 'install.packages("remotes")'
RUN R --slave -e 'remotes::install_version("optparse", "1.3.2")'
RUN R --slave -e 'remotes::install_version("tidyverse", "1.2.1")'
RUN R --slave -e 'remotes::install_version("XML", "3.99-0.9")'
RUN R --slave -e 'remotes::install_version("randomForest", "4.6-14")'

COPY ./code/ /code/

# After release 1.3.5 these lines can be used to automate docker image creation process further
#RUN wget https://github.com/TRON-Bioinformatics/EasyFuse/archive/refs/tags/v1.3.5.tar.gz && \
#tar xvfz v1.3.5.tar.gz && \
#mv EasyFuse-1.3.5/ /code/easyfuse/ && \
#cd /code/easyfuse/ && \
#mv config.ini.docker config.ini && \
#mv build_env.sh.sample build_env.sh && \
#mv blacklist.txt.sample blacklist.txt

ENV PERL5LIB "$PERL5LIB:/code/SOAPfuse/1.27/source/bin/perl_module/"
