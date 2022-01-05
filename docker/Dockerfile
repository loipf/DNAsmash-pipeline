FROM continuumio/miniconda3@sha256:7838d0ce65783b0d944c19d193e2e6232196bada9e5f3762dc7a9f07dc271179

MAINTAINER loipf

RUN conda update --all

RUN apt-get update && apt install -y procps gzip curl chrpath gcc
#RUN apt install -y build-essential aptitude libxml2-dev libssl-dev libcurl4-openssl-dev libxtst6 libgl1-mesa-glx

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

### bioinformatic tools
RUN conda install r-base==4.0.5 r-essentials
RUN R -e 'install.packages("pacman",repos = "http://cran.us.r-project.org")'
RUN R -e 'pacman::p_load(data.table, pheatmap, BiocManager)'
RUN R -e 'BiocManager::install("PoisonAlien/maftools")'

RUN conda install samtools==1.11.0

### download URMAP
RUN wget https://github.com/rcedgar/urmap/raw/master/bin/urmap
RUN chmod +x ./urmap



