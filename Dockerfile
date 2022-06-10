FROM ubuntu:20.04
MAINTAINER Lindo Nkambule (lindonkambule116@gmail.com)

ARG SAMTOOLS_VERSION=1.13

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        autoconf \
        automake \
        bzip2 \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        g++ \
        gcc \
        git \
        gzip \
        libboost-all-dev \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        make \
        python3 \
        python3-pip \
        r-mathlib \
        sudo \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

# HTSLIB
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xf htslib-${SAMTOOLS_VERSION}.tar.bz2 && rm htslib-${SAMTOOLS_VERSION}.tar.bz2 && cd htslib-${SAMTOOLS_VERSION} && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && make install && make clean

COPY makefile_shapeit4 /opt

# SHAPEIT4
RUN git clone https://github.com/odelaneau/shapeit4.git && \
    cd shapeit4 && \
    mv makefile makefile.old && cp /opt/makefile_shapeit4 . && mv makefile_shapeit4 makefile && \
    make && \
    cd /shapeit4/maps && mkdir b37 b38 && gunzip *.gz && \
    tar -xf genetic_maps.b37.tar -C b37/ && \
    tar -xf genetic_maps.b38.tar -C b38/ && \
    rm *.tar

ENV PATH /shapeit4/bin/:${PATH}

# BCFTOOLS
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/${SAMTOOLS_VERSION}/bcftools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xf bcftools-${SAMTOOLS_VERSION}.tar.bz2 && rm bcftools-${SAMTOOLS_VERSION}.tar.bz2 && cd bcftools-${SAMTOOLS_VERSION} && \
    ./configure --with-htslib=/opt/htslib-${SAMTOOLS_VERSION} && make && make install && make clean

# EAGLE
RUN cd /opt && \
    wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz && \
    gunzip Eagle_v2.4.1.tar.gz && \
    tar xvf Eagle_v2.4.1.tar && \
    cp /opt/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz /opt && \
    cp /opt/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz /opt && \
    mv Eagle_v2.4.1/eagle /usr/local/bin/ && \
    rm -rf Eagle_v2.4.1*

# IMPUTE5
COPY  impute5_v1.1.5.zip /opt
RUN cd /opt && \
    unzip impute5_v1.1.5.zip && cd impute5_v1.1.5 && \
    mv *_static /usr/local/bin/ && cd /opt && rm -rf impute5_v1.1.5*

# makeScaffold for building haplotype scaffolds for phasing
RUN pip3 install Cython --install-option="--no-cython-compile"

RUN git clone https://github.com/sinanshi/makeScaffold.git && \
    cd makeScaffold && rm makefile && \
    cmake . && \
    make && \
    mv src/scaffold /usr/local/bin/


