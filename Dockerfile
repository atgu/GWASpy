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

COPY makefile /opt

# SHAPEIT4
RUN git clone https://github.com/odelaneau/shapeit4.git && \
    cd shapeit4 && \
    mv makefile makefile.old && cp /opt/makefile . && \
    make

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
    mv Eagle_v2.4.1/eagle /usr/local/bin/

COPY makefile split_maps.sh /opt

# genetic maps files
RUN cd /opt && \
    mkdir genetic_maps_eagle && cd genetic_maps_eagle && \
    mkdir hg38 hg17 hg19 && cd /opt && \
    mkdir genetic_maps_shapeit && cd genetic_maps_shapeit && \
    mkdir hg38 hg17 hg19 && cd /opt &&  \
    bash split_maps.sh && \
    rm -rf Eagle_v2.4.1*



