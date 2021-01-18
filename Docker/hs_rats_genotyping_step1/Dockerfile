#################################################################
# Dockerfile
#
# Version:          dev
# Software:         HS Rats Genotyping Pipeline Step 1
# Software Version: dev
# about.summary:    HS Rats Genotyping Pipeline Step 1
# about.tags:       Genomics
# BASE IMAGE :      alpine:latest
#################################################################
# Base image
FROM alpine:latest

# Install bash in order to run bash script
RUN apk add --no-cache --upgrade bash

ENV GLIBC_VER=2.32-r0

# Install glibc compatibility and miniconda3 for alpine
RUN apk --no-cache add \
        binutils \
        curl \
    && curl -sL https://alpine-pkgs.sgerrand.com/sgerrand.rsa.pub -o /etc/apk/keys/sgerrand.rsa.pub \
    && curl -sLO https://github.com/sgerrand/alpine-pkg-glibc/releases/download/${GLIBC_VER}/glibc-${GLIBC_VER}.apk \
    && curl -sLO https://github.com/sgerrand/alpine-pkg-glibc/releases/download/${GLIBC_VER}/glibc-bin-${GLIBC_VER}.apk \
    && apk add --no-cache \
        glibc-${GLIBC_VER}.apk \
        glibc-bin-${GLIBC_VER}.apk \
    && curl -sL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda3 \
    && rm Miniconda3-latest-Linux-x86_64.sh \
    && apk --no-cache del \
        binutils \
        curl \
    && rm glibc-${GLIBC_VER}.apk \
    && rm glibc-bin-${GLIBC_VER}.apk \
    && rm -rf /var/cache/apk/*

# Make directories to hold code and data
RUN mkdir /pipeline \
    && mkdir /pipeline/code /pipeline/data

# Use COPY command to copy code into the container
# (there are other ways to do this)
COPY ./code /pipeline/code

# Create conda environment for the pipeline
ENV PATH=/miniconda3/bin:$PATH
RUN conda env create --file=/pipeline/code/environment.yml \
    && conda update --all && conda clean -p

# Change permission in order to run the bash script
RUN chmod 770 /pipeline/code/step1_prep.sh

# Set working directory
WORKDIR /pipeline

# Run bash script
CMD ["/pipeline/code/step1_prep.sh"]