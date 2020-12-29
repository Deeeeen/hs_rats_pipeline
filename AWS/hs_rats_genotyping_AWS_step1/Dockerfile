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

# Install glibc compatibility ans AWS-CLI V2 for alpine
ENV GLIBC_VER=2.31-r0

RUN apk --no-cache add \
        binutils \
        curl \
    && curl -sL https://alpine-pkgs.sgerrand.com/sgerrand.rsa.pub -o /etc/apk/keys/sgerrand.rsa.pub \
    && curl -sLO https://github.com/sgerrand/alpine-pkg-glibc/releases/download/${GLIBC_VER}/glibc-${GLIBC_VER}.apk \
    && curl -sLO https://github.com/sgerrand/alpine-pkg-glibc/releases/download/${GLIBC_VER}/glibc-bin-${GLIBC_VER}.apk \
    && apk add --no-cache \
        glibc-${GLIBC_VER}.apk \
        glibc-bin-${GLIBC_VER}.apk \
    && curl -sL https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip -o awscliv2.zip \
    && unzip awscliv2.zip \
    && aws/install \
    && rm -rf \
        awscliv2.zip \
        aws \
        /usr/local/aws-cli/v2/*/dist/aws_completer \
        /usr/local/aws-cli/v2/*/dist/awscli/data/ac.index \
        /usr/local/aws-cli/v2/*/dist/awscli/examples \
    && apk --no-cache del \
        binutils \
        curl \
    && rm glibc-${GLIBC_VER}.apk \
    && rm glibc-bin-${GLIBC_VER}.apk \
    && rm -rf /var/cache/apk/*

RUN aws --version

# Make directories to hold code and data
RUN mkdir /hs_genotyping_step1

RUN mkdir /hs_genotyping_step1/code /hs_genotyping_step1/data

# Use COPY command to copy code into the container
# (there are other ways to do this)
COPY ./code /hs_genotyping_step1/code

# Change permission in order to run the bash script
RUN chmod u+x /hs_genotyping_step1/code/run_code.sh
RUN chmod u+x /hs_genotyping_step1/code/step1_prep.sh

# Set working directory
WORKDIR /hs_genotyping_step1

# Run bash script
CMD ["/hs_genotyping_step1/code/run_code.sh"]