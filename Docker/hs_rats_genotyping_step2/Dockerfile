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

# Install java SE8
RUN { \
		echo '#!/bin/sh'; \
		echo 'set -e'; \
		echo; \
		echo 'dirname "$(dirname "$(readlink -f "$(which javac || which java)")")"'; \
	} > /usr/local/bin/docker-java-home \
	&& chmod +x /usr/local/bin/docker-java-home
ENV JAVA_HOME /usr/lib/jvm/java-1.8-openjdk
ENV PATH $PATH:/usr/lib/jvm/java-1.8-openjdk/jre/bin:/usr/lib/jvm/java-1.8-openjdk/bin

ENV JAVA_VERSION 8u212
ENV JAVA_ALPINE_VERSION 8.272.10-r4

RUN set -x \
	&& apk add --no-cache \
		openjdk8="$JAVA_ALPINE_VERSION" \
	&& [ "$JAVA_HOME" = "$(docker-java-home)" ]

# Make directories to hold code and data
RUN mkdir /pipeline \
    && mkdir /pipeline/code /pipeline/data /pipeline/software \
    && wget https://github.com/fulcrumgenomics/fgbio/releases/download/1.3.0/fgbio-1.3.0.jar -P /pipeline/software 

# Use COPY command to copy code into the container
# (there are other ways to do this)
COPY ./code /pipeline/code

# Change permission in order to run the bash script
RUN chmod 770 /pipeline/code/step2_demux_array_jobs.sh

# Set working directory
WORKDIR /pipeline

# Run bash script
CMD ["/pipeline/code/step2_demux_array_jobs.sh"]
