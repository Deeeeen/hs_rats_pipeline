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

# Make directories to hold code and data
RUN mkdir /hs_genotyping_step1

RUN mkdir /hs_genotyping_step1/code /hs_genotyping_step1/data

# Use COPY command to copy code into the container
# (there are other ways to do this)
COPY ./code /hs_genotyping_step1/code

# Change permission in order to run the bash script
RUN chmod 777 /hs_genotyping_step1/code/step1_prep.sh

# Set working directory
WORKDIR /hs_genotyping_step1

# Run bash script
CMD ["/hs_genotyping_step1/code/step1_prep.sh"]