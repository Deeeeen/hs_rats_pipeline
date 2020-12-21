![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# Docker
## Source code for HS rats pipeline on [Docker](https://hub.docker.com/u/deeeeen)  

## Links
Github: https://github.com/Deeeeen/hs_rats_pipeline/tree/master/Docker  
Docker Hub: https://hub.docker.com/u/deeeeen   

## Contents
**[hs_rats_genotyping_step1](hs_rats_genotyping_step1)**  
Preparation for the pipeline, which constructs the basic structure of the directory and splits the sample sheet for Fgbio demultiplex based on each Riptide library preparation.  

## Requirements for [pipeline_arguments](pipeline_arguments) file
Line 1: Flow cell directory  
Line 2: Flow cell metadata
Line 3: Sequencing data directory  
Line 4: Reference genome  
Line 5: Reference panels for STITCH  
Line 6: Genetic map for BEAGLE  
Line 7: Directory where you keep the code for the pipeline  
Line 8: The general name of this run  