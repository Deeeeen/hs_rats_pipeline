![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# [Palmer Lab](https://palmerlab.org/), HS Rats Genotyping Pipeline Step1  
## Overview  
Preparation for the pipeline, which constructs the basic structure of the directory and splits the sample sheet for Fgbio demultiplex based on each Riptide library preparation.  

## Links
Github: https://github.com/Deeeeen/hs_rats_pipeline/tree/master/Docker/hs_rats_genotyping_step1  
Docker Hub: https://hub.docker.com/r/deeeeen/hs_rats_pipeline_step1  

## Documentation
### Run Docker image on your local PC:
1. Install [Docker](https://docs.docker.com/get-docker/).
2. Make sure Docker is running with running the following code on terminal.
```
docker version
```
3. Download the docker image from Docker Hub.
```
docker pull deeeeen/hs_rats_pipeline_step1:latest
```
4. Run following command to get the image ID.
```
docker images
```
5. Run docker image. Here we use shared volumes ```-v``` flag to transfer data from host machine to Docker. Please specify a directory on your host machine here: ```<host directory>``` and specify the docker image ID you got on step 4 here: ```<IMAGE ID>```.   
```
docker run -v "<host directory>:/hs_genotyping_step1/data" <IMAGE ID>
```

### Run Docker image on TSCC (PBS):
1. Download the docker image from Docker Hub.
```
singularity pull docker://deeeeen/hs_rats_pipeline_step1:latest
```
2. Run docker image. Here we use shared volumes ```-B``` flag to transfer data from host machine to Docker. Please specify a directory on your host machine here: ```<host directory>``` and specify the path to the docker image you got on step 1 here: ```<hs_rats_pipeline_step1_latest.sif>```.   
```
singularity run -B <host directory>:/hs_genotyping_step1/data <hs_rats_pipeline_step1_latest.sif> 
```

## Requirements
1. ```<host directory>``` must have ```pipeline_arguments``` and sample sheet for Fgbio.
2. ```pipeline_arguments``` must have first 2 required lines. Please refer to [here](https://github.com/Deeeeen/hs_rats_pipeline/blob/master/Docker/README.md) for required details on ```pipeline_arguments```.