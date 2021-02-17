![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# Software 
All the software installed on TSCC are [here](https://aapalmer-lab.slack.com/files/T0JULRU14/FPS2923NU).  
This github page will instruct you to install all the require software for HS rats genotyping pipeline.  

## Contents
**[hs_rats_conda_env.yml](hs_rats_conda_env.yml)**  
yml file with dependencies and reuqired packages to create a conda environment for the pipeline.  

**[stitch_conda_env.yml](stitch_conda_env.yml)**  
yml file with dependencies and reuqired packages to create a conda environment for STITCH.  

## Documentation  
### Before running the pipeline on TSCC:
1. Install [miniconda3](https://docs.conda.io/en/latest/miniconda.html)
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
2. Create conda environment with yml file for the genotyping pipeline
```
conda env create -n hs_rats --file software/hs_rats_conda_env.yml
```  
3. Create conda environment with yml file for STITCH
```
conda env create -n stitch --file software/stitch_conda_env.yml
```  