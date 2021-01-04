![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# [Palmer Lab](https://palmerlab.org/), HS Rats Genotyping Pipeline AWS Step1  
## Overview  
Preparation for the pipeline, which constructs the basic structure of the directory and splits the sample sheet for Fgbio demultiplex based on each Riptide library preparation.  
This is the first step of the HS rats genotyping pipeline. The Docker image utilizes AWS Batch and AWS S3. It fetches data from AWS S3 bucket, processes the data on EC2 instance by using AWS Batch service, and outputs the results to the AWS S3 bucket.

## Links
Github: https://github.com/Deeeeen/hs_rats_pipeline/tree/master/AWS/hs_rats_genotyping_AWS_step1  
Docker Hub: https://hub.docker.com/r/deeeeen/hs_rats_genotyping_aws_step1  

## Documentation
### Before running the Docker image on AWS Batch:
1. Follow [Setting Up Amazon S3](https://docs.aws.amazon.com/AmazonS3/latest/gsg/SigningUpforS3.html) to create an AWS S3 bucket.

### Run Docker image on AWS Batch:
1. Follow [AWS Batch Setting Up tutorial](https://docs.aws.amazon.com/batch/latest/userguide/get-set-up-for-aws-batch.html) to set up your AWS account.
2. Install [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html) on your local PC.
3. Build a Compute Environment for your AWS Batch job with ```AWS_config/batch_computeEnv.json```.
```
aws batch create-compute-environment --cli-input-json file://computeEnv.json --region <region>
```
4. Build a Job Queue for your AWS Batch job with ```AWS_config/batch_jobQueue.json```.
```
aws batch create-job-queue --cli-input-json file://batch_jobQueue.json --region <region>
```
5. Register a Job Definition for your AWS Batch job with ```AWS_config/batch_jobDef.json```.
```
aws batch register-job-definition --cli-input-json file://batch_jobDef.json --region <region>
```
6. Submit a AWS Batch Job.
```
aws batch submit-job --job-name <job name> --job-queue <job queue>  --job-definition <job definition name> --region <regio>
```
## Requirements
1. ```<host directory>``` must have ```pipeline_arguments``` and sample sheet for Fgbio.
2. ```pipeline_arguments``` must have first 2 required lines. Please refer to [here](https://github.com/Deeeeen/hs_rats_pipeline/blob/master/Docker/README.md) for required details on ```pipeline_arguments```.