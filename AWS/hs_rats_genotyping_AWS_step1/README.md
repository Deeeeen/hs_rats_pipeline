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
2. Move the ```metadata``` that Fgbio need and the ```pipeline_argument``` file to the S3 bucket.
3. Please update the ```AWS_config/batch_computeEnv.json``` file with the resources you want to use for this job. IMPORTANT: ```instanceRole``` should be the AWS IAM Role that has the premission for ```AmazonEC2ContainerServiceforEC2Role```, and ```serviceRole``` should have premission for ```AWSBatchServiceRole```. If you want to use Spot instances or anything else, please edit accordingly. You can find more information about AWS Batch Compute Environment Parameters [here](https://docs.aws.amazon.com/batch/latest/userguide/compute_environment_parameters.html).
4. Please update the ```AWS_config/batch_jobQueue.json``` file with the resources you want to use for this job. IMPORTANT: ```computeEnvironment``` should be the compute environment you created on step 3. You can find more information about AWS Batch Job Queue Parameters [here](https://docs.aws.amazon.com/batch/latest/userguide/job_queue_parameters.html).
5. Please update the ```AWS_config/batch_jobDef.json``` file with the resources you want to use for this job. IMPORTANT: ```jobRoleArn``` should be the AWS IAM Role that has the premission for ```AmazonS3FullAccess```, and ```environment``` should be the URL for the S3 bucket you created on step 1. You can find more information about AWS Batch Job Definition Parameters [here](https://docs.aws.amazon.com/batch/latest/userguide/job_definition_parameters.html).

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
