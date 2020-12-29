#!/bin/bash
aws s3 cp ${BATCH_FILE_S3_URL} /hs_genotyping_step1/data --recursive || error_exit "Failed to download S3 folder."

/hs_genotyping_step1/code/step1_prep.sh

pipeline_arguments=/hs_genotyping_step1/data/pipeline_arguments
dir_path=/hs_genotyping_step1/data/$(head -n 1 ${pipeline_arguments} | tail -n 1)

aws s3 cp ${dir_path} ${BATCH_FILE_S3_URL} --recursive
