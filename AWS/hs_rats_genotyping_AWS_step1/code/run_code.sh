#!/bin/bash
aws s3 cp ${BATCH_FILE_S3_URL} /hs_genotyping_step1/data --recursive || error_exit "Failed to download S3 folder."

/hs_genotyping_step1/code/step1_prep.sh

pipeline_arguments=/hs_genotyping_step1/data/pipeline_arguments

aws s3 sync /hs_genotyping_step1/data ${BATCH_FILE_S3_URL}
