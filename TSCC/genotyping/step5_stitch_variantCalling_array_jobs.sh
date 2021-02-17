#!/bin/bash

pipeline_arguments=${ARG}
previous_flow_cells_bams=${PREV_BAMS}

home=$(head -n 1 ${pipeline_arguments})
dir_path=$(head -n 2 ${pipeline_arguments} | tail -n 1)
reference_panels=$(head -n 6 ${pipeline_arguments} | tail -n 1)
code=$(head -n 8 ${pipeline_arguments} | tail -n 1)
code=${code}/genotyping/util
bams_data=${dir_path}/bams
stitch_path=${dir_path}/stitch
more_bam=$(awk "{print;}" ${previous_flow_cells_bams})
#### read in arguments for the pipeline

cd ${home}

################ Variant Calling STITCH #######################
START=$(date +%s)

source activate stitch  
Rscript ${code}/variant_calling_STITCH.r \
    ${PBS_ARRAYID} \
    ${stitch_path} \
    ${reference_panels} \
    ${bams_data} \
    ${more_bam}

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "STITCH Time elapsed: $(( $END - $START )) seconds"