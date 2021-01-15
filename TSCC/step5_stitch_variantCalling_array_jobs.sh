#!/bin/bash
#PBS -q hotel
#PBS -N stitch
#PBS -l nodes=1:ppn=24
#PBS -l walltime=100:00:00
#PBS -t 1-22%4
#PBS -V
#PBS -j oe
#PBS -k oe
#PBS -M dec037@health.ucsd.edu
#PBS -m ae

#### Make sure to change the compute node, job name, num of node,
#### ppn, walltime, forward email address for notifications. (above)

pipeline_arguments=$ARG
previous_flow_cells_bams=$PREV_BAMS

home=$(head -n 1 ${pipeline_arguments})
dir_path=$(head -n 2 ${pipeline_arguments} | tail -n 1)
reference_panels=$(head -n 6 ${pipeline_arguments} | tail -n 1)
code=$(head -n 8 ${pipeline_arguments} | tail -n 1)
bams_data=${dir_path}/bams
more_bam=$(awk "{print;}" ${previous_flow_cells_bams})
#### read in arguments for the pipeline

cd ${home}

################ Variant Calling STITCH #######################
START=$(date +%s)

source activate stitch  
Rscript ${code}/variant_calling_STITCH.r \
    $PBS_ARRAYID \
    ${dir_path}/stitch \
    ${reference_panels} \
    ${bams_data} \
    ${more_bam}

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "STITCH Time elapsed: $(( $END - $START )) seconds"