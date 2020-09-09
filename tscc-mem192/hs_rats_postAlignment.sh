#!/bin/bash
#PBS -q hotel
#PBS -N hs_rats_postMap
#PBS -l nodes=tscc-4-4:ppn=1
#PBS -l walltime=8:00:00
#PBS -V
#PBS -j oe
#PBS -k oe
#PBS -M dec037@health.ucsd.edu
#PBS -m abe

cd /home/dec037
#### Make sure to change the compute node, job name, num of node,
#### ppn, walltime, forward email address for notifications. (above)
dir_path=/projects/ps-palmer/hs_rats/20200831/tscc-4-1

out_path=${dir_path}
#### all the output goes to this directory

################# Map sequence reads against reference genome #################
#### Map to reference
START=$(date +%s)

if [ -d "${out_path}/rn6/bams" ]; then
   echo "clean folder: ${out_path}/rn6/bams"
   rm -rf ${out_path}/rn6/bams/*
else
   echo "create folder: ${out_path}/rn6/bams"
   mkdir ${out_path}/rn6/bams
fi

mkdir ${out_path}/rn6/bams/metrics

if [ -d "${out_path}/rn7/bams" ]; then
   echo "clean folder: ${out_path}/rn7/bams"
   rm -rf ${out_path}/rn7/bams/*
else
   echo "create folder: ${out_path}/rn7/bams"
   mkdir ${out_path}/rn7/bams
fi

mkdir ${out_path}/rn7/bams/metrics

END=$(date +%s)
echo "Pre BWA map to reference Time elapsed: $(( $END - $START )) seconds"
