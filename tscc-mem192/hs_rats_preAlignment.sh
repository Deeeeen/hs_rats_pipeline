#!/bin/bash
#PBS -q hotel
#PBS -N hs_rats_preMap
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

if [ -d "${out_path}/rn7" ]; then
   echo "clean folder: ${out_path}/rn7"
   rm -rf ${out_path}/rn7/*
else
   echo "create folder: ${out_path}/rn7"
   mkdir ${out_path}/rn7
fi

mkdir ${out_path}/rn7/sams

if [ -d "${out_path}/rn6" ]; then
   echo "clean folder: ${out_path}/rn6"
   rm -rf ${out_path}/rn6/*
else
   echo "create folder: ${out_path}/rn6"
   mkdir ${out_path}/rn6
fi

mkdir ${out_path}/rn6/sams


END=$(date +%s)
echo "Pre BWA map to reference Time elapsed: $(( $END - $START )) seconds"
