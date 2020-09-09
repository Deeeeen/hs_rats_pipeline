#!/bin/bash
#PBS -q hotel
#PBS -N hs_rats_preQc
#PBS -l nodes=tscc-4-4:ppn=6
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
###################### Check quality of sequence reads ########################
#### FastQC reports
START=$(date +%s)
ncpu=5

if [ -d "${out_path}/qc" ]; then
   echo "clean folder: ${out_path}/qc"
   rm -rf ${out_path}/qc/*
else
   echo "create folder: ${out_path}/qc"
   mkdir ${out_path}/qc
fi

mkdir ${out_path}/qc/fastqc_demux

echo "Pre QC time elapsed: $(( $END - $START )) seconds"