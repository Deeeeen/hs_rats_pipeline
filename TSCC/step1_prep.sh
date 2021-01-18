#!/bin/bash
#PBS -q condo
#PBS -N preDemux
#PBS -l nodes=1:ppn=2
#PBS -l walltime=8:00:00
#PBS -V
#PBS -j oe
#PBS -k oe
#PBS -M dec037@health.ucsd.edu
#PBS -m ae

#### Make sure to change the compute node, job name, num of node,
#### ppn, walltime, forward email address for notifications. (above)

pipeline_arguments=$ARG
home=$(head -n 1 ${pipeline_arguments})
dir_path=$(head -n 2 ${pipeline_arguments} | tail -n 1)
code=$(head -n 8 ${pipeline_arguments} | tail -n 1)
original_sample_sheet=$(head -n 3 ${pipeline_arguments} | tail -n 1)
#### read in arguments for the pipeline

cd ${home}
################################### Pre Demux ###################################
#### This part constructs the directory structure.
mkdir ${dir_path}

#### Make directories to keep demux results
if [ -d "${dir_path}/demux" ]; then
   echo "clean folder: ${dir_path}/demux"
   rm -rf ${dir_path}/demux/*
else
   echo "create folder: ${dir_path}/demux"
   mkdir ${dir_path}/demux
fi
mkdir ${dir_path}/demux/fastq
mkdir ${dir_path}/demux/metrics

#### Make directories to keep qc results
if [ -d "${dir_path}/qc" ]; then
   echo "clean folder: ${dir_path}/qc"
   rm -rf ${dir_path}/qc/*
else
   echo "create folder: ${dir_path}/qc"
   mkdir ${dir_path}/qc
fi

#### Make a directories to keep sam files
if [ -d "${dir_path}/sams" ]; then
   echo "clean folder: ${dir_path}/sams"
   rm -rf ${dir_path}/sams/*
else
   echo "create folder: ${dir_path}/sams"
   mkdir ${dir_path}/sams
fi

#### Make a directories to keep bam files
if [ -d "${dir_path}/bams" ]; then
   echo "clean folder: ${dir_path}/bams"
   rm -rf ${dir_path}/bams/*
else
   echo "create folder: ${dir_path}/bams"
   mkdir ${dir_path}/bams
fi
mkdir ${dir_path}/bams/metrics

#### Make a directories to keep stitch results
if [ -d "${dir_path}/stitch" ]; then
   echo "clean folder: ${dir_path}/stitch"
   rm -rf ${dir_path}/stitch/*
else
   echo "create folder: ${dir_path}/stitch"
   mkdir ${dir_path}/stitch
fi

#### Make a directories to keep beagle results
if [ -d "${dir_path}/beagle" ]; then
   echo "clean folder: ${dir_path}/beagle"
   rm -rf ${dir_path}/beagle/*
else
   echo "create folder: ${dir_path}/beagle"
   mkdir ${dir_path}/beagle
fi

#### Make a directories to keep all results
if [ -d "${dir_path}/results" ]; then
   echo "clean folder: ${dir_path}/results"
   rm -rf ${dir_path}/results/*
else
   echo "create folder: ${dir_path}/results"
   mkdir ${dir_path}/results
fi

################################### Pre Demux ###################################
#### This part separates and extracts the big sample sheet that Fgbio needs into
#### several small sample sheets by combination of "pcr_barcode", "library",
#### and "full_run_id"
#### !!!!!!!!!!!!!!!!!!!!!!
#### The following part needs modifications
#### since original sample sheet always
#### comes in with DIFFERENT format.
#### !!!!!!!!!!!!!!!!!!!!!!
START=$(date +%s)
#### This block handles the original sample sheet format.
#### extract the corresponding sample barcode metadata
source activate hs_rats
python3 ${code}/separate_metadata.py \
    ${original_sample_sheet} \
    ${dir_path}/demux
conda deactivate

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Pre Demux time elapsed: $(( $END - $START )) seconds"