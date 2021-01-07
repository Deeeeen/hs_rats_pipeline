#!/bin/bash
#PBS -q home
#PBS -N qc
#PBS -l nodes=1:ppn=6
#PBS -l walltime=8:00:00
#PBS -t 1-5
#PBS -V
#PBS -j oe
#PBS -k oe
#PBS -M dec037@health.ucsd.edu
#PBS -m ae

#### Make sure to change the compute node, job name, num of node,
#### ppn, walltime, forward email address for notifications. (above)
#### !!!!!!!!!!!!!!!!!!!!!!
#### Number of array jobs needs modifications.
#### Here we have 5 array jobs since the flow cell 
#### contains 5 different library preparation.
#### !!!!!!!!!!!!!!!!!!!!!!

pipeline_arguments=$ARG
home=$(head -n 1 ${pipeline_arguments})
dir_path=$(head -n 2 ${pipeline_arguments} | tail -n 1)
demux_data=${dir_path}/demux/fastq
bams_data=${dir_path}/bams
#### read in arguments for the pipeline

cd ${home}

sample_sheet=$(ls ${dir_path}/demux/SampleSheet_*.csv | head -n ${PBS_ARRAYID} | tail -n 1)
samples=$(awk -F',' 'NR>1 {print $1}' ${sample_sheet})
library=$(echo ${sample_sheet} | rev | cut -d '/' -f1 | rev | cut -d '_' -f3)
out_path=${dir_path}/qc/${library}
#### choose a library preparation to process

if [ -d "${out_path}" ]; then
   echo "clean folder: ${out_path}"
   rm -rf ${out_path}/*
else
   echo "create folder: ${out_path}"
   mkdir ${out_path}
fi
mkdir ${out_path}/fastqc_demux
mkdir ${out_path}/qualimap
mkdir ${out_path}/picard

########################## FastQC, Qualimap, MultiQC ###########################
START=$(date +%s)
ncpu=6
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
cnt=0
for sample in ${samples[@]}
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))

   fastqs=$(ls ${demux_data}/${sample}*.fastq.gz)
   for f in ${fastqs}
   do
    echo -e "\n-----run FastQC on ${cnt}-th file: $f-----"
    /projects/ps-palmer/software/local/src/FastQC/fastqc $f --outdir=${out_path}/fastqc_demux/ &
   done
   #### FastQC

   bams=$(ls ${bams_data}/${sample}*_sorted_mkDup.bam)
   for f in ${bams}
   do
    echo -e "\n-----run Qualimap on ${cnt}-th file: $f-----"
    /projects/ps-palmer/software/local/src/qualimap_v2.2.1/qualimap bamqc -bam $f -outdir ${out_path}/qualimap/${sample}
   done
   #### Qualimap

   mkDups=$(ls ${bams_data}/metrics/${sample}*_sorted_mkDup_metrics.txt)
   for mkDup in ${mkDups}
   do
    cp ${mkDup} ${out_path}/picard/${sample}_mkDup_metrics.txt
   done
   #### makeDuplicates metrics
done
multiqc ${out_path}/ -o ${out_path}/multiqc

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "FastQC, Qualimap, MultiQC time elapsed: $(( $END - $START )) seconds"
