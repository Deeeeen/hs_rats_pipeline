#!/bin/bash
#PBS -q condo
#PBS -N markDup
#PBS -l nodes=1:ppn=12
#PBS -l walltime=8:00:00
#PBS -t 1-24%5
#PBS -V
#PBS -j oe
#PBS -k oe
#PBS -M dec037@health.ucsd.edu
#PBS -m ae

#### Make sure to change the compute node, job name, num of node,
#### ppn, walltime, forward email address for notifications. (above)
#### !!!!!!!!!!!!!!!!!!!!!!
#### Number of array jobs needs modifications.
#### Here we have 24 array jobs
#### !!!!!!!!!!!!!!!!!!!!!!

pipeline_arguments=pipeline_arguments
home=$(head -n 1 ${pipeline_arguments})
dir_path=$(head -n 2 ${pipeline_arguments} | tail -n 1)
reference_data=$(head -n 5 ${pipeline_arguments} | tail -n 1)
sample_sheet=${dir_path}/demux/sample_sheet.csv
mapped_data=${dir_path}/sams
#### read in arguments for the pipeline

cd ${home}

################# Convert SAM to BAM, sort BAM by coordinates #################
START=$(date +%s)
#### find all aligned sam files and construct a list of their prefix
#### organize a group of SAM files for this array job to process
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the number 24 based on the number of array jobs requested
#### !!!!!!!!!!!!!!!!!!!!!!
sample_IDs=$(awk -F ',' 'NR>1 {print $1}' ${sample_sheet})
num_sample=$(awk -F ',' 'NR>1 {print $1}' ${sample_sheet} | wc -l)
((num_sample_temp=num_sample/24+1,temp_sample_end=num_sample_temp*PBS_ARRAYID))
if [[ ${temp_sample_end} -gt ${num_sample} ]]; then
  ((num_sample_temp=num_sample_temp-(temp_sample_end-num_sample)))
  if [[ ${num_sample_temp} -le 0 ]]; then
      exit 0
  fi
  sample_temp=$(echo "${sample_IDs}" | head -$temp_sample_end | tail -$num_sample_temp)
else
  sample_temp=$(echo "${sample_IDs}" | head -$temp_sample_end | tail -$num_sample_temp)
fi

ncpu=12
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
cnt=0
for sample in ${sample_temp}
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   sam_prefix=$(ls ${mapped_data}/${sample}*.sam | rev | cut -d'/' -f 1 | rev | cut -d '_' -f 1)
   #### covert SAM to BAM
   echo -e "\n-----run ${cnt}-th file: ${mapped_data}/${sam_prefix}.sam > ${dir_path}/bams/${sam_prefix}.bam-----"
   /projects/ps-palmer/software/local/src/samtools-1.10/samtools view -h -b \
      -t ${reference_data} -o ${dir_path}/bams/${sam_prefix}.bam ${mapped_data}/${sam_prefix}.sam
   #### sort the BAM file
   echo -e "\n-----run ${cnt}-th file: ${dir_path}/bams/${sam_prefix}.bam > ${dir_path}/bams/${sam_prefix}_sorted.bam-----"
   /projects/ps-palmer/software/local/src/samtools-1.10/samtools sort -m 30G \
      -o ${dir_path}/bams/${sam_prefix}_sorted.bam ${dir_path}/bams/${sam_prefix}.bam
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

for sample in ${sample_temp}
do
  sam_prefix=$(ls ${mapped_data}/${sample}*.sam | rev | cut -d'/' -f 1 | rev | cut -d '_' -f 1)
  #### remove the unsorted BAM file
  if [ -f "${dir_path}/bams/${sam_prefix}.bam" ]; then
    rm ${dir_path}/bams/${sam_prefix}.bam
  else
    echo -e "ERROR: ${dir_path}/bams/${sam_prefix}.bam not exist"
  fi
done

END=$(date +%s)
echo "SAM -> BAM -> sorted BAM time elapsed: $(( $END - $START )) seconds"

########################### Mark PCR Duplicates ###############################
START=$(date +%s)
ncpu=3
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### Keep it relatively small for more memory pre process on MarkDuplicates
#### !!!!!!!!!!!!!!!!!!!!!!
cnt=0
for sample in ${sample_temp}
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   sam_prefix=$(ls ${mapped_data}/${sample}*.sam | rev | cut -d'/' -f 1 | rev | cut -d '_' -f 1)
   #### Mark duplicates
   echo -e "\n-----run ${cnt}-th file: ${dir_path}/bams/${sam_prefix}_sorted.bam > ${dir_path}/bams/${sam_prefix}_sorted_mkDup.bam-----"
   java -Xmx20G -XX:+AggressiveOpts -XX:+AggressiveHeap\
      -jar /projects/ps-palmer/software/local/src/picard-2.23.3/picard.jar MarkDuplicates \
      --INPUT ${dir_path}/bams/${sam_prefix}_sorted.bam \
      --REMOVE_DUPLICATES false \
      --ASSUME_SORTED true \
      --METRICS_FILE ${dir_path}/bams/metrics/${sam_prefix}_sorted_mkDup_metrics.txt \
      --OUTPUT ${dir_path}/bams/${sam_prefix}_sorted_mkDup.bam &
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

END=$(date +%s)
echo "Mark PCR duplicates time elapsed: $(( $END - $START )) seconds"

############### Index sorted alignment for fast random access ##################
START=$(date +%s)
ncpu=12
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
cnt=0
for sample in ${sample_temp}
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   sam_prefix=$(ls ${mapped_data}/${sample}*.sam | rev | cut -d'/' -f 1 | rev | cut -d '_' -f 1)
   echo -e "\n----- run ${cnt}-th file: ${dir_path}/bams/${sam_prefix}_sorted_mkDup.bam > ${dir_path}/bams/${sam_prefix}_sorted_mkDup.bai-----"
   /projects/ps-palmer/software/local/src/samtools-1.10/samtools index ${dir_path}/bams/${sam_prefix}_sorted_mkDup.bam ${dir_path}/bams/${sam_prefix}_sorted_mkDup.bai
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Index alignments time elapsed: $(( $END - $START )) seconds"

############################# clean up directory ###############################
for sample in ${sample_temp}
do
  sam_prefix=$(ls ${mapped_data}/${sample}*.sam | rev | cut -d'/' -f 1 | rev | cut -d '_' -f 1)
  #### clean up directory
  if [ -f "${mapped_data}/${sam_prefix}.sam" ] && [ -f "${dir_path}/bams/${sam_prefix}_sorted.bam" ] && [ -f "${dir_path}/bams/${sam_prefix}_sorted_mkDup.bam" ] && [ -f "${dir_path}/bams/${sam_prefix}_sorted_mkDup.bai" ]; then
    rm ${mapped_data}/${sam_prefix}.sam
    rm ${dir_path}/bams/${sam_prefix}_sorted.bam
  else
    echo -e "ERROR: something went wrong during SAM->BAM, BAM->sorted_BAM, or sorted_BAM->sorted_mkDup_BAM fro ${sam_prefix}"
  fi
done