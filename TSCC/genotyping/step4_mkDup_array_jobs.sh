#!/bin/bash

pipeline_arguments=${ARG}
home=$(head -n 1 ${pipeline_arguments})
dir_path=$(head -n 2 ${pipeline_arguments} | tail -n 1)
reference_data=$(head -n 5 ${pipeline_arguments} | tail -n 1)
sample_sheet=${dir_path}/demux/sample_sheet.csv
mapped_data=${dir_path}/sams
bams_data=${dir_path}/bams
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
((num_sample_temp=num_sample/num_jobs+1,temp_sample_end=num_sample_temp*PBS_ARRAYID))
if [[ ${temp_sample_end} -gt ${num_sample} ]]; then
  ((num_sample_temp=num_sample_temp-(temp_sample_end-num_sample)))
  if [[ ${num_sample_temp} -le 0 ]]; then
      exit 0
  fi
  sample_temp=$(echo "${sample_IDs}" | head -${temp_sample_end} | tail -${num_sample_temp})
else
  sample_temp=$(echo "${sample_IDs}" | head -${temp_sample_end} | tail -${num_sample_temp})
fi

ncpu=${ppn}
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
cnt=0
for sample in ${sample_temp}
do
   while [ "$(jobs -rp | wc -l)" -ge ${ncpu} ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   sam_prefix=$(ls ${mapped_data}/${sample}*.sam | rev | cut -d'.' -f 2 |cut -d'/' -f 1 | rev)
   #### covert SAM to BAM
   echo -e "\n-----run ${cnt}-th file: ${mapped_data}/${sam_prefix}.sam > ${bams_data}/${sam_prefix}.bam-----"
   /projects/ps-palmer/software/local/src/samtools-1.10/samtools view -h -b \
      -t ${reference_data} -o ${bams_data}/${sam_prefix}.bam ${mapped_data}/${sam_prefix}.sam
   #### sort the BAM file
   echo -e "\n-----run ${cnt}-th file: ${bams_data}/${sam_prefix}.bam > ${bams_data}/${sam_prefix}_sorted.bam-----"
   /projects/ps-palmer/software/local/src/samtools-1.10/samtools sort -m 30G \
      -o ${bams_data}/${sam_prefix}_sorted.bam ${bams_data}/${sam_prefix}.bam
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

for sample in ${sample_temp}
do
   sam_prefix=$(ls ${mapped_data}/${sample}*.sam | rev | cut -d'.' -f 2 |cut -d'/' -f 1 | rev)
  #### remove the unsorted BAM file
  if [ -f "${bams_data}/${sam_prefix}.bam" ]; then
    rm ${bams_data}/${sam_prefix}.bam
  else
    echo -e "ERROR: ${bams_data}/${sam_prefix}.bam not exist"
  fi
done

END=$(date +%s)
echo "SAM -> BAM -> sorted BAM time elapsed: $(( $END - $START )) seconds"

########################### Mark PCR Duplicates ###############################
START=$(date +%s)
((ncpu=ppn/4))
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### Keep it relatively small for more memory pre process on MarkDuplicates
#### !!!!!!!!!!!!!!!!!!!!!!
cnt=0
for sample in ${sample_temp}
do
   while [ "$(jobs -rp | wc -l)" -ge ${ncpu} ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   sam_prefix=$(ls ${mapped_data}/${sample}*.sam | rev | cut -d'.' -f 2 |cut -d'/' -f 1 | rev)
   #### Mark duplicates
   echo -e "\n-----run ${cnt}-th file: ${bams_data}/${sam_prefix}_sorted.bam > ${bams_data}/${sam_prefix}_sorted_mkDup.bam-----"
   /usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin/java -Xmx30G -XX:+AggressiveOpts -XX:+AggressiveHeap\
      -jar /projects/ps-palmer/software/local/src/picard-2.23.3/picard.jar MarkDuplicates \
      --INPUT ${bams_data}/${sam_prefix}_sorted.bam \
      --REMOVE_DUPLICATES false \
      --ASSUME_SORTED true \
      --METRICS_FILE ${bams_data}/metrics/${sam_prefix}_sorted_mkDup_metrics.txt \
      --OUTPUT ${bams_data}/${sam_prefix}_sorted_mkDup.bam &
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

END=$(date +%s)
echo "Mark PCR duplicates time elapsed: $(( $END - $START )) seconds"

############### Index sorted alignment for fast random access ##################
START=$(date +%s)
ncpu=${ppn}
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
cnt=0
for sample in ${sample_temp}
do
   while [ "$(jobs -rp | wc -l)" -ge ${ncpu} ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   sam_prefix=$(ls ${mapped_data}/${sample}*.sam | rev | cut -d'.' -f 2 |cut -d'/' -f 1 | rev)
   echo -e "\n----- run ${cnt}-th file: ${bams_data}/${sam_prefix}_sorted_mkDup.bam > ${bams_data}/${sam_prefix}_sorted_mkDup.bai-----"
   /projects/ps-palmer/software/local/src/samtools-1.10/samtools index ${bams_data}/${sam_prefix}_sorted_mkDup.bam ${bams_data}/${sam_prefix}_sorted_mkDup.bai
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Index alignments time elapsed: $(( $END - $START )) seconds"

############################# clean up directory ###############################
for sample in ${sample_temp}
do
   sam_prefix=$(ls ${mapped_data}/${sample}*.sam | rev | cut -d'.' -f 2 |cut -d'/' -f 1 | rev)
   #### clean up directory
   if [ -f "${mapped_data}/${sam_prefix}.sam" ] && [ -f "${bams_data}/${sam_prefix}_sorted.bam" ] && [ -f "${bams_data}/${sam_prefix}_sorted_mkDup.bam" ] && [ -f "${bams_data}/${sam_prefix}_sorted_mkDup.bai" ]; then
      rm ${mapped_data}/${sam_prefix}.sam
      rm ${bams_data}/${sam_prefix}_sorted.bam
   else
      echo -e "ERROR: something went wrong during SAM->BAM, BAM->sorted_BAM, or sorted_BAM->sorted_mkDup_BAM for ${sam_prefix}"
   fi
done