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
mapped_data=${dir_path}/sams
#### read in arguments for the pipeline

cd ${home}

################# Convert SAM to BAM, sort BAM by coordinates #################
START=$(date +%s)
#### find all aligned sam files and construct a list of their prefix
fs_in=$(ls ${mapped_data}/*.sam)
sam_prefix=()
cnt=0
for f in $fs_in
do
   (( cnt += 1 ))
   temp=$(echo ${f} | rev | cut -d '/' -f 1 | rev)
   sam_prefix[cnt]=$(echo ${temp} | rev | cut -d '.' -f 2 | rev)
done
sam_prefix=$(for f in ${sam_prefix[@]}; do echo $f; done | sort -u)

#### organize a group of SAM files for this array job to process
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the number 24 based on the number of array jobs requested
#### !!!!!!!!!!!!!!!!!!!!!!
num_sam=0
for f in ${sam_prefix}; 
do 
  (( num_sam += 1))
done
((temp=num_sam/24+1,num_sam_temp=temp*PBS_ARRAYID))
if [[ ${num_sam_temp} -gt ${num_sam} ]]; then
  ((temp=temp-(num_sam_temp-num_sam)))
  if [[ ${temp} -le 0 ]]; then
      exit 0
  fi
  sam_temp=$(echo "${sam_prefix}" | head -$num_sam_temp | tail -$temp)
else
  sam_temp=$(echo "${sam_prefix}" | head -$num_sam_temp | tail -$temp)
fi

ncpu=12
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
cnt=0
for f in ${sam_temp}
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   #### covert SAM to BAM
   echo -e "\n-----run ${cnt}-th file: ${mapped_data}/${f}.sam > ${dir_path}/bams/${f}.bam-----"
   /projects/ps-palmer/software/local/src/samtools-1.10/samtools view -h -b \
      -t ${reference_data} -o ${dir_path}/bams/${f}.bam ${mapped_data}/${f}.sam
   #### sort the BAM file
   echo -e "\n-----run ${cnt}-th file: ${dir_path}/bams/${f}.bam > ${dir_path}/bams/${f}_sorted.bam-----"
   /projects/ps-palmer/software/local/src/samtools-1.10/samtools sort -m 30G \
      -o ${dir_path}/bams/${f}_sorted.bam ${dir_path}/bams/${f}.bam
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

for f in ${sam_temp}
do
  #### remove the unsorted BAM file
  if [ -f "${dir_path}/bams/${f}.bam" ]; then
    rm ${dir_path}/bams/${f}.bam
  else
    echo -e "ERROR: ${dir_path}/bams/${f}.bam not exist"
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
for f in ${sam_temp}
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   #### Mark duplicates
   echo -e "\n-----run ${cnt}-th file: ${dir_path}/bams/${f}_sorted.bam > ${dir_path}/bams/${f}_sorted_mkDup.bam-----"
   java -Xmx20G -XX:+AggressiveOpts -XX:+AggressiveHeap\
      -jar /projects/ps-palmer/software/local/src/picard-2.23.3/picard.jar MarkDuplicates \
      --INPUT ${dir_path}/bams/${f}_sorted.bam \
      --REMOVE_DUPLICATES false \
      --ASSUME_SORTED true \
      --METRICS_FILE ${dir_path}/bams/metrics/${f}_sorted_mkDup_metrics.txt \
      --OUTPUT ${dir_path}/bams/${f}_sorted_mkDup.bam &
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
for f in ${sam_temp}
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   echo -e "\n----- run ${cnt}-th file: ${dir_path}/bams/${f}_sorted_mkDup.bam > ${dir_path}/bams/${f}_sorted_mkDup.bai-----"
   /projects/ps-palmer/software/local/src/samtools-1.10/samtools index ${dir_path}/bams/${f}_sorted_mkDup.bam ${dir_path}/bams/${f}_sorted_mkDup.bai
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Index alignments time elapsed: $(( $END - $START )) seconds"

############################# clean up directory ###############################
for f in ${sam_temp}
do
  if [ -f "${mapped_data}/${f}.sam" ] && [ -f "${dir_path}/bams/${f}_sorted.bam" ] && [ -f "${dir_path}/bams/${f}_sorted_mkDup.bam" ] && [ -f "${dir_path}/bams/${f}_sorted_mkDup.bai" ]; then
    rm ${mapped_data}/${f}.sam
    rm ${dir_path}/bams/${f}_sorted.bam
  else
    echo -e "ERROR: something went wrong during SAM->BAM, BAM->sorted_BAM, or sorted_BAM->sorted_mkDup_BAM fro ${f}"
  fi
done