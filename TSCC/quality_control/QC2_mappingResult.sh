#!/bin/bash

pipeline_arguments=${ARG}
home=$(head -n 1 ${pipeline_arguments})
code=$(head -n 8 ${pipeline_arguments} | tail -n 1)
code=${code}/quality_control/util
dir_path=$(head -n 2 ${pipeline_arguments} | tail -n 1)
sample_sheet=${dir_path}/demux/sample_sheet.csv
reference_data=$(head -n 5 ${pipeline_arguments} | tail -n 1)
bams_path=${dir_path}/bams
out_path=${dir_path}/results
demux_result=${out_path}/demux_result
mapping_result=${out_path}/mapping_result
#### read in arguments for the pipeline

cd ${home}

################################ Demux Results ################################
mkdir ${demux_result}
source activate hs_rats
Rscript ${code}/demuxed_reads.r \
   ${dir_path}/demux/metrics \
   ${sample_sheet} \
   ${demux_result}
conda deactivate

################ Mapping stats and mark Duplicates Results ####################
###########  organize MarkDuplicates metrics
mkdir ${mapping_result}
metrics_dir=${bams_path}/metrics

fs_in=$(ls ${metrics_dir}/*mkDup_metrics.txt)
if [ -f "${mapping_result}/mkDup_metrics" ]; then
   echo "rm file: ${mapping_result}/mkDup_metrics"
   rm ${mapping_result}/mkDup_metrics
fi
echo "create file: ${mapping_result}/mkDup_metrics"
touch ${mapping_result}/mkDup_metrics

cnt=0
for f in $fs_in
do
    (( cnt += 1 ))
    if [ "${cnt}" == "1" ]; then
        header=$(head -n 7 ${f} | tail -n 1)
        echo -e "rfid\t${header}" >> ${mapping_result}/mkDup_metrics
    fi
    content=$(head -n 8 ${f} | tail -n 1)
    sample=$(echo ${f} | rev | cut -d '/' -f 1 | rev | cut -d '-' -f 1)
    echo -e "$sample\t$content" >> ${mapping_result}/mkDup_metrics
done

########### Uniquely mapping
ncpu=${ppn}
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
START=$(date +%s)
if [ -f "${mapping_result}/uniq_mapped" ]; then
   echo "rm file: ${mapping_result}/uniq_mapped"
   rm ${mapping_result}/uniq_mapped
fi
echo "create file: ${mapping_result}/uniq_mapped"
touch ${mapping_result}/uniq_mapped

cnt=0
fs_in=$(ls ${bams_path}/*_mkDup.bam)
for f in ${fs_in}
do
   while [ "$(jobs -rp | wc -l)" -ge ${ncpu} ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   if [ "${cnt}" == "1" ]; then
        echo -e "rfid\tuniq_mapped" >> ${mapping_result}/uniq_mapped
   fi
   sample=$(echo ${f} | rev | cut -d '/' -f 1 | rev | cut -d '-' -f 1)
   content=$(/projects/ps-palmer/software/local/src/samtools-1.10/samtools view -c -q 40 $f)
   echo -e "${sample}\t${content}" >> ${mapping_result}/uniq_mapped
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Uniquely mapping Time elapsed: $(( $END - $START )) seconds"

########### mapping stats plot Rscript
source activate hs_rats
Rscript ${code}/mapped_reads.r \
   ${sample_sheet} \
   ${mapping_result}/mkDup_metrics \
   ${mapping_result}/uniq_mapped \
   ${mapping_result}
conda deactivate

#################### Mapping stats on each chromosome ##########################
########### mapped reads on each chromosome
START=$(date +%s)
fs_in=$(ls ${bams_path}/*_mkDup.bam)
chrs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chrX chrY)
for chr in ${chrs[@]}
do
  if [ -f "${mapping_result}/mapped_${chr}" ]; then
     echo "rm file: ${mapping_result}/mapped_${chr}"
     rm ${mapping_result}/mapped_${chr}
  fi
  echo "create file: ${mapping_result}/mapped_${chr}"
  touch ${mapping_result}/mapped_${chr}
done
ncpu=${ppn}
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
for f in ${fs_in}
do
  while [ "$(jobs -rp | wc -l)" -ge ${ncpu} ]; do
      sleep 60
  done
  sleep 5
  sampleName=$(echo ${f} | rev | cut -d'/' -f 1 | rev | cut -d '-' -f 1)
  fn=$(echo ${f} | rev | cut -d '.' -f 2 | rev)
  /projects/ps-palmer/software/local/src/samtools-1.10/samtools idxstats ${f} | \
    cut -f 1,3,4 > ${fn}.readCount
  sum_reads=$(awk '{sum = sum+$2+$3} END {print sum}' ${fn}.readCount)
  for chr in ${chrs[@]}
  do
    temp_reads=$(awk -v chr=$chr '{if ($1 == chr) print $2;}' ${fn}.readCount)
    echo "${sampleName} ${temp_reads} ${sum_reads}" >> ${mapping_result}/mapped_${chr}
  done
  rm ${fn}.readCount
done
while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60 
done

source activate hs_rats
Rscript ${code}/mapped_reads_chr.r \
   ${mapping_result} \
   ${sample_sheet} \
   ${mapping_result}
conda deactivate 

END=$(date +%s)
echo "Mapping stats on each chr time elapsed: $(( $END - $START )) seconds"

########## mapped reads GC content on each chromosome
START=$(date +%s)
fs_in=$(ls ${bams_path}/*_mkDup.bam)
chrs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chrX chrY)
for chr in ${chrs[@]}
do
  if [ -f "${mapping_result}/mapped_GC_${chr}" ]; then
     echo "rm file: ${mapping_result}/mapped_GC_${chr}"
     rm ${mapping_result}/mapped_GC_${chr}
  fi
  echo "create file: ${mapping_result}/mapped_GC_${chr}"
  touch ${mapping_result}/mapped_GC_${chr}
done

ncpu=${ppn}
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
module load bedtools
for f in ${fs_in}
do
  while [ "$(jobs -rp | wc -l)" -ge ${ncpu} ]; do
      sleep 60
  done
  sleep 5

  sampleName=$(echo ${f} | rev | cut -d'/' -f 1 | rev | cut -d '-' -f 1)
  fn=$(echo ${f} | rev | cut -d '.' -f 2 | rev)
  bedtools bamtobed -i ${f} > ${fn}.bed
  bedtools nuc -fi ${reference_data} -bed ${fn}.bed > ${fn}.nuc
  for chr in ${chrs[@]}
  do
    gc_content=$(awk -v chr=$chr '{if ($1 == chr) {num_C=num_C+$10; num_G=num_G+$11; sum=sum+$15;}} END {print num_C,num_G,sum}' ${fn}.nuc)
    echo "${sampleName} ${gc_content}" >> ${mapping_result}/mapped_GC_${chr}
  done
  rm ${fn}.bed
  rm ${fn}.nuc
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60 
done

END=$(date +%s)
echo "GC content stats on each chr time elapsed: $(( $END - $START )) seconds"

source activate hs_rats
Rscript ${code}/mapped_GC_chr.r \
   ${mapping_result} \
   ${sample_sheet} \
   ${mapping_result}
conda deactivate