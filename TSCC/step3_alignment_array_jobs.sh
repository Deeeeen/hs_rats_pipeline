#!/bin/bash

pipeline_arguments=$ARG
home=$(head -n 1 ${pipeline_arguments})
dir_path=$(head -n 2 ${pipeline_arguments} | tail -n 1)
reference_data=$(head -n 5 ${pipeline_arguments} | tail -n 1)
sample_sheet=${dir_path}/demux/sample_sheet.csv
demux_data=${dir_path}/demux/fastq
#### read in arguments for the pipeline

cd ${home}

################# Map sequence reads against reference genome #################
START=$(date +%s)
#### find all demuxed fastq.gz files and construct a list of their prefix
#### organize a group of fastq.gz files for this array job to process
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

#### a customized command to extract the header of the fastq.gz file
#### this to gather the register group info that bwa mem needs
zhead_py=$(cat <<'EOF'
import sys, gzip
gzf = gzip.GzipFile(sys.argv[1], 'rb')
outFile = sys.stdout.buffer if hasattr(sys.stdout, 'buffer') else sys.stdout
numLines = 0
maxLines = int(sys.argv[2])
for line in gzf:
    if numLines >= maxLines:
        sys.exit(0)
    outFile.write(line)
    numLines += 1
EOF
)
zhead() { python -c "$zhead_py" "$@"; }

ncpu=3
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### the -t flag on bwa mem command specifies the number of threads
#### !!!!!!!!!!!!!!!!!!!!!!
cnt=0
for sample in ${sample_temp}
do
  while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
    sleep 60
  done
  sleep 5
  (( cnt += 1 ))
  #### construct the register group for bwa
  fastq_prefix=$(ls ${demux_data}/${sample}*_R1.fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d '_' -f 1)
  fastq_header=$(zhead ${demux_data}/${fastq_prefix}_R1.fastq.gz 1)
  instrument_name=$(cut -d ':' -f 1 <<< $fastq_header | cut -d '@' -f 2)
  run_id=$(cut -d ':' -f 2 <<< $fastq_header)
  flowcell_id=$(cut -d ':' -f 3 <<< $fastq_header)
  flowcell_lane=$(cut -d ':' -f 4 <<< $fastq_header)
  library_id=$(grep -w "^${sample}" $sample_sheet | cut -d ',' -f 3)
  sample_barcode=$(cut -d '-' -f 3 <<< $fastq_prefix)

  echo -e "\n-----run ${cnt}-th file: ${demux_data}/${fastq_prefix} > ${dir_path}/sams/${fastq_prefix}.sam-----"
  /projects/ps-palmer/software/local/src/bwa-0.7.12/bwa mem -aM -t 2\
  -R "@RG\tID:${instrument_name}.${run_id}.${flowcell_id}.${flowcell_lane}\tLB:${library_id}\tPL:ILLUMINA\tSM:${sample}\tPU:${flowcell_id}.${flowcell_lane}.${sample_barcode}" \
  ${reference_data} ${demux_data}/${fastq_prefix}_R1.fastq.gz \
  ${demux_data}/${fastq_prefix}_R2.fastq.gz > ${dir_path}/sams/${fastq_prefix}.sam &
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

END=$(date +%s)
echo "BWA map to reference Time elapsed: $(( $END - $START )) seconds"
