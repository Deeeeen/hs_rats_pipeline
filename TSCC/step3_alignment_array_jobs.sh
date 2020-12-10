#!/bin/bash
#PBS -q condo
#PBS -N mapping
#PBS -l nodes=1:ppn=6
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
sample_sheet=${dir_path}/demux/sample_barcode_lib.csv
demux_data=${dir_path}/demux/fastq
#### read in arguments for the pipeline

cd ${home}

################# Map sequence reads against reference genome #################
START=$(date +%s)
#### find all demuxed fastq.gz files and construct a list of their prefix
fastq_prefix=()
fastq_fs=$(find ${demux_data}/*.fastq.gz ! -name '*unmatched*')
cnt=0
for f in ${fastq_fs}
do
   (( cnt += 1 ))
   temp=$(echo ${f} | rev | cut -d'/' -f 1 | rev)
   fastq_prefix[cnt]=$(cut -d '_' -f 1 <<< $temp)
done
fastq_prefix=$(for f in ${fastq_prefix[@]}; do echo $f; done | sort -u)

#### organize a group of fastq.gz files for this array job to process
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the number 24 based on the number of array jobs requested
#### !!!!!!!!!!!!!!!!!!!!!!
num_fastq=0
for f in ${fastq_prefix}; 
do 
  (( num_fastq += 1))
done
((temp=num_fastq/24+1,num_fastq_temp=temp*PBS_ARRAYID))
if [[ ${num_fastq_temp} -gt ${num_fastq} ]]; then
  ((temp=temp-(num_fastq_temp-num_fastq)))
  if [[ ${temp} -le 0 ]]; then
      exit 0
  fi
  fastq_temp=$(echo "${fastq_prefix}" | head -$num_fastq_temp | tail -$temp)
else
  fastq_temp=$(echo "${fastq_prefix}" | head -$num_fastq_temp | tail -$temp)
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
for f in ${fastq_temp}
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   #### construct the register group for bwa
   fastq_header=$(zhead ${demux_data}/${f}_R1.fastq.gz 1)
   instrument_name=$(cut -d ':' -f 1 <<< $fastq_header | cut -d '@' -f 2)
   run_id=$(cut -d ':' -f 2 <<< $fastq_header)
   flowcell_id=$(cut -d ':' -f 3 <<< $fastq_header)
   flowcell_lane=$(cut -d ':' -f 4 <<< $fastq_header)
   sample_id=$(cut -d '-' -f 1 <<< $f)
   library_id=$(grep -w "^${sample_id}" $sample_barcode_lib | cut -d ',' -f 5)
   sample_barcode=$(cut -d '-' -f 3 <<< $f)
   RG="@RG\tID:${instrument_name}.${run_id}.${flowcell_id}.${flowcell_lane}\tLB:${library_id}\tPL:ILLUMINA\tSM:${sample_id}\tPU:${flowcell_id}.${flowcell_lane}.${sample_barcode}"

   echo -e "\n-----run ${cnt}-th file: ${demux_data}/${f} > ${dir_path}/sams/${f}.sam-----"
   /projects/ps-palmer/software/local/src/bwa-0.7.12/bwa mem -aM -t 2\
     -R ${RG} \
     ${reference_data} ${demux_data}/${f}_R1.fastq.gz \
     ${demux_data}/${f}_R2.fastq.gz > ${dir_path}/sams/${f}.sam &
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

END=$(date +%s)
echo "BWA map to reference Time elapsed: $(( $END - $START )) seconds"
