#!/bin/bash
#PBS -q hotel
#PBS -N hs_rats_qc
#PBS -l nodes=tscc-4-4:ppn=6
#PBS -l walltime=8:00:00
#PBS -t 1-10
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
demuxed_data=${out_path}/demux/fastq
###################### Check quality of sequence reads ########################
#### FastQC reports
START=$(date +%s)
ncpu=5

num_fastq=$(find ${demuxed_data}/*.fastq.gz ! -name '*unmatched*' | wc -l)
((temp=num_fastq/10+1,num_fastq_temp=temp*PBS_ARRAYID))
if [[ ${num_fastq_temp} -gt ${num_fastq} ]]; then
  ((temp=temp-(num_fastq_temp-num_fastq)))
  if [[ ${temp} -le 0 ]]; then
      exit 0
  fi
  fastq_temp=$(find ${demuxed_data}/*.fastq.gz ! -name '*unmatched*' | head -$num_fastq_temp | tail -$temp)
else
  fastq_temp=$(find ${demuxed_data}/*.fastq.gz ! -name '*unmatched*' | head -$num_fastq_temp | tail -$temp)
fi

cnt=0
for f in ${fastq_temp}
do
   while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))
   echo -e "\n-----run FastQC on ${cnt}-th file: $f-----"
   /projects/ps-palmer/software/local/src/FastQC/fastqc $f --outdir=${out_path}/qc/fastqc_demux/ &
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "FastQC time elapsed: $(( $END - $START )) seconds"