#!/bin/bash
#PBS -q hotel
#PBS -N hs_rats_postQc
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
#### MultiQC reports
START=$(date +%s)

if [ -d "${out_path}/qc/multiqc_demux" ]; then
   echo "clean folder: ${out_path}/qc/multiqc_demux"
   rm -rf ${out_path}/qc/multiqc_demux/*
else
   echo "create folder: ${out_path}/qc/multiqc_demux"
   mkdir ${out_path}/qc/multiqc_demux
fi

multiqc ${out_path}/qc/fastqc_demux/*.zip -o ${out_path}/qc/multiqc_demux

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "MultiQC time elapsed: $(( $END - $START )) seconds"