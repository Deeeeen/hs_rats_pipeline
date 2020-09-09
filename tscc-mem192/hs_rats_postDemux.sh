#!/bin/bash
#PBS -q hotel
#PBS -N hs_rats_postDemux
#PBS -l nodes=tscc-4-2:ppn=3
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
#### where you keep all the input and output files
out_path=${dir_path}
#### all the output goes to this directory

#### Demux stats plots (# of reads per sample, demux unmatched reads)
START=$(date +%s)
module load R

Rscript ${dir_path}/code/demux_reads.r \
    ${dir_path}/demux/metrics \
    ${out_path}/demux/metrics/demux_reads.png \
    ${out_path}/demux/metrics/demux_unmatched.png

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

END=$(date +%s)
module unload R
echo "Post Demux stats plots time elapsed: $(( $END - $START )) seconds"