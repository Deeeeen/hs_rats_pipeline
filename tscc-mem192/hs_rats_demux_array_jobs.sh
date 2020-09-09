#!/bin/bash
#PBS -q hotel
#PBS -N hs_rats_demux
#PBS -l nodes=tscc-4-2:ppn=6
#PBS -l walltime=8:00:00
#PBS -t 1-7
#PBS -V
#PBS -j oe
#PBS -k oe
#PBS -M dec037@health.ucsd.edu
#PBS -m abe
START=$(date +%s)
cd /home/dec037
#### Make sure to change the compute node, job name, num of node,
#### ppn, walltime, forward email address for notifications. (above)

dir_path=/projects/ps-palmer/hs_rats/20200831/tscc-4-1
#### where you keep all the input and output files
out_path=${dir_path}
#### all the output goes to this directory

pre_demux_fastq=`head -$PBS_ARRAYID ${dir_path}/code/preDemux_fastq_list | tail -1`
full_run_id=$(rev <<< ${pre_demux_fastq}| cut -d'/' -f 2 | rev)
fastq_temp=$(rev <<< ${pre_demux_fastq}| cut -d'/' -f 1 | rev)

sample_sheets=$(ls ${dir_path}/demux/SampleSheet_*.csv)
sample_sheet=false
for sample_sheet_temp in ${sample_sheets}
do
  sample_sheet_n=$(rev <<< ${sample_sheet_temp}| cut -d'/' -f 1 | rev)
  if [[ "${sample_sheet_n}" == *"${full_run_id}"* ]]; then
    if [[ "${fastq_temp}" == *"Riptide"* ]]; then
      library_id=$(cut -d'_' -f 3 <<< ${sample_sheet_n})
      if [[ "${fastq_temp}" == *"${library_id}"* ]]; then
        sample_sheet=${sample_sheet_temp}
        break
      fi
    else
      pcr_barcode=$(cut -d'_' -f 2 <<< ${sample_sheet_n})
      if [[ "${fastq_temp}" == *"_S${pcr_barcode}_"* ]]; then
        sample_sheet=${sample_sheet_temp}
        break
      fi
    fi
  fi
done
#### to find the corresponding sample sheet file

if [ "${sample_sheet}" = false ]; then
  echo -e "No corresponding sample sheet for ${pre_demux_fastq}."
  exit 1
fi

if [[ -f "${pre_demux_fastq}_R1_001.fastq.gz" && -f "${pre_demux_fastq}_R2_001.fastq.gz" ]]; then
  echo "\n-----run demux on ${pre_demux_fastq} with ${sample_sheet}-----"
  ncpu=4
  java -Xmx40G -XX:+AggressiveOpts -XX:+AggressiveHeap \
       -jar /projects/ps-palmer/software/local/src/fgbio-1.2.0/fgbio-1.2.0.jar DemuxFastqs \
       --inputs ${pre_demux_fastq}_R1_001.fastq.gz \
                ${pre_demux_fastq}_R2_001.fastq.gz \
       --metadata ${sample_sheet} \
       --read-structures 8B12M+T 8M+T \
       --output-type=Fastq \
       --threads $ncpu \
       --output ${out_path}/demux/fastq \
       --metrics ${out_path}/demux/metrics/${fastq_temp}demux_barcode_metrics.txt
  #### Use the path where you keep fgbio after flag -jar
  #### -Xmx40G memory request for JVM
else
    echo -e "\n----- sequence data ${pre_demux_fastq}_R1_001.fastq.gz or ${pre_demux_fastq}_R2_001.fastq.gz don't exist-----"
fi

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Demux time elapsed: $(( $END - $START )) seconds"
#### START and END keep track of how much time were used for thie process.
#### The while loop makes sure all steps in this process are done before
#### entering next process.