#!/bin/bash

pipeline_arguments=/pipeline/data/pipeline_arguments
dir_path=/pipeline/data/$(head -n 1 ${pipeline_arguments} | tail -n 1)
fastq_dir=$(head -n 3 ${pipeline_arguments} | tail -n 1)
#### read in arguments for the pipeline

cd 

################################### Demux ###################################
#### !!!!!!!!!!!!!!!!!!!!!!
#### The following part may need modifications
#### since original sample sheet always
#### comes in with DIFFERENT format.
#### !!!!!!!!!!!!!!!!!!!!!!
START=$(date +%s)

sample_sheet=$(ls ${dir_path}/demux/SampleSheet_*.csv | head -n ${JOB_NUM} | tail -n 1)
flow_cell=$(echo ${sample_sheet} | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | cut -d '_' -f4-)
pre_demux_fastqs=$(head -n 2 ${sample_sheet} | tail -n 1 | cut -d ',' -f6)
pre_demux_fastq_R1=$(cut -d ';' -f1 <<< ${pre_demux_fastqs})
pre_demux_fastq_R2=$(cut -d ';' -f2 <<< ${pre_demux_fastqs} | sed 's/^ *//g')
metrics_name=$(echo ${sample_sheet} | rev | cut -d '/' -f1 | rev | cut -d '.' -f1)

echo "\n-----run demux on ${fastq_dir}/${flow_cell}/${pre_demux_fastq_R1} with ${sample_sheet}-----"
ncpu=4
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
#### Use the path where you keep fgbio after flag -jar
#### -Xmx40G memory request for JVM
java -Xmx40G -XX:+AggressiveOpts -XX:+AggressiveHeap \
     -jar /pipeline/software/fgbio-1.3.0.jar DemuxFastqs \
     --inputs ${fastq_dir}/${flow_cell}/${pre_demux_fastq_R1} \
              ${fastq_dir}/${flow_cell}/${pre_demux_fastq_R2} \
     --metadata ${sample_sheet} \
     --read-structures 8B12M+T 8M+T \
     --output-type=Fastq \
     --threads $ncpu \
     --output ${dir_path}/demux/fastq \
     --metrics ${dir_path}/demux/metrics/${metrics_name}_demux_barcode_metrics.txt

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Demux time elapsed: $(( $END - $START )) seconds"
#### START and END keep track of how much time were used for thie process.
#### The while loop makes sure all steps in this process are done before
#### entering next process.