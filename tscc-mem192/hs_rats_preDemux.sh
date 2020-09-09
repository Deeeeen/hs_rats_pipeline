#!/bin/bash
#PBS -q hotel
#PBS -N hs_rats_preDemux
#PBS -l nodes=tscc-4-2:ppn=3
#PBS -l walltime=168:00:00
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
sample_barcode_lib=/projects/ps-palmer/hs_rats/20200831/data/sequencing01.csv
#### .csv sample barcode metadata
out_path=${dir_path}
#### all the output goes to this directory

################################### Pre Demux ###################################
#### This part separate and extract the big sample sheet that Fgbio needs into
#### several small sample sheets by combination of "pcr_barcode", "library", 
#### and "full_run_id"
START=$(date +%s)

if [ -d "${out_path}/demux" ]; then
   echo "clean folder: ${out_path}/demux"
   rm -rf ${out_path}/demux/*
else
   echo "create folder: ${out_path}/demux"
   mkdir ${out_path}/demux
fi
#### Make a directory to keep all demux results

if [ -d "${out_path}/demux/fastq" ]; then
   echo "clean folder: ${out_path}/demux/fastq"
   rm -rf ${out_path}/demux/fastq/*
else
   echo "create folder: ${out_path}/demux/fastq"
   mkdir ${out_path}/demux/fastq
fi
#### Make a directory to keep all demux fastq results

if [ -d "${out_path}/demux/metrics" ]; then
   echo "clean folder: ${out_path}/demux/metrics"
   rm -rf ${out_path}/demux/metrics/*
else
   echo "create folder: ${out_path}/demux/metrics"
   mkdir ${out_path}/demux/metrics
fi
#### Make a directory to keep all demux stats results

sed 's/\([^,]*\),\(.*\)/\2/' ${sample_barcode_lib} | \
sed 's/\"//g' | \
awk 'BEGIN{FS=OFS=","} {$1 = $1 OFS $1} 1' | \
sed '1 c\
Sample_ID,Sample_Name,Sample_Project,PCR_Barcode,Sample_Barcode,Library_ID,Description,flag,full_run_id
' > ${out_path}/demux/sample_barcode_lib.csv
#### remove the first column, remove double quotation marks,
#### duplicate the rfid as sample name, rename the header
sample_barcode_lib=${out_path}/demux/sample_barcode_lib.csv
#### Sample barcode lib metadata modification for better fit in Fgbio

cnt=0
for n in $(cut -d "," -f4,6,9 < ${sample_barcode_lib} |grep -v "full_run_id" | uniq )
do
    (( cnt += 1 ))
    pcr_barcode=$(cut -d "," -f1 <<< ${n})
    library_id=$(cut -d "," -f2 <<< ${n})
    full_run_id=$(cut -d "," -f3 <<< ${n})
    if [ -f "${out_path}/demux/SampleSheet_${pcr_barcode}_${library_id}_${full_run_id}.csv" ]; then
        rm ${out_path}/demux/SampleSheet_${pcr_barcode}_${library_id}_${full_run_id}.csv
    fi
    head -1 $sample_barcode_lib > ${out_path}/demux/SampleSheet_${pcr_barcode}_${library_id}_${full_run_id}.csv
    awk -v pcr_barcode=$pcr_barcode -v full_run_id=$full_run_id -v library_id=$library_id -F ',' \
    '{ if($4 == pcr_barcode && $6 == library_id && $9 == full_run_id){ print } }' $sample_barcode_lib >> ${out_path}/demux/SampleSheet_${pcr_barcode}_${library_id}_${full_run_id}.csv
    sample_barcode_lib_f=${out_path}/demux/SampleSheet_${pcr_barcode}_${library_id}_${full_run_id}.csv
    #### extract the corresponding sample barcode metadata
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Pre Demux time elapsed: $(( $END - $START )) seconds"