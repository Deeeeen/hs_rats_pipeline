#!/bin/bash
pipeline_arguments=/hs_genotyping_step1/data/pipeline_arguments
dir_path=/hs_genotyping_step1/data/$(head -n 1 ${pipeline_arguments} | tail -n 1)
original_sample_sheet=/hs_genotyping_step1/data/$(head -n 2 ${pipeline_arguments} | tail -n 1)
#### read in arguments for the pipeline

################################### Pre Demux ###################################
#### This part constructs the directory structure.
mkdir ${dir_path}

#### Make directories to keep demux results
if [ -d "${dir_path}/demux" ]; then
   echo "clean folder: ${dir_path}/demux"
   rm -rf ${dir_path}/demux/*
else
   echo "create folder: ${dir_path}/demux"
   mkdir ${dir_path}/demux
fi
mkdir ${dir_path}/demux/fastq
mkdir ${dir_path}/demux/metrics

#### Make directories to keep qc results
if [ -d "${dir_path}/qc" ]; then
   echo "clean folder: ${dir_path}/qc"
   rm -rf ${dir_path}/qc/*
else
   echo "create folder: ${dir_path}/qc"
   mkdir ${dir_path}/qc
fi

#### Make a directories to keep sam files
if [ -d "${dir_path}/sams" ]; then
   echo "clean folder: ${dir_path}/sams"
   rm -rf ${dir_path}/sams/*
else
   echo "create folder: ${dir_path}/sams"
   mkdir ${dir_path}/sams
fi

#### Make a directories to keep bam files
if [ -d "${dir_path}/bams" ]; then
   echo "clean folder: ${dir_path}/bams"
   rm -rf ${dir_path}/bams/*
else
   echo "create folder: ${dir_path}/bams"
   mkdir ${dir_path}/bams
fi
mkdir ${dir_path}/bams/metrics

#### Make a directories to keep stitch results
if [ -d "${dir_path}/stitch" ]; then
   echo "clean folder: ${dir_path}/stitch"
   rm -rf ${dir_path}/stitch/*
else
   echo "create folder: ${dir_path}/stitch"
   mkdir ${dir_path}/stitch
fi

#### Make a directories to keep beagle results
if [ -d "${dir_path}/beagle" ]; then
   echo "clean folder: ${dir_path}/beagle"
   rm -rf ${dir_path}/beagle/*
else
   echo "create folder: ${dir_path}/beagle"
   mkdir ${dir_path}/beagle
fi

#### Make a directories to keep all results
if [ -d "${dir_path}/results" ]; then
   echo "clean folder: ${dir_path}/results"
   rm -rf ${dir_path}/results/*
else
   echo "create folder: ${dir_path}/results"
   mkdir ${dir_path}/results
fi

################################### Pre Demux ###################################
#### This part separates and extracts the big sample sheet that Fgbio needs into
#### several small sample sheets by combination of "pcr_barcode", "library",
#### and "full_run_id"
#### !!!!!!!!!!!!!!!!!!!!!!
#### The following part needs modifications
#### since original sample sheet always
#### comes in with DIFFERENT format.
#### !!!!!!!!!!!!!!!!!!!!!!
START=$(date +%s)
#### This block handles the original sample sheet format.
#### 1. remove double quotes
#### 2. keep only hs rats data
#### 3. duplicate the rfid column to be Sample_ID and Sample_Name
#### 4. inplace the header row with new header row
sed 's/\"//g' ${original_sample_sheet} | \
awk -F',' '{ if($14 =="Heterogenous stock") { print } }' | \
awk 'BEGIN{FS=OFS=","} {$1 = $5 OFS $5 OFS $1} 1' | \
sed '1 i\
Sample_ID,Sample_Name,runid,fastq_files,Library_ID,PCR_Barcode,rfid,Sample_Project,Sample_Barcode,filename,comments,flag,sex,coatcolor,organism,strain
' > ${dir_path}/demux/sample_sheet.csv
sample_sheet=${dir_path}/demux/sample_sheet.csv

#### The for loop separates out the flow cell by Library_ID
#### and create separated sample sheet for each Library_ID
for n in $(cut -d "," -f5,6 < ${sample_sheet} | grep -v "Library_ID" | sort| uniq)
do
    library_id=$(cut -d "," -f1 <<< ${n})
    pcr_barcode=$(cut -d "," -f2 <<< ${n})
    full_run_id=$(echo ${dir_path} | rev | cut -d "/" -f1 | rev)
    if [ -f "${dir_path}/demux/SampleSheet_${pcr_barcode}_${library_id}_${full_run_id}.csv" ]; then
        rm ${dir_path}/demux/SampleSheet_${pcr_barcode}_${library_id}_${full_run_id}.csv
    fi
    head -1 ${sample_sheet} > ${dir_path}/demux/SampleSheet_${pcr_barcode}_${library_id}_${full_run_id}.csv
    awk -v pcr_barcode=$pcr_barcode -v library_id=$library_id -F ',' \
    '{ if($6 == pcr_barcode && $5 == library_id){ print } }' ${sample_sheet} >> ${dir_path}/demux/SampleSheet_${pcr_barcode}_${library_id}_${full_run_id}.csv
    #### extract the corresponding sample barcode metadata
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Pre Demux time elapsed: $(( $END - $START )) seconds"