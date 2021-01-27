#!/bin/bash

pipeline_arguments=${ARG}
home=$(head -n 1 ${pipeline_arguments})
dir_path=$(head -n 2 ${pipeline_arguments} | tail -n 1)
genetic_map=$(head -n 7 ${pipeline_arguments} | tail -n 1)
stitch_path=${dir_path}/stitch
beagle_path=${dir_path}/beagle
#### read in arguments for the pipeline

cd ${home}

################ imputation STITCH #######################
START=$(date +%s)

if [ "${PBS_ARRAYID}" == "21" ]; then
    ch="X"
else
    ch=${PBS_ARRAYID}
fi

if [ -f "${beagle_path}/chr${ch}_hs_bgl.vcf.gz" ]; then
    rm ${beagle_path}/chr${ch}_hs_bgl*
fi

java -Xss200M -Xmx100G -XX:+AggressiveOpts -XX:+AggressiveHeap \
    -jar /home/dec037/applications/beagle.27Jan18.7e1.jar \
    gt=${stitch_path}/chr${ch}_hs_stitch.vcf.gz \
    gprobs=true \
    ne=${ne} \
    nthreads=${ppn} \
    map=${genetic_map}/chr${ch}_beagle.map \
    out=${beagle_path}/chr${ch}_hs_bgl

END=$(date +%s)
echo "BEAGLE Time elapsed: $(( $END - $START )) seconds"
