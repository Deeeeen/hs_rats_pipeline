#!/bin/bash

ARG=/projects/ps-palmer/hs_rats/201218_A00953_0203_BHN5HYDSXY/code/pipeline_arguments
METADATA=/projects/ps-palmer/hs_rats/201218_A00953_0203_BHN5HYDSXY/code/previous_flow_cells_metadata
PEDIGREE=/projects/ps-palmer/hs_rats/201218_A00953_0203_BHN5HYDSXY/code/pedigree_data
PREV_BAMS=/projects/ps-palmer/hs_rats/201218_A00953_0203_BHN5HYDSXY/code/previous_flow_cells_bams
code=$(head -n 8 ${ARG} | tail -n 1)

email=dec037@health.ucsd.edu

#### Make sure to change the compute node, job name, num of node,
#### ppn, walltime, forward email address for notifications. (below)

#### Genotyping pipeline
STEP1_PREP_JOB=$(qsub -q condo -N preDemux -l nodes=1:ppn=2,walltime=8:00:00 \
                      -j oe -k oe -m oe -M ${email} \
                      -V -v ARG="$ARG" \
                      ${code}/step1_prep.sh)
echo "step1_prep: $STEP1_PREP_JOB"

STEP2_DEMUX_JOB_ARR=$(qsub -q hotel -N demux -l nodes=1:ppn=4,walltime=24:00:00 -t 1-5 \
                           -j oe -k oe -m oe -M ${email} \
                           -V -v ARG="$ARG" \
                           -W depend=afterok:$STEP1_PREP_JOB \
                           ${code}/step2_demux_array_jobs.sh)
echo "step2_demux: $STEP2_DEMUX_JOB_ARR"
STEP2_DEMUX_JOB_ARR_id=$(echo "${STEP2_DEMUX_JOB_ARR}" | cut -d '.' -f 1 )
#### !!!!!!!!!!!!!!!!!!!!!!
#### Number of array jobs needs modifications.
#### Here we have 5 array jobs since the last flow cell 
#### contains 5 different library preparation.
#### !!!!!!!!!!!!!!!!!!!!!!

STEP3_ALIGNMENT_JOB_ARR=$(qsub -q hotel -N mapping -l nodes=1:ppn=6,walltime=24:00:00 -t 1-24%5 \
                               -j oe -k oe -m oe -M ${email} \
                               -V -v ARG="$ARG" \
                               -W depend=afterokarray:${STEP2_DEMUX_JOB_ARR_id} \
                               ${code}/step3_alignment_array_jobs.sh)
echo "step3_alignment: $STEP3_ALIGNMENT_JOB_ARR"
STEP3_ALIGNMENT_JOB_ARR_id=$(echo "${STEP3_ALIGNMENT_JOB_ARR}" | cut -d '.' -f 1 )
#### !!!!!!!!!!!!!!!!!!!!!!
#### Number of array jobs needs modifications.
#### Here we have 24 array jobs
#### !!!!!!!!!!!!!!!!!!!!!!

STEP4_MKDUP_JOB_ARR=$(qsub -q hotel -N markDup -l nodes=1:ppn=12,walltime=24:00:00 -t 1-24%5 \
                           -j oe -k oe -m oe -M ${email} \
                           -V -v ARG="$ARG" \
                           -W depend=afterokarray:${STEP3_ALIGNMENT_JOB_ARR_id} \
                           ${code}/step4_mkDup_array_jobs.sh)
echo "step4_markDuplicates: $STEP4_MKDUP_JOB_ARR"
STEP4_MKDUP_JOB_ARR_id=$(echo "${STEP4_MKDUP_JOB_ARR}" | cut -d '.' -f 1 )
#### !!!!!!!!!!!!!!!!!!!!!!
#### Number of array jobs needs modifications.
#### Here we have 24 array jobs
#### !!!!!!!!!!!!!!!!!!!!!!

STEP5_STITCH_JOB_ARR=$(qsub -q hotel -N stitch -l nodes=1:ppn=24,walltime=100:00:00 -t 1-22%4 \
                            -j oe -k oe -m oe -M ${email} \
                            -V -v ARG="$ARG",PREV_BAMS="$PREV_BAMS" \
                            -W depend=afterokarray:${STEP4_MKDUP_JOB_ARR_id} \
                            ${code}/step5_stitch_variantCalling_array_jobs.sh)
echo "step5_stitch: $STEP5_STITCH_JOB_ARR"
STEP5_STITCH_JOB_ARR_id=$(echo "${STEP5_STITCH_JOB_ARR}" | cut -d '.' -f 1 )

STEP6_BEAGLE_JOB_ARR=$(qsub -q hotel -N beagle -l nodes=1:ppn=24,walltime=168:00:00 -t 1-21%4 \
                            -j oe -k oe -m oe -M ${email} \
                            -V -v ARG="$ARG" \
                            -W depend=afterokarray:${STEP5_STITCH_JOB_ARR_id} \
                            ${code}/step6_beagle_imputation_array_jobs.sh)
echo "step6_beagle: $STEP6_BEAGLE_JOB_ARR"
STEP6_BEAGLE_JOB_ARR_id=$(echo "${STEP6_BEAGLE_JOB_ARR}" | cut -d '.' -f 1 )

#### QC pipeline
QC1_MULTIQC_JOB=$(qsub -q condo -N qc -l nodes=1:ppn=6,walltime=8:00:00 -t 1-5 \
                       -j oe -k oe -m oe -M ${email} \
                       -V -v ARG="$ARG" \
                       -W depend=afterokarray:${STEP4_MKDUP_JOB_ARR_id} \
                       ${code}/QC1_multiqc_array_jobs.sh)
echo "QC1_multiQC: $QC1_MULTIQC_JOB"
#### !!!!!!!!!!!!!!!!!!!!!!
#### Number of array jobs needs modifications.
#### Here we have 5 array jobs since the flow cell 
#### contains 5 different library preparation.
#### !!!!!!!!!!!!!!!!!!!!!!

QC2_MAPPINGRESULT_JOB=$(qsub -q hotel -N mapping_stat -l nodes=1:ppn=12,walltime=168:00:00 \
                             -j oe -k oe -m oe -M ${email} \
                             -V -v ARG="$ARG" \
                             -W depend=afterokarray:${STEP4_MKDUP_JOB_ARR_id} \
                             ${code}/QC2_mappingResult.sh)
echo "QC2_mapping_results: $QC2_MAPPINGRESULT_JOB"

QC3_GENOTYPERESULT_JOB=$(qsub -q hotel -N genotype_stat -l nodes=1:ppn=12,walltime=168:00:00 \
                              -j oe -k oe -m oe -M ${email} \
                              -V -v ARG="$ARG",METADATA="$METADATA",PEDIGREE="$PEDIGREE" \
                              -W depend=afterokarray:${STEP6_BEAGLE_JOB_ARR_id} \
                              ${code}/QC3_genotypeResult.sh)
echo "QC3_genotype_results: $QC3_GENOTYPERESULT_JOB"