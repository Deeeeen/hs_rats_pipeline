#!/bin/bash

pipeline_arguments=pipeline_arguments
code=$(head -n 8 ${pipeline_arguments} | tail -n 1)

#### submit job array
STEP1_PREP_JOB=$(qsub ${code}/step1_prep.sh)
echo "step1_prep: $STEP1_PREP_JOB"

STEP2_DEMUX_JOB_ARR=$(qsub -W depend=afterok:$STEP1_PREP_JOB ${code}/step2_demux_array_jobs.sh)
echo "step2_demux: $STEP2_DEMUX_JOB_ARR"
STEP2_DEMUX_JOB_ARR_id=$(echo "${STEP2_DEMUX_JOB_ARR}" | cut -d '.' -f 1 )

STEP3_ALIGNMENT_JOB_ARR=$(qsub -W depend=afterokarray:${STEP2_DEMUX_JOB_ARR_id} ${code}/step3_alignment_array_jobs.sh)
echo "step3_alignment: $STEP3_ALIGNMENT_JOB_ARR"
STEP3_ALIGNMENT_JOB_ARR_id=$(echo "${STEP3_ALIGNMENT_JOB_ARR}" | cut -d '.' -f 1 )

STEP4_MKDUP_JOB_ARR=$(qsub -W depend=afterokarray:${STEP3_ALIGNMENT_JOB_ARR_id} ${code}/step4_mkDup_array_jobs.sh)
echo "step4_markDuplicates: $STEP4_MKDUP_JOB_ARR"
STEP4_MKDUP_JOB_ARR_id=$(echo "${STEP4_MKDUP_JOB_ARR}" | cut -d '.' -f 1 )

STEP5_STITCH_JOB_ARR=$(qsub -W depend=afterokarray:${STEP4_MKDUP_JOB_ARR_id} ${code}/step5_stitch_variantCalling_array_jobs.sh)
echo "step5_stitch: $STEP5_STITCH_JOB_ARR"
STEP5_STITCH_JOB_ARR_id=$(echo "${STEP5_STITCH_JOB_ARR}" | cut -d '.' -f 1 )

STEP6_BEAGLE_JOB_ARR=$(qsub -W depend=afterokarray:${STEP5_STITCH_JOB_ARR_id} ${code}/step6_beagle_imputation_array_jobs.sh)
echo "step6_beagle: $STEP6_BEAGLE_JOB_ARR"
STEP6_BEAGLE_JOB_ARR_id=$(echo "${STEP6_BEAGLE_JOB_ARR}" | cut -d '.' -f 1 )
