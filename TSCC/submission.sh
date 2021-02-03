#!/bin/bash

ARG=/projects/ps-palmer/hs_rats/201218_A00953_0203_BHN5HYDSXY/code/pipeline_arguments
METADATA=/projects/ps-palmer/hs_rats/201218_A00953_0203_BHN5HYDSXY/code/previous_flow_cells_metadata
PEDIGREE=/projects/ps-palmer/hs_rats/201218_A00953_0203_BHN5HYDSXY/code/pedigree_data
PREV_BAMS=/projects/ps-palmer/hs_rats/201218_A00953_0203_BHN5HYDSXY/code/previous_flow_cells_bams
current_metadata=$(head -n 3 ${ARG} | tail -n 1)
code=$(head -n 8 ${ARG} | tail -n 1)

#### Make sure to change the compute node, job name, num of node,
#### ppn, walltime, forward email address for notifications. (below)
email=dec037@health.ucsd.edu

######################## Genotyping pipeline ########################
STEP1_PREP_JOB=$(qsub -q condo -N preDemux -l nodes=1:ppn=1,walltime=8:00:00 \
                      -j oe -k oe -m ae -M ${email} \
                      -V -v ARG="${ARG}" \
                      ${code}/step1_prep.sh)
echo "step1_prep: ${STEP1_PREP_JOB}"

#### a customized command to extract the number of library prep from metadata
source activate hs_rats
num_lib_py=$(cat <<'EOF'
import pandas as pd
import sys
origial_metadata = pd.read_csv(sys.argv[1])
metadata_cols = origial_metadata.columns.tolist()
origial_metadata = origial_metadata[origial_metadata['strain'] == 'Heterogenous stock'].reset_index(drop=True)
Library_ID = ""
for col in metadata_cols:
    if 'library_name' == col.lower():
        Library_ID = col
        break
sys.stdout.write(str(len(set(origial_metadata[Library_ID]))))
EOF
)
num_lib() { python -c "${num_lib_py}" "$@"; }

#### submit demux array jobs based on the number of library prep
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=4
num_jobs=$(num_lib ${current_metadata})
STEP2_DEMUX_JOB_ARR=$(qsub -q hotel -N demux -l nodes=1:ppn=${ppn},walltime=24:00:00 -t 1-${num_jobs} \
                           -j oe -k oe -m ae -M ${email} \
                           -V -v ARG="${ARG}",ppn="${ppn}" \
                           -W depend=afterok:${STEP1_PREP_JOB} \
                           ${code}/step2_demux_array_jobs.sh)
echo "step2_demux: ${STEP2_DEMUX_JOB_ARR}"
STEP2_DEMUX_JOB_ARR_id=$(echo "${STEP2_DEMUX_JOB_ARR}" | cut -d '.' -f 1 )

#### submit mapping array jobs
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### (optional) change num_jobs for the number of array jobs based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=6
num_jobs=24
STEP3_ALIGNMENT_JOB_ARR=$(qsub -q hotel -N mapping -l nodes=1:ppn=${ppn},walltime=24:00:00 -t 1-${num_jobs}%5 \
                               -j oe -k oe -m ae -M ${email} \
                               -V -v ARG="${ARG}",num_jobs="${num_jobs}",ppn="${ppn}" \
                               -W depend=afterokarray:${STEP2_DEMUX_JOB_ARR_id} \
                               ${code}/step3_alignment_array_jobs.sh)
echo "step3_alignment: ${STEP3_ALIGNMENT_JOB_ARR}"
STEP3_ALIGNMENT_JOB_ARR_id=$(echo "${STEP3_ALIGNMENT_JOB_ARR}" | cut -d '.' -f 1 )

#### submit marking duplicates array jobs
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### (optional) change num_jobs for the number of array jobs based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=12
num_jobs=24
STEP4_MKDUP_JOB_ARR=$(qsub -q hotel -N markDup -l nodes=1:ppn=${ppn},walltime=24:00:00 -t 1-${num_jobs}%5 \
                           -j oe -k oe -m ae -M ${email} \
                           -V -v ARG="${ARG}",num_jobs="${num_jobs}",ppn="${ppn}" \
                           -W depend=afterokarray:${STEP3_ALIGNMENT_JOB_ARR_id} \
                           ${code}/step4_mkDup_array_jobs.sh)
echo "step4_markDuplicates: ${STEP4_MKDUP_JOB_ARR}"
STEP4_MKDUP_JOB_ARR_id=$(echo "${STEP4_MKDUP_JOB_ARR}" | cut -d '.' -f 1 )

#### submit STITCH variant calling array jobs
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=24
STEP5_STITCH_JOB_ARR=$(qsub -q hotel -N stitch -l nodes=1:ppn=${ppn},walltime=100:00:00 -t 1-22%4 \
                            -j oe -k oe -m ae -M ${email} \
                            -V -v ARG="${ARG}",PREV_BAMS="${PREV_BAMS}" \
                            -W depend=afterokarray:${STEP4_MKDUP_JOB_ARR_id} \
                            ${code}/step5_stitch_variantCalling_array_jobs.sh)
echo "step5_stitch: ${STEP5_STITCH_JOB_ARR}"
STEP5_STITCH_JOB_ARR_id=$(echo "${STEP5_STITCH_JOB_ARR}" | cut -d '.' -f 1 )

#### submit BEAGLE imputation array jobs
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=24
ne=150
STEP6_BEAGLE_JOB_ARR=$(qsub -q hotel -N beagle -l nodes=1:ppn=${ppn},walltime=168:00:00 -t 1-21%4 \
                            -j oe -k oe -m ae -M ${email} \
                            -V -v ARG="${ARG}",ppn="${ppn}",ne="${ne}" \
                            -W depend=afterokarray:${STEP5_STITCH_JOB_ARR_id} \
                            ${code}/step6_beagle_imputation_array_jobs.sh)
echo "step6_beagle: ${STEP6_BEAGLE_JOB_ARR}"
STEP6_BEAGLE_JOB_ARR_id=$(echo "${STEP6_BEAGLE_JOB_ARR}" | cut -d '.' -f 1 )


######################## QC pipeline ########################
#### submit multiQC array jobs based on the number of library prep
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=6
num_jobs=$(num_lib ${current_metadata})
QC1_MULTIQC_JOB=$(qsub -q home -N qc -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
                       -j oe -k oe -m ae -M ${email} \
                       -V -v ARG="${ARG}",ppn="${ppn}" \
                       -W depend=afterokarray:${STEP4_MKDUP_JOB_ARR_id} \
                       ${code}/QC1_multiqc_array_jobs.sh)
echo "QC1_multiQC: ${QC1_MULTIQC_JOB}"

#### submit mapping stats array jobs based on the number of library prep
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=12
QC2_MAPPINGRESULT_JOB=$(qsub -q hotel -N mapping_stat -l nodes=1:ppn=${ppn},walltime=168:00:00 \
                             -j oe -k oe -m ae -M ${email} \
                             -V -v ARG="${ARG}",ppn="${ppn}" \
                             -W depend=afterokarray:${STEP4_MKDUP_JOB_ARR_id} \
                             ${code}/QC2_mappingResult.sh)
echo "QC2_mapping_results: ${QC2_MAPPINGRESULT_JOB}"

#### submit genotyping stats array jobs based on the number of library prep
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=12
QC3_GENOTYPERESULT_JOB=$(qsub -q hotel -N genotype_stat -l nodes=1:ppn=${ppn},walltime=168:00:00 \
                              -j oe -k oe -m ae -M ${email} \
                              -V -v ARG="${ARG}",METADATA="${METADATA}",PEDIGREE="${PEDIGREE}",ppn="${ppn}" \
                              -W depend=afterokarray:${STEP6_BEAGLE_JOB_ARR_id} \
                              ${code}/QC3_genotypeResult.sh)
echo "QC3_genotype_results: ${QC3_GENOTYPERESULT_JOB}"