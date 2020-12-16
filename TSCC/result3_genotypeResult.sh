#!/bin/bash
#PBS -q hotel
#PBS -N genotype_stat
#PBS -l nodes=1:ppn=4
#PBS -l walltime=168:00:00
#PBS -V
#PBS -j oe
#PBS -k oe
#PBS -M dec037@health.ucsd.edu
#PBS -m ae

#### Make sure to change the compute node, job name, num of node,
#### ppn, walltime, forward email address for notifications. (above)

pipeline_arguments=pipeline_arguments
previous_flow_cells_metadata=previous_flow_cells_metadata
pedigree_data=pedigree_data

home=$(head -n 1 ${pipeline_arguments})
code=$(head -n 8 ${pipeline_arguments} | tail -n 1)
dir_path=$(head -n 2 ${pipeline_arguments} | tail -n 1)
original_sample_sheet=$(head -n 3 ${pipeline_arguments} | tail -n 1)
stitch_path=${dir_path}/stitch
beagle_path=${dir_path}/beagle
out_path=${dir_path}/results
genotype_result=${out_path}/genotype_result
vcf_prefix=$(head -n 9 ${pipeline_arguments} | tail -n 1)
current_date=$(date +%m%d%Y)
#### read in arguments for the pipeline

cd ${home}
mkdir ${genotype_result}

####!!!!!!!!!!!!!!!!!!
####TO-DO: code to output log file
####       code to plot PCA, code to plot het vs missing,
####       code to plot mssing vs number of reads, code for coat color QC
####       code for pairwise concordance 
####!!!!!!!!!!!!!!!!!!
################ Output log file, and combine all metadata ####################
source activate hs_rats
Rscript ${code}/combine_csv.r \
  ${vcf_prefix}_${current_date} \
  ${original_sample_sheet} \
  ${previous_flow_cells_metadata} \
  ${pedigree_data} \
  ${out_path}
conda deactivate

##################### genotypes results after STITCH ##########################
mkdir ${genotype_result}/stitch_result
mkdir ${genotype_result}/stitch_result/plink
########### index stitch vcf files
fs_in=$(ls ${stitch_path}/*.vcf.gz)
ncpu=3
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
for f in $fs_in
do
    while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   /projects/ps-palmer/software/local/src/bcftools-1.10.2/bcftools index -f -t ${f}
done
while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
      sleep 60
done

########### concat stitch .vcf.gz
START=$(date +%s)
fs_in=$(ls ${stitch_path}/*.vcf.gz)
/projects/ps-palmer/software/local/src/bcftools-1.10.2/bcftools concat \
  --no-version -a -d none -O z -o ${stitch_path}/${vcf_prefix}_stitch_${current_date}.vcf.gz ${fs_in}
END=$(date +%s)
echo "Concat stitch vcfs Time elapsed: $(( $END - $START )) seconds"

#### imputation INFO scores
START=$(date +%s)
/projects/ps-palmer/software/local/src/bcftools-1.10.2/bcftools index -f -t ${stitch_path}/${vcf_prefix}_stitch_${current_date}.vcf.gz
/projects/ps-palmer/software/local/src/bcftools-1.10.2/bcftools query -f '%POS\t%INFO/INFO_SCORE\n' \
  ${stitch_path}/${vcf_prefix}_stitch_${current_date}.vcf.gz > ${genotype_result}/stitch_result/stitch_INFO

source activate hs_rats
Rscript ${code}/stitch_info_score.r \
   ${genotype_result}/stitch_result/stitch_INFO \
   ${genotype_result}/stitch_result
conda deactivate

END=$(date +%s)
echo "Imputation INFO score Time elapsed: $(( $END - $START )) seconds"

#### STITCH VCF to plink binary file
START=$(date +%s)
/projects/ps-palmer/software/local/src/plink2 --vcf ${stitch_path}/${vcf_prefix}_stitch_${current_date}.vcf.gz \
  --set-missing-var-ids @:# --make-bed --out ${genotype_result}/stitch_result/plink/${vcf_prefix}_stitch_${current_date}
END=$(date +%s)
echo "Stitch VCF -> plink Time elapsed: $(( $END - $START )) seconds"

#### missing rate vs number of reads
START=$(date +%s)
/projects/ps-palmer/software/local/src/plink-1.90/plink --bfile ${genotype_result}/stitch_result/plink/${vcf_prefix}_stitch_${current_date} \
  --missing --out ${genotype_result}/stitch_result/plink/${vcf_prefix}_stitch_${current_date}
END=$(date +%s)
echo "Missing rate calculation Time elapsed: $(( $END - $START )) seconds"

####!!!!!!!!!!!!!!!!!!
####TO-DO: code to plot the missing rate vs. number of reads
####!!!!!!!!!!!!!!!!!!

#### heterozygosity vs missing rate
START=$(date +%s)
/projects/ps-palmer/software/local/src/plink-1.90/plink --bfile ${genotype_result}/stitch_result/plink/${vcf_prefix}_stitch_${current_date} \
  --het --out ${genotype_result}/stitch_result/plink/${vcf_prefix}_stitch_${current_date}
END=$(date +%s)
echo "Heterozygosity rate calculation Time elapsed: $(( $END - $START )) seconds"

####!!!!!!!!!!!!!!!!!!
####TO-DO: code to plot the heterozygosity rate vs. missing rate
####!!!!!!!!!!!!!!!!!!

##################### genotypes results after BEAGLE ##########################
mkdir ${genotype_result}/beagle_result
mkdir ${genotype_result}/beagle_result/plink
########### index beagle vcf files
fs_in=$(ls ${beagle_path}/*.vcf.gz)
ncpu=3
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
for f in $fs_in
do
    while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   /projects/ps-palmer/software/local/src/bcftools-1.10.2/bcftools index -f -t ${f}
done
while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
      sleep 60
done

########### BEAGLE stats
START=$(date +%s)
ncpu=3
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
fs_in=$(ls ${beagle_path}/*.vcf.gz)
cnt=0
for f in $fs_in
do
    while [ "$(jobs -rp | wc -l)" -ge $ncpu ]; do
      sleep 60
   done
   sleep 5
   (( cnt += 1 ))

   temp=$(echo ${f} | rev | cut -d '.' -f 3 | cut -d '/' -f 1 | rev)
   /projects/ps-palmer/software/local/src/bcftools-1.10.2/bcftools stats -v ${f} > ${genotype_result}/beagle_result/${temp}_stats
done

source activate hs_rats
python3 ${code}/variant_stats.py \
    ${genotype_result}/beagle_result \
    ${genotype_result}/beagle_result/beagle
conda deactivate

END=$(date +%s)
echo "BEAGLE variant calling stats Time elapsed: $(( $END - $START )) seconds"

########### concat beagle .vcf.gz
START=$(date +%s)

fs_in=$(ls ${beagle_path}/*.vcf.gz)
/projects/ps-palmer/software/local/src/bcftools-1.10.2/bcftools concat \
  --no-version -a -d none -O z -o ${beagle_path}/${vcf_prefix}_beagle_${current_date}.vcf.gz ${fs_in}
END=$(date +%s)
echo "Concat beagle vcfs Time elapsed: $(( $END - $START )) seconds"

#### SNP density
START=$(date +%s)
module load bcftools
bgzip -d -c ${beagle_path}/${vcf_prefix}_beagle_${current_date}.vcf.gz > ${beagle_path}/${vcf_prefix}_beagle_${current_date}.vcf
module unload bcftools

source activate hs_rats
python3 ${code}/snp_density.py -b 1000 \
    -i ${beagle_path}/${vcf_prefix}_beagle_${current_date}.vcf \
    -o ${genotype_result}/beagle_result/beagle_snp_density
conda deactivate
END=$(date +%s)
echo "SNP density Time elapsed: $(( $END - $START )) seconds"

#### convert vcf to bed, HWE, FRQ, PCA
START=$(date +%s)
ncpu=3
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
fs_in=$(ls ${beagle_path}/*.vcf.gz)
/projects/ps-palmer/software/local/src/plink2 --vcf ${beagle_path}/${vcf_prefix}_beagle_${current_date}.vcf.gz \
  --set-missing-var-ids @:# --make-bed --out ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle_${current_date}
/projects/ps-palmer/software/local/src/plink2 --bfile ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle_${current_date} \
  --hardy --out ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle_${current_date}
/projects/ps-palmer/software/local/src/plink-1.90/plink --bfile ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle_${current_date} \
  --freq --out ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle_${current_date}
/projects/ps-palmer/software/local/src/plink2 --bfile ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle_${current_date} \
  --pca --out ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle_${current_date}

source activate hs_rats
python3 ${code}/hwe_maf.py \
    ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle_${current_date} \
    ${genotype_result}/beagle_result/beagle
conda deactivate
END=$(date +%s)
echo "VCF -> BED, HWE, FRQ PCA variant calling stats Time elapsed: $(( $END - $START )) seconds"

