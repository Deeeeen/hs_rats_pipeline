#!/bin/bash

pipeline_arguments=${ARG}
previous_flow_cells_metadata=${METADATA}
pedigree_data=${PEDIGREE}

home=$(head -n 1 ${pipeline_arguments})
code=$(head -n 8 ${pipeline_arguments} | tail -n 1)
dir_path=$(head -n 2 ${pipeline_arguments} | tail -n 1)
sample_sheet=${dir_path}/demux/sample_sheet.csv
stitch_path=${dir_path}/stitch
beagle_path=${dir_path}/beagle
out_path=${dir_path}/results
genotype_result=${out_path}/genotype_result
#### read in arguments for the pipeline

cd ${home}
mkdir ${genotype_result}

####!!!!!!!!!!!!!!!!!!
####TO-DO: code to output log file
####       code for brown coat color QC
####!!!!!!!!!!!!!!!!!!
########################## Combine all metadata ##############################
source activate hs_rats
Rscript ${code}/combine_csv.r \
  ${sample_sheet} \
  ${previous_flow_cells_metadata} \
  ${pedigree_data} \
  ${out_path}
conda deactivate

vcf_prefix=$(ls ${out_path}/hs_rats_n*_metadata.csv | rev | cut -d'_' -f 2- |cut -d'/' -f 1 | rev)

##################### genotypes results after STITCH ##########################
mkdir ${genotype_result}/stitch_result
mkdir ${genotype_result}/stitch_result/plink
########### index stitch vcf files
fs_in=$(ls ${stitch_path}/*.vcf.gz)
ncpu=${ppn}
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
  --no-version -a -d none -O z -o ${stitch_path}/${vcf_prefix}_stitch.vcf.gz ${fs_in}
END=$(date +%s)
echo "Concat stitch vcfs Time elapsed: $(( $END - $START )) seconds"

#### imputation INFO scores
START=$(date +%s)
/projects/ps-palmer/software/local/src/bcftools-1.10.2/bcftools index -f -t ${stitch_path}/${vcf_prefix}_stitch.vcf.gz
/projects/ps-palmer/software/local/src/bcftools-1.10.2/bcftools query -f '%POS\t%INFO/INFO_SCORE\n' \
  ${stitch_path}/${vcf_prefix}_stitch.vcf.gz > ${genotype_result}/stitch_result/stitch_INFO

source activate hs_rats
Rscript ${code}/stitch_info_score.r \
   ${genotype_result}/stitch_result/stitch_INFO \
   ${genotype_result}/stitch_result
conda deactivate

END=$(date +%s)
echo "Imputation INFO score Time elapsed: $(( $END - $START )) seconds"

#### STITCH VCF to plink binary file
START=$(date +%s)
/projects/ps-palmer/software/local/src/plink2 --vcf ${stitch_path}/${vcf_prefix}_stitch.vcf.gz \
  --set-missing-var-ids @:# --make-bed --out ${genotype_result}/stitch_result/plink/${vcf_prefix}_stitch
END=$(date +%s)
echo "Stitch VCF -> plink Time elapsed: $(( $END - $START )) seconds"

#### missing rate vs number of reads
START=$(date +%s)
/projects/ps-palmer/software/local/src/plink2 --bfile ${genotype_result}/stitch_result/plink/${vcf_prefix}_stitch \
  --missing --out ${genotype_result}/stitch_result/plink/${vcf_prefix}_stitch
END=$(date +%s)

temp_dir=$(echo ${dir_path} | rev | cut -d '/' -f 2- | rev)
source activate hs_rats
Rscript ${code}/missing_vs_reads.r \
   ${temp_dir} \
   ${dir_path} 
conda deactivate
echo "Missing rate calculation Time elapsed: $(( $END - $START )) seconds"

#### heterozygosity vs missing rate
START=$(date +%s)
/projects/ps-palmer/software/local/src/plink-1.90/plink --bfile ${genotype_result}/stitch_result/plink/${vcf_prefix}_stitch \
  --het --out ${genotype_result}/stitch_result/plink/${vcf_prefix}_stitch

source activate hs_rats
Rscript ${code}/het_vs_missing.r \
  ${genotype_result}/stitch_result/plink/${vcf_prefix}_stitch.smiss \
  ${genotype_result}/stitch_result/plink/${vcf_prefix}_stitch.het \
  ${genotype_result}/stitch_result/plink/ \
  ${vcf_prefix}
conda deactivate

END=$(date +%s)
echo "Heterozygosity rate calculation Time elapsed: $(( $END - $START )) seconds"

##################### genotypes results after BEAGLE ##########################
mkdir ${genotype_result}/beagle_result
mkdir ${genotype_result}/beagle_result/plink
########### index beagle vcf files
fs_in=$(ls ${beagle_path}/*.vcf.gz)
ncpu=12
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
for f in ${fs_in}
do
    while [ "$(jobs -rp | wc -l)" -ge ${ncpu} ]; do
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
ncpu=12
#### !!!!!!!!!!!!!!!!!!!!!!
#### May need to change the ncpu based on the ppn requested
#### !!!!!!!!!!!!!!!!!!!!!!
fs_in=$(ls ${beagle_path}/*.vcf.gz)
cnt=0
for f in ${fs_in}
do
    while [ "$(jobs -rp | wc -l)" -ge ${ncpu} ]; do
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
  --no-version -a -d none -O z -o ${beagle_path}/${vcf_prefix}_beagle.vcf.gz ${fs_in}
END=$(date +%s)
echo "Concat beagle vcfs Time elapsed: $(( $END - $START )) seconds"

#### SNP density
START=$(date +%s)
module load bcftools
bgzip -d -c ${beagle_path}/${vcf_prefix}_beagle.vcf.gz > ${beagle_path}/${vcf_prefix}_beagle.vcf
module unload bcftools

source activate hs_rats
python3 ${code}/snp_density.py -b 1000 \
    -i ${beagle_path}/${vcf_prefix}_beagle.vcf \
    -o ${genotype_result}/beagle_result/beagle_snp_density
conda deactivate
END=$(date +%s)
echo "SNP density Time elapsed: $(( $END - $START )) seconds"

#### convert vcf to bed, HWE, FRQ, PCA
START=$(date +%s)
fs_in=$(ls ${beagle_path}/*.vcf.gz)
/projects/ps-palmer/software/local/src/plink2 --vcf ${beagle_path}/${vcf_prefix}_beagle.vcf.gz \
  --set-missing-var-ids @:# --make-bed --out ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle
/projects/ps-palmer/software/local/src/plink2 --bfile ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle \
  --hardy --out ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle
/projects/ps-palmer/software/local/src/plink-1.90/plink --bfile ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle \
  --freq --out ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle
/projects/ps-palmer/software/local/src/plink2 --bfile ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle \
  --pca --out ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle

source activate hs_rats
python3 ${code}/hwe_maf.py \
  ${genotype_result}/beagle_result/plink \
  ${genotype_result}/beagle_result/beagle

Rscript ${code}/pca.r \
  ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle.eigenvec \
  ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle.eigenval \
  ${out_path}/${vcf_prefix}_metadata.csv \
  ${genotype_result}/beagle_result
conda deactivate
END=$(date +%s)
echo "VCF -> BED, HWE, FRQ PCA variant calling stats Time elapsed: $(( $END - $START )) seconds"

#### Albino coat color QC based on SNP 1:151097606
START=$(date +%s)
/projects/ps-palmer/software/local/src/plink-1.90/plink --bfile ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle \
  --alleleACGT --snp 1:151097606 --recode  --out ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle_albino_1_151097606

source activate hs_rats
Rscript ${code}/coat_color_albino.r \
  ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle_albino_1_151097606.ped \
  ${out_path}/${vcf_prefix}_metadata.csv \
  ${genotype_result}/beagle_result/plink
conda deactivate
END=$(date +%s)
echo "Albino coat color QC Time elapsed: $(( $END - $START )) seconds"

#### Brown coat color QC based on chr3:150285633​ chr3:150288295​ chr3:150432118​ 
#### chr3:150449245​ chr3:150488934​ chr3:150530733​ chr3:150574499​ chr3:150584522
START=$(date +%s)
/projects/ps-palmer/software/local/src/plink-1.90/plink --bfile ${genotype_result}/beagle_result/plink/hs_rats_n1912_01082021_beagle \
  --alleleACGT --snps 3:150285633, 3:150288295, 3:150432118, 3:150449245, 3:150488934, 3:150530733, 3:150584522\
  --recode  --out ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle_brown

source activate hs_rats
Rscript ${code}/coat_color_brown.r \
  ${genotype_result}/beagle_result/plink/${vcf_prefix}_beagle_brown.ped \
  ${out_path}/${vcf_prefix}_metadata.csv \
  ${genotype_result}/beagle_result/plink
conda deactivate
END=$(date +%s)
echo "Brown coat color QC Time elapsed: $(( $END - $START )) seconds"

#### Pairwise concordance check
START=$(date +%s)
/projects/ps-palmer/software/local/src/bcftools-1.10.2/bcftools gtcheck ${beagle_path}/${vcf_prefix}_beagle.vcf.gz > ${genotype_result}/beagle_result/genotypes_bcftools_gtcheck
grep '^ERR' ${genotype_result}/beagle_result/genotypes_bcftools_gtcheck > ${genotype_result}/beagle_result/genotypes_bcftools_gtcheck_ERR
grep '^CN' ${genotype_result}/beagle_result/genotypes_bcftools_gtcheck > ${genotype_result}/beagle_result/genotypes_bcftools_gtcheck_CN
grep '^CLUSTER' ${genotype_result}/beagle_result/genotypes_bcftools_gtcheck > ${genotype_result}/beagle_result/genotypes_bcftools_gtcheck_CLUSTER
grep '^TH' ${genotype_result}/beagle_result/genotypes_bcftools_gtcheck > ${genotype_result}/beagle_result/genotypes_bcftools_gtcheck_TH
grep '^DOT' ${genotype_result}/beagle_result/genotypes_bcftools_gtcheck > ${genotype_result}/beagle_result/genotypes_bcftools_gtcheck_DOT
rm ${genotype_result}/beagle_result/genotypes_bcftools_gtcheck

source activate hs_rats
Rscript ${code}/pairwise_concordance.r \
  ${genotype_result}/beagle_result/genotypes_bcftools_gtcheck_ERR \
  ${out_path}/${vcf_prefix}_metadata.csv \
  ${genotype_result}/beagle_result
conda deactivate
END=$(date +%s)
echo "Pairwise concordance check Time elapsed: $(( $END - $START )) seconds"

