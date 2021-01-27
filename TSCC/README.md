![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# TSCC
## Source code for HS rats pipeline on [TSCC](https://www.sdsc.edu/support/user_guides/tscc.html)
:information_source: :information_source: :information_source:  **INFOMATION** :information_source: :information_source: :information_source:  
All the software installed on TSCC are [here](https://aapalmer-lab.slack.com/files/T0JULRU14/FPS2923NU). To change the software version used on this pipeline, you may need to do it manually. (Only the software for the pipeline main structure are listed. Software used on the result analysis part are NOT listed on this Github.)  

## Contents
**[HS_Rats_Genotyping_Pipeline_Summary_Report.pdf](HS_Rats_Genotyping_Pipeline_Summary_Report.pdf)**  
Graphic explanation for the pipeline structure and documentation.  

**[submission.sh](submission.sh)**  
Submission script for the pipeline.  

**[pipeline_arguments](pipeline_arguments)**  
Pipeline arguments.  
Line 1: home directory  
Line 2: Flow cell directory  
Line 3: Flow cell metadata
Line 4: Sequencing data directory  
Line 5: Reference genome  
Line 6: Reference panels for STITCH  
Line 7: Genetic map for BEAGLE  
Line 8: Directory where you keep the code for the pipeline  

**[previous_flow_cells_metadata](previous_flow_cells_metadata)**  
Paths to previous flow cells' metadata.  

**[previous_flow_cells_bams](previous_flow_cells_bams)**  
Paths to previous flow cells' BAM files.

**[pedigree_data](pedigree_data)**  
Paths to all flow cells' pedigree data. 

**[step1_prep.sh](step1_prep.sh)**  
Preparation for the pipeline, which constructs the basic structure of the directory and splits the sample sheet for Fgbio demultiplex based on each Riptide library preparation.

**[step2_demux_array_jobs.sh](step2_demux_array_jobs.sh)**  
Use array jobs feature on TSCC PBS to <ins>demultiplex the fastq files ([Fgbio 1.2.0](http://fulcrumgenomics.github.io/fgbio/))</ins> in parallel.  

**[step3_alignment_array_jobs.sh](step3_alignment_array_jobs.sh)**  
Use array jobs feature on TSCC PBS to <ins>map the sequencing reads to reference genome ([BWA 0.7.12](http://bio-bwa.sourceforge.net/index.shtml))</ins> in parallel. 

**[step4_mkDup_array_jobs.sh](step4_mkDup_array_jobs.sh)**  
Use array jobs feature on TSCC PBS to <ins>convert SAM files to BAM files ([Samtools 1.10](http://www.htslib.org/)), sort BAM files ([Samtools 1.10](http://www.htslib.org/)), mark PCR duplicates on BAM files ([Picard 2.23.3](https://broadinstitute.github.io/picard/)) and index the marked-duplicates BAM files ([Samtools 1.10](http://www.htslib.org/))</ins> in parallel.

**[step5_stitch_variantCalling_array_jobs.sh](step5_stitch_variantCalling_array_jobs.sh)**  
Use array jobs feature on TSCC PBS to <ins>do variant calling ([STITCH 1.6.3](https://github.com/rwdavies/STITCH))</ins> in parallel.  

**[step6_beagle_imputation_array_jobs.sh](step6_beagle_imputation_array_jobs.sh)**  
Use array jobs feature on TSCC PBS to <ins>do imputation ([BEAGLE 4.1](https://faculty.washington.edu/browning/beagle/b4_1.html))</ins> in parallel. 

**[QC1_multiqc_array_jobs.sh](QC1_multiqc_array_jobs.sh)**  
Run FastQC on fastq files and run Qualimap on marked-duplicates BAM files, and run MultiQC on the results of FastQC, Qualimap and Picard DuplicationMetrics. Resutls are separated by Riptide library preparation.

**[QC2_mappingResult.sh](QC2_mappingResult.sh)**  
Make plots for demultiplex results
1. Boxplot and scatter plot for number of matched reads for each library
2. Boxplot for % of unmatched reads  

Make plots for alignment results  
1. Boxplot and scatter plot for mapped read pairs
2. Boxplot and scatter plot for unmapped reads
3. Boxplot and scatter plot for duplications
4. Boxplot and scatter plot for uniquely mapped reads
5. Boxplot for mapped reads on each chromosome
6. Boxplot for GC content of mapped reads on each chromosome  

**[QC3_genotypeResult.sh](QC3_genotypeResult.sh)**  
Make plots for genotype results  
1. Histogram for STITCH info score
2. Heterozygosity rate vs. missing rate after STITCH
3. SNPs stats and density plot after BEAGLE
4. MAF histogram after BEAGLE
5. HWE histogram after BEAGLE
6. HWE vs. MAF heatmap after BEAGLE
7. PCA plots after BEAGLE
8. Albino coat color QC based on SNP 1:151097606
9. Pairwise concordance check

## Documentation
### Before running the pipeline on TSCC:
Please update the following files to suit your purpose:  
1. [pipeline_arguments](pipeline_arguments)
2. [previous_flow_cells_metadata](previous_flow_cells_metadata)
3. [previous_flow_cells_bams](previous_flow_cells_bams)
4. [pedigree_data](pedigree_data)
5. Update the PBS Torque arguments and the corresponding file locations on [submission.sh](submission.sh).


### Run the pipeline on TSCC:
1. Change the permission of the submission script
```
chmod u+x submission.sh
```
2. Run the submission script
```
./submission.sh
```  
**TO-DO**  
- [ ] Output log file
- [ ] Brown coat color QC