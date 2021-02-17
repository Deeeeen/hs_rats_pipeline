![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# TSCC Genotyping
## Source code for HS rats pipeline on [TSCC](https://www.sdsc.edu/support/user_guides/tscc.html)
:information_source: :information_source: :information_source:  **INFOMATION** :information_source: :information_source: :information_source:  
All the software installed on TSCC are [here](https://aapalmer-lab.slack.com/files/T0JULRU14/FPS2923NU). To change the software version used on this pipeline, you may need to do it manually. (Only the software for the pipeline main structure are listed. Software used on the result analysis part are NOT listed on this Github.)  

## Contents
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
