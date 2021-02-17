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
