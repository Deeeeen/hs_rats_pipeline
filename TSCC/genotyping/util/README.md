![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# TSCC Genotyping Utilities
## Source code for HS rats pipeline on [TSCC](https://www.sdsc.edu/support/user_guides/tscc.html)
:information_source: :information_source: :information_source:  **INFOMATION** :information_source: :information_source: :information_source:  
All the software installed on TSCC are [here](https://aapalmer-lab.slack.com/files/T0JULRU14/FPS2923NU). To change the software version used on this pipeline, you may need to do it manually. (Only the software for the pipeline main structure are listed. Software used on the result analysis part are NOT listed on this Github.)  

## Contents
**[separate_metadata.py](separate_metadata.py)**  
Separated and organize sample metadata based on the needs for Fgbio.  

**[variant_calling_STITCH.r](variant_calling_STITCH.r)**  
Rscript to call STITCH.  