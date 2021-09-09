# EMT_Scoring_RNASeq

This code takes raw read counts from the user and calculated TPM expression data and calculates EMT score using methods mentioned in Byers et al.,2013, Tan et al., 2014, George et al., 2017 and Chakraborty., et al 2020

### Preprocessing of the raw data

After QC filtering, align using STAR-aligner using appropriate reference genome (hg38/mm10) and then calculate raw read counts using htseq-count.

### Running the code
####  ***Make Data, Data_generated, Plot files and Output folders***

Using the above said raw counts, source the all_GSE_EMT.R code

This code will calculate the TPM expression data, EMT scores of three metrics: 76GS, KS and MLR. Also the correlations between these three metrics.


***If you have any problem running this code, please let me know***

