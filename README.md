# A Clinically Annotated Transcriptomic Atlas of Nervous System Tumors
While DNA methylation signatures are distinct across various nervous system neoplasms, it has not been comprehensively demonstrated whether transcriptomic signatures exhibit similar uniqueness. Additionally, there is a lack of a single, large-scale dataset for comparative gene expression analyses across these neoplasms. To address these gaps, we generated a large-scale, clinically annotated transcriptomic dataset comprising 5,402 neoplastic and 1,973 non-neoplastic nervous system samples from the public domain. We reprocessed, integrated, and reclassified samples with questionable diagnoses. Visualization using machine learning tools showed clustering primarily based on diagnosis, confirming that the bulk transcriptomic signature of nervous system neoplasms is unique across the diagnostic spectrum. The dataset’s broad coverage of diagnoses, including rarely studied entities, spans all ages and includes individuals from diverse ethnic backgrounds, enhancing its utility for comparative gene expression analyses. Our methods can be used to integrate and harmonize existing raw transcriptomic data of rare conditions, increasing their utility.

<p align="center">
     <img src="https://github.com/axitamm/BrainTumorAtlas/test.jpeg?raw=true" width="750px" align="center", class="only-light" >
</p>

## Instructions
### Gene Expression Data
The ExpData.rds files (zipped and segmented into 11 files) in the Gene Expression Folder contains the gene expression data for each sample. Once downloaded and unzipped, the file (and its data) can be loaded into R with the following simple command. 

ExpData<-read.RDS(“ExpData.rds”)

### Clinical Data
The clinical data associated with the samples are available in Atlas_Data.csv.

## Reproducibilty
All the code files needed to reproduce the results and figures in the original manuscript are available in the R code folder above.