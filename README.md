# Final Project Outline

The makeup for this project starts on Feb 10. The expected time of completion would be 1 month from the start, ending on March 9. There will be a checkpoint or milestone discussion every Thursday during the period, through Google Calender invitation. The meeting time is subject to change upon negotiation in advance (either in person or online).

## Title
Differential Gene Expression in TCGA within Stage II gastric cancer between asian and white using DeSEQ2.

## Author
Linghao Zou (Loe)

## Contact
zoul@usc.edu

## Overview of project

I identified differentially expressed genes between Stomach Cancer Adenocarcomas and Adenocarcinomas between white vs. asians. This analysis will utilize the package DeSEQ2 and follow the specific [vignette and bioconductor package](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)(DESEQ). For this analysis, I'll use the TCGA cohort and have identified 205 star_counts open files for tumors that fit within my cohort with 144 whites and 46 asians. Within the analysis, I will control for gender, age and disease type. 

## Data

I will use the data from [GDC Portal](https://portal.gdc.cancer.gov/repository). By examining the clinical data, there are 190 tumor samples (but 205 files in total, as some samples have relapses and thus two associated files), among which 144 are whites and 46 are asians. The specific files are available under this repository ([click to view](https://github.com/Loe-zou/Final-Project/blob/main/clinical.tsv)).

*********

## Milestone 1

> Define the milestone based on initial view of the data.

**Check Point:** Feb 16

**Data Wrangling and Loading (Table and Sheets formed

**Due Date:** Feb 23

**Data fully loaded into vignette through DESeq2 steps (Star-Count section) ** Initial analysis through the vignette starts here.


## Milestone 2 

> Define the milestone 2 which includes initial running of the entire vignette for feedback.

**Due Date:** March 2

**An initial completion of the vignette.** I will complete an entire first draft of the data analysis. Data will be loaded into the vignette (through star_count), marked with errors and bugs, to be discussed with the supervisor. 

## Deliverable

**Due Date:** March 9

A complete repository with clear documentation and description of the analysis with results will be delivered by the due date.

## Method - Data Wrangling and Inputing (Pre-Process), all in bash/unix

Downloaded from gdc portal, following the data selection for comparation as mentioned above under "Data", the folder labeled "gdc_download_20221115_212727.200657" is found in the local.

For each sample file inside the folder, extract the 3rd column (unstranded counts) and convert it into a txt file, which is named for the first 4 characters of the original file.
```{bash}
awk'{print$4}' sample_genes.tsv > sample_genes.txt
```
For instance,
```{bash}
awk'{print$4}' e6fd3732-7e90-46ec-82fc-d1dcb8849e67.rna_seq.augmented_star_gene_counts.tsv > E6FD.txt
```

Replace the header of each sample gene file to match the file name (23 samples selected in total for analysis)
 **Credited to EunkSung**
```
./replace_header.sh
```

Extract the gene_id and remove the 1st column from one of the sample files, for instance:

```
awk '{print $1}'e6fd3732-7e90-46ec-82fc-d1dcb8849e67.rna_seq.augmented_star_gene_counts.tsv > ID.txt
tail -n +2 ID.txt > newID.txt
```
Remove the row 2-4 (if printed, remove them in the folder):
```
awk '!/^N_*/' newID.txt > newID.txt
```

Paste the new gene_id（newID.txt) and each sample counts (txt) together to obtain the final_genes for package import
```
paste gene_id.txt E6FD.txt F6FE.txt 0D85.txt 0FDB.txt 1BA6.txt 1D58.txt 2C5F.txt 5FFC.txt 49A0.txt 59EF.txt 77EE.txt 84CE.txt 69EE.txt 6272.txt 3468.txt 5497.txt 3641.txt B69A.txt B120.txt BD06.txt B502.txt D6BB.txt E1A1.txt > final_genes.txt
```

## Generate Race Table 

Follow the same logic, use awk to extract the column $16 from clinical.tsv > txt and $0 of each sample gene txt to paste them together to redirect into > race_table.csv and race_table.tsv (found in the repo)


## R Studio Data Importing

```{r}
install.packages('BiocManager')
BiocManager::install(c("DESeq2")
                     
library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ batch + condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

BiocManager::install("rnaseqGene")

library(BiocManager)
BiocManager::install("Named The Package You Need")

setwd('~/Desktop/510makeup')

# matrix matching
count_matrix <- read.delim("~/Desktop/510makeup/raw/sample_genes.txt", header=T, sep="\t")

# Column change and sample sheet
row.names(count_matrix) <- count_matrix$gene_id
count_matrix <- count_matrix[-c(1)]

sampletable <- read_tsv('~/Desktop/510makeup/raw/race_table.tsv')
```




