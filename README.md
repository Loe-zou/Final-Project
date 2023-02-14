# Final Project Outline

The makeup for this project starts on Feb 10. The expectaion of completion would be 1 month from the start, ending on March 9. There will be a checkpoint or milestone discussion every Thursday during the period, through Google Calender invitation. The meeting time is subject to change upon negotiation in advance.

## Title
Differential Gene Expression in TCGA within Stage II gastric cancer between asian and white using DeSEQ2.

## Author
Linghao Zou (Loe)

## Contact
zoul@usc.edu

## Overview of project

I identified differentially expressed genes between Stomach Cancer Adenocarcomas and Adenocarcinomas between white vs. asians. This analysis will utilize the package DeSEQ2 and follow the specific [vignette and bioconductor package](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)(DESEQ). For this analysis, I'll use the TCGA cohort and have identified 205 star_counts open files for tumors that fit within my cohort with 144 whites and 46 asians. Within the analysis, I will control for gender, age and disease type. 

## Data

I will use the data from [GDC Portal](https://portal.gdc.cancer.gov/repository). By examining the clinical data, there are 190 tumor samples (but 205 files in total, as some samples have relapses and thus two associated files), and 144 are whites and 46 are asians. The specific files are available under this repository ([click to view](https://github.com/Loe-zou/Final-Project/blob/main/clinical.tsv)).

*********

## Milestone 1

> Define the milestone based on initial view of the data.

**Check Point:** Feb 16
**Due Date:** Feb 23

** Data fully loaded into vignette through DESeq2 steps (Star-Count section) ** Initial analysis through the vignette starts here.


## Milestone 2 

> Define the milestone 2 which includes initial running of the entire vignette for feedback.

**Due Date:** March 2

**An initial completion of vignette.** I will complete an entire first draft of analysis analyzed through the vignette. Data loaded into vignette (through star_count), for seeking feedback.  Not all sections in the writing will be completed, but will be final project.

## Deliverable

**Due Date:** March 9

A complete repository with clear documentation and description of the analysis and results will be delivered by the due date.

## Method - Data Wrangling and Input

```{bash}
awk'{print $4}' sample_ID.tsv > sample_ID.txt
```



