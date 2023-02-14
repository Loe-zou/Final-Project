# Final Project Outline


##Updated Dec 15: Confirmed Covid and still experiencing symptoms, new deadline is Dec 16 with double extension.

## Title
Differential Gene Expression in TCGA within Stage II gastric cancer between asian and white using DeSEQ2.

## Author
Linghao Zou (Loe)

## Overview of project

I identified differentially expressed genes between Stomach Cancer Adenocarcomas and Adenocarcinomas for white and vs. asians. This analysis will utilize the package DeSEQ2 and follow the specific [vignette](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)(DESEQ). For this analysis, I'll use the TCGA cohort and have identified 205 star_counts open files for tumors that fit within my cohort with 144 whites and 46 asians. Within the analysis, I will control for gender, age and disease type.

[Vignette and bioconductor source link](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

## Data

I will use the data from [GDC Portal](https://portal.gdc.cancer.gov/repository). Examining clinical data, there are 190 tumor samples (but 205 files in total, as some samples have relapses and thus two associated files), and 144 are whites and 46 are asians. The specific files are available in the repository [here](https://github.com/Loe-zou/Final-Project/blob/main/clinical.tsv).

## Method - Data Wrangling and Input

```{bash}
awk'{print $4}' sample_ID.tsv > sample_ID.txt
```


*********

## Milestone 1

> Define the milestone based on initial view of the data.

**Check Point:** Feb 16
**Due Date:** Feb 23

** Data fully loaded into vignette through DESeq2 steps (Star-Count section) ** I will complete an entire first draft of analysis analyzed through the vignette.


## Milestone 2 

> Define the milestone 2 which includes initial running of the entire vignette for feedback.

**Due Date:** Tuesday March 1

**An initial completion of vignette.** I will complete an entire first draft of analysis analyzed through the vignette. Data loaded into vignette (through star_count), for seeking feedback.  Not all sections in the writing will be completed, but will be final project.


## Deliverable

**Due Date:** March 7

A complete repository with clear documentation and description of the analysis and results will be delivered by the due date.
