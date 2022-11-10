# Final Project Outline

## Title
Differential Gene Expression in TCGA within Stage II gastric cancer between asian and white using DeSEQ2.

## Author
Linghao Zou (Loe)

## Overview of project

I will identify differentially expressed genes between Stomach Cancer Adenocarcomas and Adenocarcinomas for white and vs. asians. This analysis will utilize the package DeSEQ2 and follow the specific vignette: [http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html](DESEQ). For this analysis, I'll use the TCGA cohort and have identified 205 star_counts open files for tumors that fit within my cohort with 144 whites and 46 asians. Within the analysis, I will control for gender, age and disease type.

Vignette: [http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

## Data

I will use the data from [https://portal.gdc.cancer.gov/repository](https://portal.gdc.cancer.gov/repository). Examining clinical data, there are 190 tumor samples (but 205 files in total, as some samples have relapses and thus two associated files), and 144 are whites and 46 are asians. The specific files are available in the repository: https://github.com/Loe-zou/Final-Project/blob/main/clinical.tsv.


## Milestone 1

> Define the milestone based on initial view of the data.

**Due Date:** Tuesday November 22nd

** Data fully loaded into vignette through DESeq2 steps (Star-Count section) ** I will complete an entire first draft of analysis analyzed through the vignette.


## Milestone 2

> Define the milestone 2 which includes initial running of the entire vignette for feedback.

**Due Date:** Tuesday November 29th

**An initial completion of vignette.** I will complete an entire first draft of analysis analyzed through the vignette. Data loaded into vignette (through star_count), for seeking feedback.  Not all sections in the writing will be completed, but will be final project.


## Deliverable

**Due Date:** December 3rd 11:50pm

A complete repository with clear documentation and description of the analysis and results will be delivered by the due date.
