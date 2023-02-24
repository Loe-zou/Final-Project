510Makeup
================
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

Here is the instruction from the bioconductor to set up the filtering:

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

## Part 1

``` r
BiocManager::install("DeSeq2")
```

    ## Bioconductor version 3.16 (BiocManager 1.30.19), R 4.2.1 (2022-06-23)

    ## Installing package(s) 'DeSeq2'

    ## Warning: package 'DeSeq2' is not available for Bioconductor version '3.16'
    ## 
    ## A version of this package for your version of R might be available elsewhere,
    ## see the ideas at
    ## https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages

    ## Warning: Perhaps you meant 'DESeq2' ?

    ## Old packages: 'BiocManager', 'RcppArmadillo'

``` r
BiocManager::install("dplyr")
```

    ## Bioconductor version 3.16 (BiocManager 1.30.19), R 4.2.1 (2022-06-23)

    ## Warning: package(s) not installed when version(s) same as or greater than current; use
    ##   `force = TRUE` to re-install: 'dplyr'

    ## Old packages: 'BiocManager', 'RcppArmadillo'

``` r
BiocManager::install("apeglm")
```

    ## Bioconductor version 3.16 (BiocManager 1.30.19), R 4.2.1 (2022-06-23)

    ## Warning: package(s) not installed when version(s) same as or greater than current; use
    ##   `force = TRUE` to re-install: 'apeglm'

    ## Old packages: 'BiocManager', 'RcppArmadillo'

``` r
BiocManager::install("vsn")
```

    ## Bioconductor version 3.16 (BiocManager 1.30.19), R 4.2.1 (2022-06-23)

    ## Warning: package(s) not installed when version(s) same as or greater than current; use
    ##   `force = TRUE` to re-install: 'vsn'

    ## Old packages: 'BiocManager', 'RcppArmadillo'

``` r
BiocManager::install("pheatmap")
```

    ## Bioconductor version 3.16 (BiocManager 1.30.19), R 4.2.1 (2022-06-23)

    ## Warning: package(s) not installed when version(s) same as or greater than current; use
    ##   `force = TRUE` to re-install: 'pheatmap'

    ## Old packages: 'BiocManager', 'RcppArmadillo'

``` r
library(DESeq2)
```

    ## Warning: package 'DESeq2' was built under R version 4.2.2

    ## Loading required package: S4Vectors

    ## Warning: package 'S4Vectors' was built under R version 4.2.2

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomicRanges

    ## Warning: package 'GenomicRanges' was built under R version 4.2.2

    ## Loading required package: GenomeInfoDb

    ## Warning: package 'GenomeInfoDb' was built under R version 4.2.2

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following object is masked from 'package:matrixStats':
    ## 
    ##     count

    ## The following objects are masked from 'package:GenomicRanges':
    ## 
    ##     intersect, setdiff, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tibble)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ readr     2.1.4
    ## ✔ ggplot2   3.4.1     ✔ stringr   1.5.0
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ lubridate::%within%() masks IRanges::%within%()
    ## ✖ dplyr::collapse()     masks IRanges::collapse()
    ## ✖ dplyr::combine()      masks Biobase::combine(), BiocGenerics::combine()
    ## ✖ dplyr::count()        masks matrixStats::count()
    ## ✖ dplyr::desc()         masks IRanges::desc()
    ## ✖ tidyr::expand()       masks S4Vectors::expand()
    ## ✖ dplyr::filter()       masks stats::filter()
    ## ✖ dplyr::first()        masks S4Vectors::first()
    ## ✖ dplyr::lag()          masks stats::lag()
    ## ✖ ggplot2::Position()   masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()       masks GenomicRanges::reduce(), IRanges::reduce()
    ## ✖ dplyr::rename()       masks S4Vectors::rename()
    ## ✖ lubridate::second()   masks S4Vectors::second()
    ## ✖ lubridate::second<-() masks S4Vectors::second<-()
    ## ✖ dplyr::slice()        masks IRanges::slice()
    ## ℹ Use the ]8;;http://conflicted.r-lib.org/conflicted package]8;; to force all conflicts to become errors

``` r
library(apeglm)
library(ggplot2)
library(vsn)
library(pheatmap)
library(ReportingTools)
```

    ## Loading required package: knitr
    ## 
    ## Registered S3 method overwritten by 'GGally':
    ##   method from   
    ##   +.gg   ggplot2

``` r
setwd('~/Desktop/510makeup/')

count_matrix <- read.delim('~/Desktop/510makeup/raw/final_genes.txt')

count_matrix <- head(count_matrix, -4)
```

``` r
row.names(count_matrix) <- count_matrix$gene_id
count_matrix <- count_matrix[-c(1)]

sampletable <- read_tsv('~/Desktop/510makeup/raw/race_table.tsv')
```

    ## Rows: 23 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): sample_ID, race
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
row.names(sampletable) <- sampletable$sample_ID
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
sampletable$race <- as.factor(sampletable$race)

DES_dataset <- DESeqDataSetFromMatrix(countData = count_matrix,
                                      colData = sampletable,
                                      design = ~ race)
```

``` r
#pre-filtering
nrow(DES_dataset)
```

    ## [1] 60656

``` r
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

DES_dataset <- DES_dataset[rowSums(counts(DES_dataset)) > 10, ]

# Number of gene after filtering
nrow(DES_dataset)
```

    ## [1] 49354

``` r
DES_dataset <- DESeq(DES_dataset)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 4261 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing
