RTNsurvival example - Hepatocellular carcinoma (TCGA-LIHC)
================
Clarice Groeneveld, Gordon Robertson, Mauro Castro <br>
28 January 2019

Summary
-------

This example shows how to prepare input data for the RTN and RTNsurvival packages using the publicly available mRNA-seq data for the TCGA-LIHC cohort. We will show how to download the relevant harmonized GRCh38/hg38 data from the Genomic Data Commons (GDC) using the TCGAbiolinks package (Colaprico *et al.*, 2016). We will also update the GDC information using molecular features from the TCGA LIHC analysis publication (The Cancer Genome Atlas Research Network, 2017) and outcomes from the Pan-Cancer Atlas clinical data publication (Liu *et al.*, 2018). The preprocessing generates a `SummarizedExperiment` object that contains expression, clinical and molecular data.

Libraries
---------

Please make sure you have all libraries installed before proceeding. Also, please set a working directory and download the relevant `data`.

``` r
library(RTNsurvival)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(plyr)
library(tidyverse)
library(snow)
library(readxl)
library(caret)
library(knitr)
```

``` r
#-- Check to make sure the data directory is in place
if (!dir.exists("data") || !file.exists("data/transcriptionFactors.RData")) {
    stop("-- NOTE: Please make sure to download the relevant `data` directory and place it into the working directory.")
}
```

Using TCGAbiolinks to download data from GDC
--------------------------------------------

We'll use the Bioconductor package `TCGAbiolinks` to query and download from GDC. We are looking for the harmonized, pre-processed RNA-seq for the TCGA-LIHC cohort. `TCGAbiolinks` will create a directory called GDCdata in your working directory and save the files downloaded from GDC. The download can take a while. If you find that the download has failed, feel free to change the `files.per.chunk` argument. The files for each patient will be downloaded in a separate file. Then, the GDCprepare function will compile them into an R object of class `RangedSummarizedExperiment`.

The `RangedSummarizedExperiment` has 6 slots. The most important being the `rowRanges` (gene metadata), `colData` (patient metadata), and `assays`, which contains the gene expression matrix.

``` r
#-- Download the TCGA-LIHC data
query <- GDCquery("TCGA-LIHC", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - FPKM")
GDCdownload(query)
tcgaLIHCdata <- GDCprepare(query)
```

The object downloaded from GDC has 56716 rows for different transcripts. This list includes both coding and lincRNAs. For best RTN performance, we will filter these genes to maintain only genes annotated in the UCSC hg38 known gene list.

``` r
#-- Subset by known gene locations
geneRanges <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
tcgaLIHCdata <- subsetByOverlaps(tcgaLIHCdata, geneRanges)
```

Then, we'll select only primary tumor data and unique patients.

``` r
#-- Select unique primary tumors
tcgaLIHCdata <- subset(tcgaLIHCdata, select = definition == "Primary solid Tumor")
tcgaLIHCdata <- subset(tcgaLIHCdata, select = isUnique(patient))
```

Finally, we'll perform a one-line pre-process step for better internal pre-processing in RTN's `tni.constructor` function. Having the `"SYMBOL"` column will enable genes with the same symbol to be preprocessed by this function.

``` r
#-- Change column names for best tni.constructor performance
colnames(rowData(tcgaLIHCdata)) <- c("ENSEMBL", "SYMBOL", "OG_ENSEMBL")

#-- Get patient metadata for additional pre-processing
colAnnotation <- colData(tcgaLIHCdata)
```

The `tcgaLIHCdata` object would be ready for the RTN pipeline. But since we want to assess survival and molecular features, we will bring in other annotation sources on the TCGA cohort to complement the patient metadata available in GDC.

Download survival data from the TCGA Pan-Cancer Clinical Data paper (Liu et. al, 2018)
--------------------------------------------------------------------------------------

The current survival data from GDC is not the most recent or complete data. For survival data, we suggest using the data from the Pan-Cancer Clinical publication (Liu *et al.*, 2018).

We can download the data directly from the supplements of the publication. We'll use the `readxl` package to read the data, and `tidyverse` functions for filtering the data, remapping and renaming features and deriving new features from the provided features.

``` r
#-- Get more complete survival data (Liu et. al. 2018)
download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx",
              destfile = "data/Liu2018_survivalData.xlsx")

#-- Read survival data
nms <- names(read_excel("data/Liu2018_survivalData.xlsx", n_max = 0))
col_types <- c("skip", ifelse(grepl("residual_tumor|cause_of_death", nms), "text", "guess"))
liu_survData <- read_excel("data/Liu2018_survivalData.xlsx", sheet = 1, 
                           col_types = col_types, na = c("#N/A", "[Not Available]"))

#-- Preprocess
lihc_survData <- liu_survData %>%
    filter(type == "LIHC", 
           bcr_patient_barcode %in% colAnnotation$patient) %>%
    select(bcr_patient_barcode, gender, Age = age_at_initial_pathologic_diagnosis,
           ajcc_pathologic_tumor_stage, OS, OS.time, PFI, PFI.time) %>%
    mutate(OS.time.months = OS.time/30,
           PFI.time.months = PFI.time/30,
           Tumor_Stage = as.numeric(mapvalues(ajcc_pathologic_tumor_stage, 
                         from = c("Stage I", "Stage II", "Stage III",
                                  "Stage IIIA", "Stage IIIB", "Stage IIIC", 
                                  "Stage IV", "Stage IVA", "Stage IVB",
                                  "[Discrepancy]"),
                         to = c(1, 2, 3, 3, 3, 3, 4, 4, 4, NA))),
           Stage = as.factor(mapvalues(Tumor_Stage,
                   from = c(1,2,3,4), to = c("I", "II", "III", "IV")))) %>%
    as.data.frame()

#-- Dummy encode Stage for visualization
dummyStage <- dummyVars(~ Stage, lihc_survData, sep = "_") %>% 
    predict(lihc_survData)
lihc_survData <- cbind(lihc_survData, dummyStage)
```

Download molecular data from the TCGA-LIHC cohort paper
-------------------------------------------------------

Since we'd like to show some molecular features as examples, we'll access the TCGA-LIHC cohort publication's (The Cancer Genome Atlas Research Network, 2017) supplements. Here, we want to plot the `mRNA clusters` derived for the cohort, but we could use many other interesting molecular features.

``` r
#-- Read molecular covariates from the cohort paper
download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867417306396-mmc1.xlsx",
              destfile = "data/TCGA_coreClinicalMolecular.xlsx")

lihc_molcData <- read_excel("data/TCGA_coreClinicalMolecular.xlsx", range = "A4:CT200")

lihc_molcData <- lihc_molcData %>%
    select(Barcode, mRNA = `mRNA clusters (5 group NMF, Hoadley group)`) %>%
    transmute(barcode = str_sub(Barcode, 1, 12),
              mRNA = mapvalues(mRNA, from = "NA", to = NA))

#-- Dummy encode the mRNA cluster variable from the core samples
dummy_mRNA <- dummyVars(~ mRNA, lihc_molcData) %>% predict(lihc_molcData)
lihc_molcData <- cbind(lihc_molcData, dummy_mRNA)

#-- Join molecular features to clinical features
lihc_survData <- left_join(lihc_survData, lihc_molcData, 
                           by = c("bcr_patient_barcode" = "barcode"))
```

Join molecular and clinical features to SummarizedExperiment
------------------------------------------------------------

Finally, we'll join the features we derived from the two Cell publications to the data we got from GDC through TCGAbiolinks.

``` r
#-- Conform names
idx <- match(lihc_survData$bcr_patient_barcode, colAnnotation$patient)
rownames(lihc_survData) <- rownames(colAnnotation)[idx]

#-- Add back to SummarizedExperiment
colData(tcgaLIHCdata) <- as(lihc_survData, "DataFrame")

dir.create("results")
```

    ## Warning in dir.create("results"): 'results' already exists

``` r
save(tcgaLIHCdata, file = "results/tcgaLIHCdata_preprocessed.RData")
```

References
----------

Colaprico,A. *et al.* (2016) TCGAbiolinks: An r/bioconductor package for integrative analysis of tcga data. *Nucleic Acids Research*, **44**, e71.

Liu,J. *et al.* (2018) An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. *Cell*, **173**, 400–416.e11.

The Cancer Genome Atlas Research Network (2017) Comprehensive and integrative genomic characterization of hepatocellular carcinoma. *Cell*, **169**, 1327–1341.e23.
