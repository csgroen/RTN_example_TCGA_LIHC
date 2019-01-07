RTNsurvival example - Hepatocellular carcinoma (TCGA-LIHC)
================
Clarice Groeneveld, Gordon Robertson, Mauro Castro
3 January 2019

Summary
-------

This example uses the RTN package to compute a transcriptional network for hepatocellular carcinoma using the publically available TCGA-LIHC cohort.

We will show how to download the relevant data from the Genomic Data Commons (GDC) using the TCGA-biolinks package (Colaprico et al. 2016). We will also update the GDC information using molecular features from Ally et al. (2017) (the cohort analysis paper) and J. Liu et al. (2018) (pan-cancer clinical analysis of TCGA cohorts).

After joining data from these sources, we will compute the transcriptional regulatory network using a list of 807 transcription factors (Carro et al. 2010) using the RTN implementation of mutual-information and the ARACNE algorithm.

Finally, we will find regulons informative of 5-year Overall Survival by calculating regulon activity using the RTNsurvival package.

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

The current survival data from GDC is not the most recent or complete data. For survival data, we suggest using the data from the Pan-Cancer Clinical publication (J. Liu et al. 2018).

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

Since we'd like to show some molecular features as examples, we'll access the TCGA-LIHC cohort publication's (Ally et al. 2017) supplements. Here, we want to plot the `mRNA clusters` derived for the cohort, but we could use many other interesting molecular features.

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

Pre-process list of regulatory elements (transcription factors)
---------------------------------------------------------------

We'll use the same 807 transcription factors (TFs) used by Castro et. al. (M. A. A. Castro et al. 2016) for a breast cancer transcriptional network reconstruction. Since the annotation available does not use Ensembl Gene IDs, we'll match the TFs by HGNC symbol and find the Ensembl Gene IDs.

``` r
#-- Get list of regulatory elements
load("data/transcriptionFactors.RData")
tfSymbols <- tfs$SYMBOL

#-- Conform regulatory elements annotation 
tfEnsembls <- rowData(tcgaLIHCdata) %>%
    as.data.frame() %>%
    filter(SYMBOL %in% tfSymbols) %>%
    pull(ENSEMBL)
save(tfEnsembls, file = "results/tfEnsembls.RData")
```

References
----------

Ally, Adrian, Miruna Balasundaram, Rebecca Carlsen, Eric Chuah, Amanda Clarke, Noreen Dhalla, Robert A. Holt, et al. 2017. “Comprehensive and Integrative Genomic Characterization of Hepatocellular Carcinoma.” *Cell* 169 (7). Elsevier BV: 1327–1341.e23. doi:[10.1016/j.cell.2017.05.046](https://doi.org/10.1016/j.cell.2017.05.046).

Carro, Maria Stella, Wei Keat Lim, Mariano Javier Alvarez, Robert J. Bollo, Xudong Zhao, Evan Y. Snyder, Erik P. Sulman, et al. 2010. “The Transcriptional Network for Mesenchymal Transformation of Brain Tumours.” *Nature* 463 (7279). Springer Nature: 318–25. doi:[10.1038/nature08712](https://doi.org/10.1038/nature08712).

Castro, Mauro A A, Ines de Santiago, Thomas M Campbell, Courtney Vaughn, Theresa E Hickey, Edith Ross, Wayne D Tilley, Florian Markowetz, Bruce A J Ponder, and Kerstin B Meyer. 2016. “Regulators of Genetic Risk of Breast Cancer Identified by Integrative Network Analysis.” *Nature Genetics* 48 (1). Springer Nature: 12–21. doi:[10.1038/ng.3458](https://doi.org/10.1038/ng.3458).

Colaprico, Antonio, Tiago C. Silva, Catharina Olsen, Luciano Garofano, Claudia Cava, Davide Garolini, Thais S. Sabedot, et al. 2016. “TCGAbiolinks: An R/Bioconductor Package for Integrative Analysis of Tcga Data.” *Nucleic Acids Research* 44 (8): e71. doi:[10.1093/nar/gkv1507](https://doi.org/10.1093/nar/gkv1507).

Liu, Jianfang, Tara Lichtenberg, Katherine A. Hoadley, Laila M. Poisson, Alexander J. Lazar, Andrew D. Cherniack, Albert J. Kovatich, et al. 2018. “An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics.” *Cell* 173 (2). Elsevier BV: 400–416.e11. doi:[10.1016/j.cell.2018.02.052](https://doi.org/10.1016/j.cell.2018.02.052).
