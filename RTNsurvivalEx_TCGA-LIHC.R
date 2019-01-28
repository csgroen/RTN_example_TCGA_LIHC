################################################################################
### RTNsurvival example - TCGA-LIHC (pre-process)
################################################################################

#-- Call/Install libraries
# install with BiocManager:install("LIBRARYNAME")
# Works for packages from CRAN or Bioconductor
library(SummarizedExperiment)
library(TCGAbiolinks)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(plyr)
library(tidyverse)
library(readxl)
library(caret)
library(knitr)

#-------------------------------------------------------------------------------
#-- Download the TCGA-LIHC data
query <- GDCquery("TCGA-LIHC", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - FPKM")
GDCdownload(query)
tcgaLIHCdata <- GDCprepare(query)

#-- Subset by known gene locations
geneRanges <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
tcgaLIHCdata <- subsetByOverlaps(tcgaLIHCdata, geneRanges)

#-- Select unique primary tumors
tcgaLIHCdata <- subset(tcgaLIHCdata, select = definition == "Primary solid Tumor")
tcgaLIHCdata <- subset(tcgaLIHCdata, select = isUnique(patient))

#-- Change column names for best tni.constructor performance
colnames(rowData(tcgaLIHCdata)) <- c("ENSEMBL", "SYMBOL", "OG_ENSEMBL")

#-- Get patient metadata for additional pre-processing
colAnnotation <- colData(tcgaLIHCdata)

#-------------------------------------------------------------------------------
#-- Get more complete survival data (Liu et. al. 2018)
download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx",
              destfile = "Liu2018_survivalData.xlsx")

#-- Read survival data
nms <- names(read_excel("Liu2018_survivalData.xlsx", n_max = 0))
col_types <- c("skip", ifelse(grepl("residual_tumor|cause_of_death", nms), "text", "guess"))
liu_survData <- read_excel("Liu2018_survivalData.xlsx", sheet = 1, 
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

#-------------------------------------------------------------------------------
#-- Read molecular covariates from the cohort paper
download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867417306396-mmc1.xlsx",
              destfile = "TCGA_coreClinicalMolecular.xlsx")

lihc_molcData <- read_excel("TCGA_coreClinicalMolecular.xlsx", range = "A4:CT200")

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

#-------------------------------------------------------------------------------
#-- Conform names
rownames(lihc_survData) <- lihc_survData$bcr_patient_barcode
lihc_survData <- lihc_survData[colAnnotation$patient,]
rownames(lihc_survData) <- rownames(colAnnotation)

#-- Add back to SummarizedExperiment
colData(tcgaLIHCdata) <- as(lihc_survData, "DataFrame")

save(tcgaLIHCdata, file = "tcgaLIHCdata_preprocessed.RData")
