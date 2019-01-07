Summary
=======

This example uses the RTN package to compute a transcriptional network for hepatocellular carcinoma using the publically available TCGA-LIHC cohort.

We will show how to download the relevant data from the Genomic Data Commons (GDC) using the TCGA-biolinks package (Colaprico et al. 2016). We will also update the GDC information using molecular features from Ally et al. (2017) (the cohort analysis paper) and J. Liu et al. (2018) (pan-cancer clinical analysis of TCGA cohorts).

After joining data from these sources, we will compute the transcriptional regulatory network using a list of 807 transcription factors (Carro et al. 2010) using the RTN implementation of mutual-information and the ARACNE algorithm.

Finally, we will find regulons informative of 5-year Overall Survival by calculating regulon activity using the RTNsurvival package.

Libraries
=========

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
============================================

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
======================================================================================

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
=======================================================

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
============================================================

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
===============================================================

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

Session information
===================

    ## R version 3.5.2 (2018-12-20)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.1 LTS
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] bindrcpp_0.2.2                         
    ##  [2] knitr_1.21                             
    ##  [3] caret_6.0-81                           
    ##  [4] lattice_0.20-38                        
    ##  [5] readxl_1.1.0                           
    ##  [6] snow_0.4-3                             
    ##  [7] forcats_0.3.0                          
    ##  [8] stringr_1.3.1                          
    ##  [9] dplyr_0.7.8                            
    ## [10] purrr_0.2.5                            
    ## [11] readr_1.3.0                            
    ## [12] tidyr_0.8.2                            
    ## [13] tibble_1.4.2                           
    ## [14] ggplot2_3.1.0                          
    ## [15] tidyverse_1.2.1                        
    ## [16] plyr_1.8.4                             
    ## [17] TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0
    ## [18] GenomicFeatures_1.34.1                 
    ## [19] AnnotationDbi_1.44.0                   
    ## [20] TCGAbiolinks_2.10.0                    
    ## [21] SummarizedExperiment_1.12.0            
    ## [22] DelayedArray_0.8.0                     
    ## [23] BiocParallel_1.16.2                    
    ## [24] matrixStats_0.54.0                     
    ## [25] Biobase_2.42.0                         
    ## [26] GenomicRanges_1.34.0                   
    ## [27] GenomeInfoDb_1.18.1                    
    ## [28] IRanges_2.16.0                         
    ## [29] S4Vectors_0.20.1                       
    ## [30] BiocGenerics_0.28.0                    
    ## [31] RTNsurvival_1.6.0                      
    ## [32] RTNduals_1.6.0                         
    ## [33] RTN_2.7.1                              
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] backports_1.1.2             circlize_0.4.5             
    ##   [3] aroma.light_3.12.0          igraph_1.2.2               
    ##   [5] selectr_0.4-1               ConsensusClusterPlus_1.46.0
    ##   [7] lazyeval_0.2.1              splines_3.5.2              
    ##   [9] sva_3.30.0                  digest_0.6.18              
    ##  [11] foreach_1.4.4               htmltools_0.3.6            
    ##  [13] magrittr_1.5                memoise_1.1.0              
    ##  [15] cluster_2.0.7-1             doParallel_1.0.14          
    ##  [17] RedeR_1.30.0                mixtools_1.1.0             
    ##  [19] limma_3.38.3                recipes_0.1.4              
    ##  [21] ComplexHeatmap_1.20.0       Biostrings_2.50.1          
    ##  [23] annotate_1.60.0             gower_0.1.2                
    ##  [25] modelr_0.1.2                R.utils_2.7.0              
    ##  [27] prettyunits_1.0.2           colorspace_1.3-2           
    ##  [29] blob_1.1.1                  rvest_0.3.2                
    ##  [31] ggrepel_0.8.0               haven_2.0.0                
    ##  [33] xfun_0.4                    crayon_1.3.4               
    ##  [35] RCurl_1.95-4.11             jsonlite_1.6               
    ##  [37] genefilter_1.64.0           bindr_0.1.1                
    ##  [39] zoo_1.8-4                   survival_2.43-1            
    ##  [41] iterators_1.0.10            glue_1.3.0                 
    ##  [43] survminer_0.4.3             gtable_0.2.0               
    ##  [45] ipred_0.9-8                 zlibbioc_1.28.0            
    ##  [47] XVector_0.22.0              GetoptLong_0.1.7           
    ##  [49] shape_1.4.4                 scales_1.0.0               
    ##  [51] DESeq_1.34.0                DBI_1.0.0                  
    ##  [53] edgeR_3.24.1                ggthemes_4.0.1             
    ##  [55] Rcpp_1.0.0                  xtable_1.8-3               
    ##  [57] progress_1.2.0              cmprsk_2.2-7               
    ##  [59] bit_1.1-14                  matlab_1.0.2               
    ##  [61] km.ci_0.5-2                 lava_1.6.4                 
    ##  [63] prodlim_2018.04.18          httr_1.4.0                 
    ##  [65] RColorBrewer_1.1-2          pkgconfig_2.0.2            
    ##  [67] XML_3.98-1.16               R.methodsS3_1.7.1          
    ##  [69] nnet_7.3-12                 locfit_1.5-9.1             
    ##  [71] reshape2_1.4.3              tidyselect_0.2.5           
    ##  [73] rlang_0.3.0.1               cellranger_1.1.0           
    ##  [75] munsell_0.5.0               tools_3.5.2                
    ##  [77] cli_1.0.1                   downloader_0.4             
    ##  [79] generics_0.0.2              RSQLite_2.1.1              
    ##  [81] broom_0.5.1                 evaluate_0.12              
    ##  [83] yaml_2.2.0                  ModelMetrics_1.2.2         
    ##  [85] bit64_0.9-7                 survMisc_0.5.5             
    ##  [87] EDASeq_2.16.0               nlme_3.1-137               
    ##  [89] R.oo_1.22.0                 xml2_1.2.0                 
    ##  [91] biomaRt_2.38.0              rstudioapi_0.8             
    ##  [93] compiler_3.5.2              e1071_1.7-0                
    ##  [95] minet_3.40.0                viper_1.16.0               
    ##  [97] geneplotter_1.60.0          stringi_1.2.4              
    ##  [99] Matrix_1.2-15               KMsurv_0.1-5               
    ## [101] pillar_1.3.0                GlobalOptions_0.1.0        
    ## [103] data.table_1.11.8           bitops_1.0-6               
    ## [105] rtracklayer_1.42.1          R6_2.3.0                   
    ## [107] latticeExtra_0.6-28         hwriter_1.3.2              
    ## [109] ShortRead_1.40.0            KernSmooth_2.23-15         
    ## [111] gridExtra_2.3               codetools_0.2-16           
    ## [113] MASS_7.3-51.1               assertthat_0.2.0           
    ## [115] rjson_0.2.20                withr_2.1.2                
    ## [117] GenomicAlignments_1.18.0    Rsamtools_1.34.0           
    ## [119] GenomeInfoDbData_1.2.0      mgcv_1.8-26                
    ## [121] hms_0.4.2                   rpart_4.1-13               
    ## [123] timeDate_3043.102           grid_3.5.2                 
    ## [125] class_7.3-14                rmarkdown_1.11             
    ## [127] segmented_0.5-3.0           ggpubr_0.2                 
    ## [129] lubridate_1.7.4             rematch_1.0.1

References
==========

Ally, Adrian, Miruna Balasundaram, Rebecca Carlsen, Eric Chuah, Amanda Clarke, Noreen Dhalla, Robert A. Holt, et al. 2017. “Comprehensive and Integrative Genomic Characterization of Hepatocellular Carcinoma.” *Cell* 169 (7). Elsevier BV: 1327–1341.e23. doi:[10.1016/j.cell.2017.05.046](https://doi.org/10.1016/j.cell.2017.05.046).

Carro, Maria Stella, Wei Keat Lim, Mariano Javier Alvarez, Robert J. Bollo, Xudong Zhao, Evan Y. Snyder, Erik P. Sulman, et al. 2010. “The Transcriptional Network for Mesenchymal Transformation of Brain Tumours.” *Nature* 463 (7279). Springer Nature: 318–25. doi:[10.1038/nature08712](https://doi.org/10.1038/nature08712).

Castro, Mauro A A, Ines de Santiago, Thomas M Campbell, Courtney Vaughn, Theresa E Hickey, Edith Ross, Wayne D Tilley, Florian Markowetz, Bruce A J Ponder, and Kerstin B Meyer. 2016. “Regulators of Genetic Risk of Breast Cancer Identified by Integrative Network Analysis.” *Nature Genetics* 48 (1). Springer Nature: 12–21. doi:[10.1038/ng.3458](https://doi.org/10.1038/ng.3458).

Colaprico, Antonio, Tiago C. Silva, Catharina Olsen, Luciano Garofano, Claudia Cava, Davide Garolini, Thais S. Sabedot, et al. 2016. “TCGAbiolinks: An R/Bioconductor Package for Integrative Analysis of Tcga Data.” *Nucleic Acids Research* 44 (8): e71. doi:[10.1093/nar/gkv1507](https://doi.org/10.1093/nar/gkv1507).

Liu, Jianfang, Tara Lichtenberg, Katherine A. Hoadley, Laila M. Poisson, Alexander J. Lazar, Andrew D. Cherniack, Albert J. Kovatich, et al. 2018. “An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics.” *Cell* 173 (2). Elsevier BV: 400–416.e11. doi:[10.1016/j.cell.2018.02.052](https://doi.org/10.1016/j.cell.2018.02.052).
