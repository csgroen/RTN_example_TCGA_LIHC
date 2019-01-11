Example of data preprocessing for RTNsurvival, using the TCGA-LIHC cohort
================
Clarice Groeneveld, Gordon Robertson, Mauro Castro<br>
10 January 2019

Summary
=======

This example uses the RTN package to compute a transcriptional regulatory network for hepatocellular carcinoma using the publicly available mRNA-seq data for the TCGA-LIHC cohort.

We will show how to download the relevant harmonized GRCh38/hg38 data from the Genomic Data Commons (GDC) using the TCGAbiolinks package (Colaprico et al. 2016). We will also update the GDC information using molecular features from the TCGA LIHC analysis publication The Cancer Genome Atlas Research Network (2017) and outcomes from the Pan-Cancer Atlas clinical data publication J. Liu et al. (2018).

After joining data from these sources, we will compute the mutual information-based transcriptional regulatory network, using the RTN implementation of the ARACNE algorithm (Margolin et al. 2006). We will then assess the regulons for a list of 807 transcription factors (Carro et al. 2010) as regulators.

Finally, we will identify regulons whose activity is informative of 5-year Overall Survival, using the RTNsurvival package.

Libraries
=========

Please make sure you have all libraries installed before proceeding. Also, please set a working directory.

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
#-- Check to make sure the data directory exists in the working directory
if (!dir.exists("data") || !file.exists("data/transcriptionFactors.RData")) {
    stop("-- NOTE: Please make sure to download the relevant `data` directory and place it into the working directory.")
}
```

Use TCGAbiolinks to download harmonized data from the GDC
=========================================================

We'll use the Bioconductor package `TCGAbiolinks` to query and download from the GDC. We are looking for the harmonized, pre-processed RNA-seq for the TCGA-LIHC cohort. `TCGAbiolinks` will create a directory called GDCdata in your working directory and save the files downloaded from the GDC. The files for each patient will be downloaded as a separate file. For LIHC, there are 424 files (~200 MB total), so the download can take a while. The expression profile for each patient will be downloaded into a separate folder and file. Then, the GDCprepare function will compile them into an R object of class `RangedSummarizedExperiment`.

The `RangedSummarizedExperiment` has 6 slots. The most important slots are `rowRanges` (gene metadata), `colData` (patient metadata), and `assays` (gene expression matrix).

``` r
#-- Download the TCGA-LIHC data
query <- GDCquery("TCGA-LIHC", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - FPKM")
GDCdownload(query)
tcgaLIHCdata <- GDCprepare(query)
```

The object downloaded from the GDC has 56,716 gene records. This list includes both coding and noncoding genes, e.g. lincRNAs. We will filter these genes to retain only genes annotated in the UCSC hg38 known gene list.

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

Finally, we'll perform a one-line step for better internal pre-processing in RTN's `tni.constructor` function. Having the `"SYMBOL"` column will enable duplicated genes (those that share a HGNC symbol) to be preprocessed by this function.

``` r
#-- Change column names for best tni.constructor performance
colnames(rowData(tcgaLIHCdata)) <- c("ENSEMBL", "SYMBOL", "OG_ENSEMBL")

#-- Get patient metadata for additional pre-processing
colAnnotation <- colData(tcgaLIHCdata)
```

The `tcgaLIHCdata` object is now ready for the RTN pipeline. But, since we want to assess survival and molecular features, we will bring in other annotation sources on the TCGA cohort to complement the patient metadata available in the GDC.

Download survival data from the TCGA Pan-Cancer Clinical Data paper
===================================================================

While survival data is available from the GDC, we suggest using outcomes from the Pan-Cancer clinical data publication (J. Liu et al. 2018). Table 3 in that publication indicates that, for the LIHC cohort (n=377), OS, PFI, and DFI data are directly usable, but cautions that DSS needs a longer follow-up. Here, we will use the OS data.

We can download the data directly from the supplements of the publication. We'll use the `readxl` package to read the data, and `tidyverse` functions for filtering the data, remapping and renaming features and deriving new features from the provided features.

``` r
#-- Get Pan-Cancer survival data (Liu et. al. 2018)
download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx",
              destfile = "data/Liu2018_survivalData.xlsx")

#-- Read survival data
nms <- names(read_excel("data/Liu2018_survivalData.xlsx", n_max = 0))
```

    ## readxl works best with a newer version of the tibble package.
    ## You currently have tibble v1.4.2.
    ## Falling back to column name repair from tibble <= v1.4.2.
    ## Message displays once per session.

``` r
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
```

In order to use the tnsPlotKM method, covariate tracks must be binary variables. So, we need to transform (i.e. ‘dummy-encode’) categorical variable of interest, like Tumor\_Stage (I, II, III, IV), into a set of binary TRUE/FALSE variables.

``` r
#-- Dummy encode Stage for visualization
dummyStage <- dummyVars(~ Stage, lihc_survData, sep = "_") %>% 
    predict(lihc_survData)
lihc_survData <- cbind(lihc_survData, dummyStage)
```

Download molecular data from the TCGA-LIHC cohort publication
=============================================================

Since we'd like to show some molecular features as examples, we'll access the TCGA-LIHC cohort publication's supplements (The Cancer Genome Atlas Research Network 2017). Here, we want to plot the `mRNA clusters` derived for the cohort, but we could also assess many other interesting molecular features.

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

Join clinical and molecular features to RangedSummarizedExperiment object
=========================================================================

Finally, we'll join the features we derived from the two Cell publications to the data we got from the GDC through TCGAbiolinks.

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

We'll woth with the same 807 transcription factors (TFs) used by Castro et al. (2016) with a breast cancer transcriptional network reconstruction. Since the TCGA annotations do not use Ensembl Gene IDs, we'll match the TFs by HGNC symbol and then find the Ensembl Gene IDs.

``` r
#-- Get list of regulatory elements
load("data/transcriptionFactors.RData")
tfSymbols <- tfs$SYMBOL

#-- Add ENSEMBL information to regulatory element annotations
tfEnsembls <- rowData(tcgaLIHCdata) %>%
    as.data.frame() %>%
    filter(SYMBOL %in% tfSymbols) %>%
    pull(ENSEMBL)
save(tfEnsembls, file = "results/tfEnsembls.RData")
```

Inference of the regulatory network with RTN
--------------------------------------------

The *RTN* pipeline can be run in single- or multi-threaded modes. To run in multithreaded, we simply load the `parallel` library, detect available cores and define the cluster before initiating the pipeline functions.

``` r
#-- OPTIONAL - for parallel processing
library(parallel)
n <- detectCores() - 1
options(cluster=makeCluster(n, "SOCK"))
```

The `tni.constructor` method takes in a matrix of gene expression and metadata on the samples and genes, as well as a vector of regulatory elements to be evaluated. In this case, the expression matrix and metadata are enclosed as a `RangedSummarizedExpression` object. The ‘tni.constructor’ method also performs pre-processing to check the consistency of all given arguments and maximize algorithm performance. It returns a `TNI` (Transcriptional Network - Inference) object.

``` r
#-- Construct a TNI object
lihcTNI <- tni.constructor(tcgaLIHCdata, regulatoryElements = tfEnsembls)
```

The next step uses the `tni.permutation` method to estimate the statistical significance of calculated mutual information (MI) between each regulatory element and each of its potential targets. This method performs permutation analysis to construct a null distribution of MI values, in order to assess the significance of calculated an MI value between a regulator and a target. Using only the significant regulator-target pairs, the method constructs the **reference network**. In our example, only regulator-target pairs that have a p &lt; 10<sup>-5</sup> are considered significant.

``` r
#-- Permutation analysis
lihcTNI <- tni.permutation(lihcTNI, pValueCutoff = 10^-5, estimator = "spearman")
```

The `tni.bootstrap` takes the reference network and performs bootstrap analysis to assess the stability of the inferred regulator-target pairs. In this example, we perform 200 bootstraps and keep the pairs that are observed in at least 95% of the iterations (i.e. a 95% consensus).

``` r
#-- Bootstrap analysis
lihcTNI <- tni.bootstrap(lihcTNI, nBootstraps = 200)
```

The `tni.dpi.filter` performs the Data Processing Inequality filter on the reference network. This filter looks at relationships where two regulators have significant mutual information with each other and also with a common target. Those triplets are broken at the weakest link (i.e. the edge with the lowest MI), to minimize the number of indirect interactions present in the network. This filtered network is called the **DPI network** (or DPI-filtered network). Please refer to Margolin et al. (2006), Fletcher et al. (2013), Castro et al. (2016) and Robertson et al. (2017) for additional details.

``` r
#-- Data Processing Inequality filter
lihcTNI <- tni.dpi.filter(lihcTNI)
```

If run in parallel mode/multithreaded, you should disable the cluster now.

``` r
#-- OPTIONAL - for parallel processing
stopCluster(getOption("cluster"))
```

References
==========

Carro, Maria Stella, Wei Keat Lim, Mariano Javier Alvarez, Robert J. Bollo, Xudong Zhao, Evan Y. Snyder, Erik P. Sulman, et al. 2010. “The Transcriptional Network for Mesenchymal Transformation of Brain Tumours.” *Nature* 463 (7279). Springer Nature: 318–25. doi:[10.1038/nature08712](https://doi.org/10.1038/nature08712).

Castro, M A A, Ines de Santiago, Thomas M Campbell, Courtney Vaughn, Theresa E Hickey, Edith Ross, Wayne D Tilley, Florian Markowetz, Bruce A J Ponder, and Kerstin B Meyer. 2016. “Regulators of Genetic Risk of Breast Cancer Identified by Integrative Network Analysis.” *Nature Genetics* 48 (1). Springer Nature: 12–21. doi:[10.1038/ng.3458](https://doi.org/10.1038/ng.3458).

Colaprico, Antonio, Tiago C. Silva, Catharina Olsen, Luciano Garofano, Claudia Cava, Davide Garolini, Thais S. Sabedot, et al. 2016. “TCGAbiolinks: An R/Bioconductor Package for Integrative Analysis of Tcga Data.” *Nucleic Acids Research* 44 (8): e71. doi:[10.1093/nar/gkv1507](https://doi.org/10.1093/nar/gkv1507).

Fletcher, M N, M A Castro, X Wang, I de Santiago, M O’Reilly, S F Chin, O M Rueda, et al. 2013. “Master regulators of FGFR2 signalling and breast cancer risk.” *Nature Communications* 4: 2464. doi:[10.1038/ncomms3464](https://doi.org/10.1038/ncomms3464).

Liu, Jianfang, Tara Lichtenberg, Katherine A. Hoadley, Laila M. Poisson, Alexander J. Lazar, Andrew D. Cherniack, Albert J. Kovatich, et al. 2018. “An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics.” *Cell* 173 (2). Elsevier BV: 400–416.e11. doi:[10.1016/j.cell.2018.02.052](https://doi.org/10.1016/j.cell.2018.02.052).

Margolin, Adam A, Ilya Nemenman, Katia Basso, Chris Wiggins, Gustavo Stolovitzky, Riccardo Favera, and Andrea Califano. 2006. “ARACNE: An Algorithm for the Reconstruction of Gene Regulatory Networks in a Mammalian Cellular Context.” *BMC Bioinformatics* 7 (Suppl 1). Springer Nature: S7. doi:[10.1186/1471-2105-7-s1-s7](https://doi.org/10.1186/1471-2105-7-s1-s7).

Robertson, A G, Kim J, Al-Ahmadie H, Bellmunt J, Guo G, Cherniack A D, Hinoue T, et al. 2017. “Comprehensive Molecular Characterization of Muscle-Invasive Bladder Cancer.” *Cell* 171 (3): 540–56. doi:[10.1016/j.cell.2017.09.007](https://doi.org/10.1016/j.cell.2017.09.007).

The Cancer Genome Atlas Research Network. 2017. “Comprehensive and Integrative Genomic Characterization of Hepatocellular Carcinoma.” *Cell* 169 (7). Elsevier BV: 1327–1341.e23. doi:[10.1016/j.cell.2017.05.046](https://doi.org/10.1016/j.cell.2017.05.046).
