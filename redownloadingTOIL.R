library(UCSCXenaTools)
library(dplyr)
library(data.table)
library(TCGAbiolinks)

setwd("C:/Users/nputn/Google Drive/UCLA Research/")

# Getting the IDs for all of the UVM samples from TCGA using TCGAbiolinks
query <- GDCquery(project = "TCGA-UVM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")

results <- getResults(query)
sampleID <- results$sample.submitter_id

# Getting TOILed gene expression data using UCSCXenaTools
temp <- XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "PANCAN")

toDownload <- XenaFilter(temp, filterDatasets = "tcga_RSEM_Hugo_norm_count")


####  only run if data isn't downloaded #####################################
XenaQuery(toDownload) %>% #(only run)
  XenaDownload() -> xe_download

############################################################################


# if file exists, just read in
normExpData <- as.data.frame(fread("tcga_RSEM_Hugo_norm_count"))

# Finding intersect of all UVM gene expression samples and all of the TOILed gene expression samples
sampleID <- sapply(sampleID, substr, 1, 15)
sampIDInter <- intersect(colnames(normExpData), sampleID)

# Subsetting TOIL'd data using UVM sample names
UVMExp <- normExpData[colnames(normExpData) %in% c("sample", sampIDInter)]
colnames(UVMExp) <- c("gene", colnames(UVMExp[2:ncol(UVMExp)]))

