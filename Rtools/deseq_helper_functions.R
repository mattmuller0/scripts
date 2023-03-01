###########################################################################
#
#                                 HEADER
#
###########################################################################
# Author: Matthew Muller
# 
# Date: 2023-02-18
# 
# Script Name: DESeq2 Helper Functions
# 
# Notes:
# This script contains helper functions that can be used to assist in differential
# expression analysis using DESeq2. Most inputs are expected as DDS or SE objects.


###########################################################################
#
#                                 LIBRARIES
#
###########################################################################
packages <- c(
  "tidyverse",
  "ggplot2",
  "devtools",
  "BiocManager",
  "SummarizedExperiment",
  "DESeq2",
  "edgeR",
  "AnnotationDbi",
  "apeglm"
)

for (pkg in packages) {
  paste0(pkg)[[1]]
  library(pkg, character.only = T, quietly = T)
}

# LOAD FUNCTIONS
# space reserved for sourcing in functions
options(stringsAsFactors = FALSE)


###########################################################################
#
#                                 CODE
#
###########################################################################
# Add code here
#
ovr_deseq_results <- function(dds, column, outpath){
  # This function takes in a DDS or SE object, the condition column of interest
  # and the out directory path to push results to. It will run OVR differential
  # expression analysis on each level within the condition column.
  counts <- assay(dds)
  lvls <- levels(colData(dds)[,column])
  cond <- as.character(colData(dds)[,column])

  # loop over condition levels in a one versus rest manner
  # ideally this for loop could be an apply statement with a custom function, 
  # but it works for now so oh well lol
  for (lvl in lvls) {
    print(paste0('Testing ', column, ' ', lvl, ' versus rest'))
    # Set our OVR analysis
    cond_ <- cond
    cond_[cond_ != lvl] <- 'rest'
    
    # create temporary dds objects for analysis
    dds_ <- DESeqDataSetFromMatrix(counts, colData <- DataFrame(condition = as.factor(cond_)), 
                                   design <- ~ condition)
    dds_ <- DESeq(dds_)
    res <- results(dds_)
    write.csv(dds_, file=paste0(outpath, '/deseqDataset_', column,'__',lvl,'_v_','rest.csv'))
    write.csv(res, file=paste0(outpath, '/dge_results_', column,'__',lvl,'_v_','rest.csv'))
    }
  }