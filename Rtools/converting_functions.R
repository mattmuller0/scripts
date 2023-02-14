# HEADER --------------------------------------------------
# Author: Matthew Muller
# 
# Date: 2023-01-15
# 
# Script Name: Converting Functions Script
# 
# Notes:
# 

# SET WORKING DIRECTORY -----------------------------------
wd <- FALSE
if (wd != FALSE) {
  setwd(wd)
  cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
  }

# LOAD LIBRARIES ------------------------------------------
packages <- c(
  "tidyverse",
  "ggplot2"
)

for (pkg in packages) {
  paste0(pkg)[[1]]
  library(pkg, character.only = T, quietly = T)
}

# LOAD FUNCTIONS ------------------------------------------
# space reserved for sourcing in functions



# CODE BLOCK ----------------------------------------------
# Add code here
#

getMatrixWithSelectedIds <- function(df, type, keys){
  require("AnnotationDbi")
  require("org.Hs.eg.db")
  geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(df), column=c(type), keytype=keys, multiVals="first")
  # get the entrez ids with gene symbols i.e. remove those with NA's for gene symbols
  inds <- which(!is.na(geneSymbols))
  found_genes <- geneSymbols[inds]
  # subset your data frame based on the found_genes
  df2 <- df[names(found_genes), ]
  rownames(df2) <- found_genes
  return(df2)
}

mor_normalization <- function(data){
  require(dplyr)
  require(tibble)
  
  # take the log
  log_data = log(data) 
  
  # find the psuedo-references per sample by taking the geometric mean
  log_data = log_data %>% 
    rownames_to_column('gene') %>% 
    mutate (gene_averages = rowMeans(log_data)) %>% 
    dplyr::filter(gene_averages != "-Inf")
  
  # the last columns is the pseudo-reference column 
  pseudo_column = ncol(log_data)
  
  # where to stop before the pseudo column 
  before_pseduo = pseudo_column - 1
  
  # find the ratio of the log data to the pseudo-reference
  ratios = sweep(log_data[,2:before_pseduo], 1, log_data[,pseudo_column], "-")
  
  # find the median of the ratios
  sample_medians = apply(ratios, 2, median)
  
  # convert the median to a scaling factor
  scaling_factors = exp(sample_medians)
  
  # use scaling factors to scale the original data
  manually_normalized = sweep(data, 2, scaling_factors, "/")
  return(manually_normalized)
}
