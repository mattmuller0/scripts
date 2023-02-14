# HEADER --------------------------------------------------
# Author: Matthew Muller
# 
# Date: 2023-01-22
# 
# Script Name: Gene Ratio Functions
# 
# Notes:
# So this is all the gene ratio stuff I've been working with.
# Now I'm trying to just bring them all to one spot

# SET WORKING DIRECTORY -----------------------------------
wd <- FALSE
if (wd != FALSE) {
  setwd(wd)
  cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
  }

# LOAD LIBRARIES ------------------------------------------
packages <- c(
  "tidyverse",
  "ggplot2",
  "skimr",
  "DESeq2",
  "reshape2",
  "scales",
  "corrplot",
  "umap"
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

###############################################################
#
#             Create some things to use first
#
###############################################################

## Load in DGE data from duke and pace datasets
pace_dge <- read.csv('data/hyper_v_hypo_deseqoutput.csv')
duke_dge <- read.csv('data/duke_validation_run3/dge_analysis/comp_group1__hyper_v_nothyper_AGCONTROL_deseqout.csv')
press_genes <- t(read.csv('data/clean/genes.csv'))

# Select press genes present in both
pace_dge <- pace_dge[pace_dge$X %in% press_genes,]
pace_dge <- pace_dge[order(pace_dge$X),]
duke_dge <- duke_dge[duke_dge$X %in% press_genes,]
duke_dge <- duke_dge[order(duke_dge$X),]

dge_data <- subset(pace_dge, select=c("X", "log2FoldChange"))
colnames(dge_data) <- c("Genes", "Pace")
dge_data$"Duke" <- duke_dge$log2FoldChange


# Filter with p value
dge_filtered_sign <- dge_data[duke_dge$pvalue < 0.05,]

# Take only genes that align and write them to csvs
gene_list_up_sign <- dge_filtered_sign[dge_filtered_sign$Pace > 0 & dge_filtered_sign$Duke > 0,]
gene_list_down_sign <- dge_filtered_sign[dge_filtered_sign$Pace < 0 & dge_filtered_sign$Duke < 0,]

## Load in the correlated gene list data (data is median of ratios)
subset_up <- gene_list_up_sign$Genes
subset_down <- gene_list_down_sign$Genes

pace_genes <- read.csv('data/clean/pace/features.csv')
pace_labels <- read.csv('data/clean/pace/labels.csv')$X0

duke_genes <- read.csv('data/clean/duke/features_group1.csv')
duke_labels <- read.csv('data/clean/duke/labels_group1.csv')$compouttable...2.

subset_up <- subset_up[subset_up %in% colnames(pace_genes)]
subset_down <- subset_down[subset_down %in% colnames(pace_genes)]


###############################################################
#
#       First Set of Functions for making gene ratios
#
###############################################################





gene_ratio_iterator <- function(
    up_samples_df, down_samples_df){
  # This function brute forces the ratio of every gene in up_df (up-regulated genes) 
  # with those in down_df (down-regulated genes). It is assumed the dataframes are
  # built identically and have the same number of rows (ie. samples).
  
  # Initialize things...
  up_genes <- colnames(up_samples_df)
  down_genes <- colnames(down_samples_df)
  gene_ratio_df <- data.frame(matrix(ncol = length(down_genes)*length(up_genes),
                                     nrow = dim(up_samples_df)[1]
  )
  )
  
  # Fill in ratios
  ratio_names <- c()
  for (sample in 1:dim(up_samples_df)[1]) {
    for (down_idx in 1:length(down_genes)) {
      for (up_idx in 1:length(up_genes)) {
        # Get values we need
        up_gene <- up_genes[up_idx]
        down_gene <- down_genes[down_idx]
        if (sample == 1) {
          ratio_names <- append(ratio_names, paste0(up_gene,"/",down_gene))
        }
        
        # Add in the ratio value
        idx <- up_idx+(down_idx-1)*length(up_genes)
        gene_ratio_df[sample,idx] <- up_samples_df[sample,up_idx] / down_samples_df[sample,down_idx]
      }
    }
  }
  colnames(gene_ratio_df) <- ratio_names
  return(gene_ratio_df)
}

gene_to_ratio_table_converter <- function(
    genes, labels){
  # Subset the gene data and summarize it with skim
  pace_up_genes_hyper <- genes[!!labels,] %>%
    select(subset_up)
  
  pace_down_genes_hyper <- genes[!!labels,] %>%
    select(subset_down)
  
  # Call the gene ratio iterator function
  pace_hyper_ratios <- gene_ratio_iterator(pace_up_genes_hyper, pace_down_genes_hyper)
  is.na(pace_hyper_ratios)<-sapply(pace_hyper_ratios, is.infinite)
  pace_hyper_ratios[is.na(pace_hyper_ratios)]<-0
  
  
  
  
  # Subset the gene data and summarize it with skim
  pace_up_genes_normal <- genes[!!!labels,] %>%
    select(subset_up)
  
  pace_down_genes_normal <- genes[!!!labels,] %>%
    select(subset_down)
  
  # Call the gene ratio iterator function
  pace_normal_ratios <- gene_ratio_iterator(pace_up_genes_normal, pace_down_genes_normal)
  is.na(pace_normal_ratios)<-sapply(pace_normal_ratios, is.infinite)
  pace_normal_ratios[is.na(pace_normal_ratios)] <- 0
  
  ret <- list(pace_hyper_ratios, pace_normal_ratios)
  return(ret)
}

ratio_significance_finder <- function(genes, labels){
  # Subset the gene data and summarize it with skim
  pace_up_genes_hyper <- genes[!!labels,] %>%
    select(subset_up)
  
  pace_down_genes_hyper <- genes[!!labels,] %>%
    select(subset_down)
  
  # Call the gene ratio iterator function
  pace_hyper_ratios <- gene_ratio_iterator(pace_up_genes_hyper, pace_down_genes_hyper)
  is.na(pace_hyper_ratios)<-sapply(pace_hyper_ratios, is.infinite)
  pace_hyper_ratios[is.na(pace_hyper_ratios)]<-0
  
  
  
  
  # Subset the gene data and summarize it with skim
  pace_up_genes_normal <- genes[!!!labels,] %>%
    select(subset_up)
  
  pace_down_genes_normal <- genes[!!!labels,] %>%
    select(subset_down)
  
  # Call the gene ratio iterator function
  pace_normal_ratios <- gene_ratio_iterator(pace_up_genes_normal, pace_down_genes_normal)
  is.na(pace_normal_ratios)<-sapply(pace_normal_ratios, is.infinite)
  pace_normal_ratios[is.na(pace_normal_ratios)] <- 0
  
  pace_hyper_ratios_skim <- pace_hyper_ratios %>%skim()
  pace_normal_ratios_skim <- pace_normal_ratios %>%skim()
  # hist(pace_hyper_ratios_skim$numeric.mean - pace_normal_ratios_skim$numeric.mean, breaks="fd")
  
  pace_ratios = data.frame(pace_hyper_ratios_skim$skim_variable,
                           pace_hyper_ratios_skim$numeric.mean,
                           pace_hyper_ratios_skim$numeric.sd,
                           pace_normal_ratios_skim$numeric.mean,
                           pace_hyper_ratios_skim$numeric.sd)
  colnames(pace_ratios) <- c("Gene Ratio", "Hyper Mean", "Hyper SD", "Normal Mean", "Normal SD")
  pace_ratios$Difference <- pace_hyper_ratios_skim$numeric.mean - pace_normal_ratios_skim$numeric.mean
  
  # pace_ratios <- pace_ratios[pace_ratios$Difference > 500,]
  
  significance_testing <- list()
  for (gene in colnames(pace_hyper_ratios)) {
    tmp <- t.test(pace_hyper_ratios[gene],pace_normal_ratios[gene])
    significance_testing <- append(significance_testing, tmp$p.value)
  }
  
  pace_ratios$pvalue <- significance_testing
  return(pace_ratios)
}



###############################################################
#
#      Now some prediction and metric tools
#
###############################################################




gene_ratio_predictor <- function(
    # Input dataframes need to have more than one column right now!!!
  gene_ratio_df_condition,
  gene_ratio_df_normal,
  thresholds = NULL, # should be a vector
  method = median # a function to calculate the thresholds
  # algorithm will compute midpoint of the method on the two conditions
){
  if (is.null(thresholds)) {
    print("Setting Thresholds")
    median_condition <- sapply(gene_ratio_df_condition, method)
    median_normal <- sapply(gene_ratio_df_normal, method)
    thresholds <- (median_condition + median_normal)/2
  }
  thresholds_cond <- matrix(nrow=dim(gene_ratio_df_condition)[1], 
                            ncol=dim(gene_ratio_df_condition)[2] )
  
  thresholds_norm <- matrix(nrow=dim(gene_ratio_df_normal)[1], 
                            ncol=dim(gene_ratio_df_normal)[2] )
  for (idx in 1:length(thresholds)) {
    thresholds_cond[,idx] <- thresholds[idx]
    thresholds_norm[,idx] <- thresholds[idx]
  }
  
  true_condition <- colSums(gene_ratio_df_condition >= thresholds_cond)
  total_condition <- dim(gene_ratio_df_condition)[1]
  sensitivity <- true_condition / total_condition
  
  true_normal <- colSums(gene_ratio_df_normal <= thresholds_norm)
  total_normal <- dim(gene_ratio_df_normal)[1]
  specificity <- true_normal / total_normal
  
  true_total <- true_normal + true_condition
  total <- total_normal + total_condition
  accuracy <- true_total / total
  
  results <- t(data.frame(thresholds, sensitivity, specificity, accuracy))
  
  return(results)
}

gene_ratio_sample_scorer <- function(
    # Input dataframes need to have more than one column right now!!!
  gene_ratio_df_condition,
  gene_ratio_df_normal,
  thresholds = NULL, # should be a vector
  method = median # a function to calculate the thresholds
  # algorithm will compute midpoint of the method on the two conditions
){
  if (is.null(thresholds)) {
    print("Setting Thresholds")
    median_condition <- sapply(gene_ratio_df_condition, method)
    median_normal <- sapply(gene_ratio_df_normal, method)
    thresholds <- (median_condition + median_normal)/2
  }
  thresholds_cond <- matrix(nrow=dim(gene_ratio_df_condition)[1], 
                            ncol=dim(gene_ratio_df_condition)[2] )
  
  thresholds_norm <- matrix(nrow=dim(gene_ratio_df_normal)[1], 
                            ncol=dim(gene_ratio_df_normal)[2] )
  for (idx in 1:length(thresholds)) {
    thresholds_cond[,idx] <- thresholds[idx]
    thresholds_norm[,idx] <- thresholds[idx]
  }
  sample_scores_cond <- data.frame(scores = rowSums(gene_ratio_df_condition >= thresholds_cond))
  sample_scores_cond$label <- "positive"
  
  sample_scores_norm <- data.frame(scores = rowSums(gene_ratio_df_normal >= thresholds_norm))
  sample_scores_norm$label <- "normal"
  
  results <- rbind(sample_scores_cond, sample_scores_norm)
  
  return(results)
}

gene_ratio_majority_voting <- function(
    # Input dataframes need to have more than one column right now!!!
  gene_ratio_df,
  labels, # labels should be binary vector
  weights = NULL, # use if weighing the gene ratios
  
  thresholds = NULL, # should be a vector
  margin = NULL, # error of thresholds used to determine prediction confidence
  C = 1,
  
  method = median, # a function to calculate the thresholds
  probability = F
  # algorithm will compute midpoint of the method on the two conditions
){
  
  # Subset the data into the condition and normal groups
  label_ <- unique(labels)
  
  gene_ratio_df_condition <- gene_ratio_df[labels == label_[1],]
  gene_ratio_df_normal <- gene_ratio_df[labels == label_[2],]
  
  if (is.null(weights)) {
    
    if (!probability) {
      if (is.null(thresholds)) {
        print("Setting Thresholds")
        median_condition <- sapply(gene_ratio_df_condition, method)
        median_normal <- sapply(gene_ratio_df_normal, method)
        thresholds <- (median_condition + median_normal) / 2
      }
      
      thresholds_ <- matrix(nrow=dim(gene_ratio_df)[1], 
                            ncol=dim(gene_ratio_df)[2] )
      
      for (idx in 1:length(thresholds)) {thresholds_[,idx] <- thresholds[idx]}
      
      sample_scores <- data.frame(scores = rowSums(gene_ratio_df >= thresholds_))
      sample_scores <- sample_scores >= dim(gene_ratio_df)[2]/2
      sample_scores <- data.frame(
        scores = mapvalues(sample_scores, from=c(T, F), to=c(label_[1], label_[2]))
      )
      return(sample_scores)
    }
    
    if (probability) {
      print("Using probability function. Must have margins.")
      if (is.null(thresholds)) {
        print("Setting thresholds and margins")
        median_condition <- sapply(gene_ratio_df_condition, method)
        sd_condition <- sapply(gene_ratio_df_condition, function(x){sd(x)/sqrt(length(x))})
        
        median_normal <- sapply(gene_ratio_df_normal, method)
        sd_normal <- sapply(gene_ratio_df_normal, function(x){sd(x)/sqrt(length(x))})
        
        thresholds <- (median_condition + median_normal)
        margin <- (sd_condition + sd_normal) * C # following error propagation logic
      }
      
      thresholds_ <- matrix(nrow=dim(gene_ratio_df)[1], 
                            ncol=dim(gene_ratio_df)[2] )
      for (idx in 1:length(thresholds)) {thresholds_[,idx] = thresholds[idx]}
      
      sample_scores <- (gene_ratio_df - thresholds) / margin
      sample_scores[ sample_scores > 1 ] <- 1
      sample_scores[ sample_scores < -1] <- -1
      sample_scores <- rowSums(sample_scores) > dim(gene_ratio_df)[2]/2
      sample_scores <- data.frame(
        scores = mapvalues( sample_scores, from = c(T, F), to = c(label_[1], label_[2]) )
      )
      return(sample_scores)
    }
  }
  
  if (!is.null(weights)) {
    
    if (!probability) {
      if (is.null(thresholds)) {
        print("Setting Thresholds")
        median_condition <- sapply(gene_ratio_df_condition, method)
        median_normal <- sapply(gene_ratio_df_normal, method)
        thresholds <- (median_condition + median_normal) / 2
      }
      
      thresholds_ <- matrix(nrow=dim(gene_ratio_df)[1], 
                            ncol=dim(gene_ratio_df)[2] )
      weights_ <- matrix(nrow=dim(gene_ratio_df)[1], 
                         ncol=dim(gene_ratio_df)[2] )
      for (idx in 1:length(thresholds)) {thresholds_[,idx] <- thresholds[idx]}
      for (idx in 1:length(weights)) {weights_[,idx] = weights[idx]}
      
      sample_scores <- as.numeric(gene_ratio_df >= thresholds_) * weights_
      sample_scores <- data.frame(scores = rowSums(sample_scores))
      sample_scores <- sample_scores >= dim(gene_ratio_df)[2]/2
      sample_scores <- data.frame( 
        scores = mapvalues(sample_scores, from=c(T, F), to=c(label_[1], label_[2]))
      )
      return(sample_scores)
    }
    
    if (probability) {
      print("Using probability function. Must have margins.")
      if (is.null(thresholds)) {
        print("Setting thresholds and margins")
        median_condition <- sapply(gene_ratio_df_condition, method)
        sd_condition <- sapply(gene_ratio_df_condition, function(x){sd(x)/length(x)})
        
        median_normal <- sapply(gene_ratio_df_normal, method)
        sd_normal <- sapply(gene_ratio_df_normal, function(x){sd(x)/length(x)})
        
        thresholds <- (median_condition + median_normal)
        margin <- (sd_condition + sd_normal) * C # following error propagation logic
      }
      
      thresholds_ <- matrix(nrow=dim(gene_ratio_df)[1], 
                            ncol=dim(gene_ratio_df)[2] )
      weights_ <- matrix(nrow=dim(gene_ratio_df)[1], 
                         ncol=dim(gene_ratio_df)[2] )
      for (idx in 1:length(thresholds)) {thresholds_[,idx] = thresholds[idx]}
      for (idx in 1:length(weights)) {weights_[,idx] = weights[idx]}
      
      sample_scores <- (gene_ratio_df - thresholds) / margin
      sample_scores[ sample_scores > 1 ] <- 1
      sample_scores[ sample_scores < -1] <- -1
      sample_scores <- sample_scores * weights_
      sample_scores <- rowSums(sample_scores) > dim(gene_ratio_df)[2]/2
      sample_scores <- data.frame(
        scores = mapvalues( sample_scores, from = c(T, F), to = c(label_[1], label_[2]) )
      )
      return(sample_scores)
    }
  }
  
}

gene_ratio_threshold <- function(gene_ratio_df, labels, method = median){
  label_ <- unique(labels)
  
  gene_ratio_df_condition <- gene_ratio_df[labels == label_[1],]
  gene_ratio_df_normal <- gene_ratio_df[labels == label_[2],]
  
  print("Setting thresholds")
  median_condition <- sapply(gene_ratio_df_condition, method)
  sd_condition <- sapply(gene_ratio_df_condition, sd)
  
  median_normal <- sapply(gene_ratio_df_normal, method)
  sd_normal <- sapply(gene_ratio_df_normal, sd)
  
  thresholds <- (median_condition + median_normal) / 2
  margin <- (sd_condition + sd_normal) # following error propagation logic
  return(thresholds)
}

gene_ratio_margin <- function(gene_ratio_df, labels, method = median){
  label_ <- unique(labels)
  
  gene_ratio_df_condition <- gene_ratio_df[labels == label_[1],]
  gene_ratio_df_normal <- gene_ratio_df[labels == label_[2],]
  
  print("Setting margins")
  median_condition <- sapply(gene_ratio_df_condition, method)
  sd_condition <- sapply(gene_ratio_df_condition, sd)
  
  median_normal <- sapply(gene_ratio_df_normal, method)
  sd_normal <- sapply(gene_ratio_df_normal, sd)
  
  thresholds <- (median_condition + median_normal) / 2
  margin <- (sd_condition + sd_normal) / 2# following error propagation logic
  return(margin)
}


plot_umap <- function(x, labels, main="A UMAP Visualization") {
  
  layout <- x
  if (is(x, "umap")) {
    layout <- x$layout
  }
  
  umap_df <- data.frame(
    "UMAP.1" = layout[,1], 
    "UMAP.2" = layout[,2],
    "cohort" = labels )
  ggplot(umap_df, aes(x=UMAP.1, y=UMAP.2, col = cohort)) + 
    geom_point() + theme_minimal() + labs(title=main)
}




