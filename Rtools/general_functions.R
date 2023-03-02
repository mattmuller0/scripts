# HEADER --------------------------------------------------
# Author: Matthew Muller
# 
# Date: 2023-01-11
# 
# Script Name: General Functions
# 
# Notes:
# Lots of fun functions I've made to work with.


# LOAD LIBRARIES ------------------------------------------
packages <- c(
  "tidyverse",
  "ggplot2",
  "NMF"
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

rank_estimator <- function(
    expression_matrix,
    ranks = c(2:8),
    save_file = NULL
){
  require(NMF)
  # Runk NMF ranks
  rank_out = nmf(expression_matrix, rank = ranks,
                 method = nmf.getOption("default.algorithm"),
                 nrun = 50)
  
  # Graph cophenetic correlation plots
  cophenetic_correlation_plot <- ggplot(rank_out, 
                                        aes(x=rank, y=cophenetic)) + 
    geom_point() + labs(x="Rank", y="Cophenetic Correlation Coefficient")
  # Save data
  if (!is.null(save_file)) {
    save(rank_out, cophentic_correlation_plot,
         file = save_file)
  }
  output <- list(rank_out, cophenetic_correlation_plot)
  return(output)
}




add_missing_rows <- function(
    df, # cols = samples, rows = genes
    rows, # cols = samples, rows = genes
    sorted = TRUE ){
  missingRowNames <-  rows[which(!rows %in% rownames(df))]
  print(missingRowNames)
  df_tmp <- as.data.frame(matrix(0,
                                 nrow = length(missingRowNames),
                                 ncol = dim(df)[2]
  )
  )
  # print(dim(df_tmp))
  # print(length(missingRowNames))
  colnames(df_tmp) <- colnames(df)
  rownames(df_tmp) <- missingRowNames
  # print(head(df_tmp))
  # print(head(df))
  df <- rbind(df_tmp, df)
  if (sorted) {
    df <- df[order(rownames(df)),]
  }
  
  return(df)
}

