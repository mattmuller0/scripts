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




nmf_estimator <- function(
    counts_matr,
    outfile,
    ranks = 2:8, runs = 30,
    options = 'v',
    ){
  require(NMF)
  # Call nmf with some settings to make test the ranks chosen
  # We also are going to randomize the data to see how it compares to
  # data that is purely random in terms of residuals and other metrics
  # Ideally, the next step is to plot the NMF
  nmf_ranks_out <- nmf(counts_matr, ranks, nrun = runs, 
                       .opt = options) 
  nng_ranks_out <- nmf(randomize(counts_matr), ranks, nrun = runs, 
                       .opt = options) 
  save(nmf_ranks_out, nng_ranks_out, file = outfile)
}

  
  
  
  
  
  
  
  
  
  



