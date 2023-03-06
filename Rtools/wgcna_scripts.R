###########################################################################
#
#                                 HEADER
#
###########################################################################
# Author: Matthew Muller
# 
# Date: 2023-02-06
# 
# Script Name: This will put together the WGCNA scripts below to run them all
# 
# Notes:
# Simple sourcing script to run wgcna. Eventually have inputs to this script,
# but getting lazy for now so each section is unique.


###########################################################################
#
#                                 LIBRARIES
#
###########################################################################
packages <- c(
  "tidyverse",
  "ggplot2",
  "BiocManager",
  "SummarizedExperiment",
  "WGCNA",
  "clusterProfiler",
  "ggbiplot"
)

for (pkg in packages) {
  paste0(pkg)[[1]]
  library(pkg, character.only = T, quietly = T)
}

# LOAD FUNCTIONS
# space reserved for sourcing in functions


###########################################################################
#
#                             Helper Functions
#
###########################################################################

plot_wgcna_dendogram <- function(net){
  # open a graphics window
  # sizeGrWindow(12, 9)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
}

module_wilcoxan_test <- function(MEs, trait){
  moduleCohort_wilcox <- tibble()
  for (ME in names(MEs)){
    MEs_wilcox <- MEs
    MEs_wilcox$cohort <- trait
    formula = as.formula( paste(ME, "cohort", sep="~") )
    
    tmp <- MEs_wilcox %>%
      wilcox_test(data =., formula=formula ) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj")
    moduleCohort_wilcox <- rbind(moduleCohort_wilcox, tmp)
  }
  return(moduleCohort_wilcox)
}


###########################################################################
#
#                           Larger Functions
#
###########################################################################
# Add code here
#
# outdir <- file.path(outdir, paste0('wgcna', '__',format(Sys.time(),"%F"), '/01_initial_process/'))


wgcna_initial_process <- function(
    se, labels, 
    
    outdir,
    
    log2_read_filter = 5
    ){
  require(WGCNA)
  require(DESeq2)
  require(SummarizedExperiment)
  require(ggbiplot)
  require(gridExtra)
  # Output directory
  outdir <- file.path(outdir, paste0('01_initial_process/'))
  dir.create(outdir, showWarnings = F)
  
  #=====================================================================================
  #
  #  Data Prep
  #
  #=====================================================================================
  
  ## Begin WGCNA module creation process using all data samples 

  datExpr0 = t(assay(se))
  
  traitData = colData(se)
  
  #=====================================================================================
  #
  #  goodSamplesGenes Filtering
  #
  #=====================================================================================
  
  
  gsg = goodSamplesGenes(datExpr0, verbose = 2)
  
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }

  #=====================================================================================
  #
  #  Plot sample dendrogram clustering
  #
  #=====================================================================================
  
  sampleTree = hclust(dist(datExpr0), method = "average");
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  pdf(paste0(outdir, "sample_outliers.pdf"))
  #png(file = "Plots/sampleClustering.png", width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
       cex.axis = 1.5, cex.main = 2)
  
  # Plot a line to show the cut
  # abline(h = 8e6, col = "red");
  dev.off()
  # Determine cluster under the line
  clust = cutreeStatic(sampleTree, cutHeight = 8e6, minSize = 10)
  table(clust)
  # clust 1 contains the samples we want to keep.
  keepSamples = (clust==1)

  datExpr = datExpr0[keepSamples, ]
  datTraits <- traitData[keepSamples,]
  
  #=====================================================================================
  #
  #  Replot Agglomerative Clustering
  #
  #=====================================================================================
  pdf(paste0(outdir, "sample_dendrogram.pdf"))
  
  # Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = labels2colors(datTraits[labels])
  # Plot the sample dendrogram and the colors underneath.
  plotDendroAndColors(sampleTree2, traitColors, rowText = datTraits[labels],
                      groupLabels = names(datTraits), addTextGuide = TRUE,
                      main = "Sample dendrogram and trait heatmap",
                      cex.dendroLabels = 0.5, cex.rowText = 0.5)
  
  dev.off()
  #=====================================================================================
  #
  #  More Counts Filtering
  #
  #=====================================================================================

  #Now I move to normalize the seq data 
  #First I will examine the table for total amount of reads across samples 
  
  readTotals = as.data.frame(rowSums(datExpr))
  CountsPerSample<-ggplot(readTotals, aes(x=rowSums(datExpr))) + 
    geom_histogram(color="black", fill="lightblue")+
    ggtitle("Counts for total reads per Sample prefiltering")+
    xlab("Number of reads in Sample")+
    ylab("Count")
  CountsPerSample
  #It does not appear that there are any that need to be removed for low number reads
  
  
  #Next we look at the -log2 of the average reads for each gene. 
  #After identifying where we see a bimodal distribution (if there is one) we will delete the lower levels
  readAvg = as.data.frame(log2(colMeans(datExpr)))
  colnames(readAvg) = c("Average Reads")
  CountsPerDepth_prefilter<-ggplot(readAvg, aes(x=`Average Reads`))+
    geom_histogram(color="black", fill="pink")+
    ggtitle("Counts for rowMeans per Gene postfiltering")+
    labs(x="log2(RowMean)", y="Counts")+
    xlim(-8,10)+
    geom_vline(xintercept=log2_read_filter, color = "black")
  CountsPerDepth_prefilter

  #Based on the result of the graph above, We will cut off any that do not have greater than exp(5) reads in 1/2 the samples
  datExpr = datExpr[,-which(readAvg < log2_read_filter)]
  print(paste0('A total of ', length(which(readAvg < log2_read_filter)), ' have library sizes lower than the cutoff.'))
  
  #Re-examining the histogram to view the new distribution of genes. 
  readAvg = as.data.frame(log2(colMeans(datExpr)))
  colnames(readAvg) = c("Average Reads")
  CountsPerDepth_postfilter<-ggplot(readAvg, aes(x=`Average Reads`))+
    geom_histogram(color="black", fill="pink")+
    ggtitle("Counts for rowMeans per Gene")+
    labs(x="log2(RowMean)", y="Counts")+
    xlim(-8,10)+
    geom_vline(xintercept=5, color = "black")
  CountsPerDepth_postfilter
  
  g <- arrangeGrob(CountsPerDepth_prefilter, CountsPerDepth_postfilter)
  ggsave(paste0(outdir, "count_filtering.pdf"), g)

  #=====================================================================================
  #
  #  Initial DESeq2 Prep and PCA
  #
  #=====================================================================================
  
  se_filtered <- SummarizedExperiment(assays = list(counts = t(datExpr)), colData = datTraits)

  dds <- DESeqDataSet(se_filtered, design = as.formula(paste("~", labels)))
  dds <- estimateSizeFactors(dds)
  vsd <- vst(dds)
  
  datExprNorm <- counts(dds, normalized = T)
  datExprNorm <- assay(vsd)
  
  #Let's examine the samples remaining with a PCA to see if the cluster on the right hand side of the above graph cluster there too
  sample.pca <- prcomp(t(datExprNorm), center = T, scale = T)
  
  pca.plot <- ggbiplot(sample.pca, groups = colData(dds)[[labels]],
                       var.axes = F,ellipse = T) + theme_classic()
  ggsave(paste0(outdir, "pca_plot.pdf"), pca.plot)
  
  ##saving files and work for next analysis point
  save(dds, file=paste0(outdir, "initial_processing_dds.rdata"))
}


wgcna_running <- function(
    se, labels, 
    
    outdir,
    
    sft_RsquaredCut = 0.8
    ){
  # Output directory
  outdir <- file.path(outdir, paste0('02_wgcna_running/'))
  dir.create(outdir, showWarnings = F)
  
  ## Wrap everything in a pdf output
  #=====================================================================================
  #
  #  Data Prep
  #
  #=====================================================================================
  
  ## Begin WGCNA module creation process using all data samples 
  
  datExpr <- t(assay(se))
  datTraits <- colData(se)
  
  
  #=====================================================================================
  #
  #  Thresholding Determination
  #
  #=====================================================================================
  
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = sft_RsquaredCut)
  
  plot_sft_threshold <- function(sft){
    # Plot the results:
    # sizeGrWindow(9, 5)
    par(mfrow = c(1,2));
    cex1 = 0.85;
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red");
    # this line corresponds to using an R^2 cut-off of h
    abline(h=sft_RsquaredCut,col="red")
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  }
  
  pdf(paste0(outdir, "sft_threshold_plots.pdf"))
  plot_sft_threshold(sft)
  dev.off()
  
  
  #=====================================================================================
  #
  #  WGCNA Blockwise Module Formation
  #
  #=====================================================================================
  
  dir.create(paste0(outdir,"/TOMs"), showWarnings = F)
  net = blockwiseModules(datExpr, power = sft[[1]],
                         TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = paste0(outdir, "TOMs/harp_TOM"),
                         verbose = 3)
  
  
  #=====================================================================================
  #
  #  Plot Dendrogram of WGCNA
  #
  #=====================================================================================
  
  pdf(paste0(outdir, "wgcna_dendrogram.pdf"))
  plot_wgcna_dendogram(net)
  dev.off()
  
  #=====================================================================================
  #
  #  Save Modules
  #
  #=====================================================================================
  
  
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  
  save(MEs, moduleLabels, moduleColors, geneTree, 
       file = paste0(outdir, "wgcna_modules.rdata"))
}


wgcna_traits <- function(
    se, labels, 
    
    moduleColors,
    
    outdir
    ){
  require(WGCNA)
  require(SummarizedExperiment)
  require(rstatix)
  require(org.Hs.eg.db)
  require(AnnotationDbi)
  # Output directory
  outdir <- file.path(outdir, paste0('03_wgcna_traits/'))
  dir.create(outdir, showWarnings = F)
  
  ## TEMPLATE LABELING
  datExpr = t(assay(se))
  datTraits = colData(se)
  
  
  #=====================================================================================
  #
  #  Determine module correlations to significant traits
  #
  #=====================================================================================

  # Define numbers of genes and samples
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  
  
  ## Custom Code to check the significance in cohort
  module_wilcoxans <- module_wilcoxan_test(MEs, datTraits[[labels]])
  
  #=====================================================================================
  #
  #  Clinical Traits Matrix
  #
  #=====================================================================================

  pdf(paste0(outdir, "module_trait_correlations.pdf"), 12,8)
  # Will display correlations and their p-values
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5, cex.lab = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  
  
  
  #=====================================================================================
  #
  #  Check Module Importance versus a clinical feature
  #
  #=====================================================================================
  
  
  # Define variable weight containing the weight column of datTrait
  weight = datTraits[[labels]];
  weight = mapvalues(weight, from=unique(datTraits[[labels]]),
                     to=1:length(unique(datTraits[[labels]])) )
  names(weight) = "label"
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  
  ## These two lines need to be fixed
  # names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
  # names(GSPvalue) = paste("p.GS.", names(weight), sep="");
  
  
  #=====================================================================================
  #
  #  Look at Modules versus clinical feature chosen
  #
  #=====================================================================================
  
  
  module = modNames[2]
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  
  pdf(paste0(outdir, "module_membership_significance.pdf"), 7,7)
  # sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for Cohort",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  
  #=============================================================================
  
  
  gene_list <- str_replace(colnames(datExpr), pattern = ".[0-9]+$", replacement = "")
  symbols <- mapIds(org.Hs.eg.db, keys = gene_list,
                    column = c('ENTREZID'), keytype = 'ENSEMBL')
  
  
  #=====================================================================================
  #
  #  Code chunk 9
  #
  #=====================================================================================
  
  
  # Create the starting data frame
  geneInfo0 = data.frame(substanceBXH = gene_list,
                         geneSymbol = symbols,
                         moduleColor = moduleColors,
                         geneSig  = geneTraitSignificance$V1,
                         GVal = GSPvalue$V1)
  # Order modules by their significance for weight
  modOrder = order(-abs(cor(MEs, weight, use = "p")));
  # Add module membership information in the chosen order
  for (mod in 1:ncol(geneModuleMembership)){
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                           MMPvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }
  # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GVal));
  geneInfo = geneInfo0[geneOrder, ]
  
  
  #=====================================================================================
  #
  #  Code chunk 10
  #
  #=====================================================================================
  
  
  write.csv(geneInfo, file = paste0(outdir, "wgcna_geneInfo.csv"))
  write.csv(module_wilcoxans, file = paste0(outdir, "module_wilcoxan_significance.csv"))
}


wgcna_enrichment <- function(
    se,
    
    moduleColors,
    
    outdir
){
  require(WGCNA)
  require(SummarizedExperiment)
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(AnnotationDbi)
  # Output directory
  outdir <- file.path(outdir, paste0('04_wgcna_enrichment/'))
  dir.create(outdir, showWarnings = F)
  
  ## TEMPLATE LABELING
  datExpr = t(assay(se))
  datTraits = colData(se)
  
  
  #=====================================================================================
  #
  #  Determine module correlations to significant traits
  #
  #=====================================================================================
  
  ## TEMPLATE LABELING
  datExpr = t(assay(se_harp_wb))
  datTraits = colData(se_harp_wb)
  
  #=====================================================================================
  #
  #  GO Enrichment analysis
  #
  #=====================================================================================
  
  gene_list <- str_replace(colnames(datExpr), pattern = ".[0-9]+$", replacement = "")
  symbols <- mapIds(org.Hs.eg.db, keys = gene_list,
                    column = c('ENTREZID'), keytype = 'ENSEMBL')
  
  ###########################################################################
  #
  #                               Plot Enrichment data
  #
  ###########################################################################
  moduleGenes <- lapply(unique(moduleColors), 
                        function(module) {symbols[which(moduleColors == module)]
  })
  
  # perform gene set enrichment analysis for each module
  for (i in 1:length(moduleGenes)) {
    module <- na.omit(moduleGenes[[i]])
    result <- enrichGO(module, org.Hs.eg.db, keyType = "ENTREZID")
    if (length(result$ID) > 0) {
      write.table(result, paste0(outdir, "module_", i, "_enrichment_results.txt"), 
                  sep = "\t", quote = FALSE, row.names = FALSE)
      dotplot(result)
      ggsave(paste0(outdir, "module_", i, "_enrichment_dotplot.pdf"))
    }
  }
  
  #=====================================================================================
  #
  #  Save Data
  #
  #=====================================================================================
  
  # write.table(tab, file = "analysis/wgcna/temp/harp_wb_GOEnrichmentTable.csv", 
  #             sep = ",", quote = TRUE, row.names = FALSE)
}


wgcna_visualization <- function(
    se,
    
    moduleColors,
    
    outdir
){
  require(WGCNA)
  require(SummarizedExperiment)
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(AnnotationDbi)
  # Output directory
  outdir <- file.path(outdir, paste0('05_wgcna_visualization/'))
  dir.create(outdir, showWarnings = F)
  
  ## TEMPLATE LABELING
  datExpr = t(assay(se))
  datTraits = colData(se)
  
  
  #=====================================================================================
  #
  #  Dissimilarity Plot
  #
  #=====================================================================================
  
  # Calculate topological overlap anew: this could be done more efficiently by saving the TOM
  # calculated during module detection, but let us do it again here.
  dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
  
  nSelect = 2000 # We use this just to make things faster
  # For reproducibility, we set the random seed
  set.seed(10);
  select = sample(dim(datExpr)[2], size = nSelect);
  selectTOM = dissTOM[select, select];
  # There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select];
  # Open a graphical window
  pdf(paste0(outdir, "TOMplot.pdf"), 9,9)
  # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
  # the color palette; setting the diagonal to NA also improves the clarity of the plot
  plotDiss = selectTOM^10;
  diag(plotDiss) = NA;
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
  dev.off()
  
  
  #=====================================================================================
  #
  #  Plot Eigengene adjacency / similarity
  #
  #=====================================================================================
  
  # Recalculate module eigengenes
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MET = orderMEs(MEs)
  # Plot the relationships among the eigengenes and the trait
  pdf(paste0(outdir, "eigenegene_adjacency_plot.pdf"), 5,7.5)
  par(cex = 0.6)
  plotEigengeneNetworks(MET, "Eigengene Adjacency Plot", marDendro = c(0,4,1,2), 
                        marHeatmap = c(3,4,1,2), xLabelsAngle= 90)
  dev.off()
}