################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Read in each file from output of htseq count, and concat together naming each column by the filename, then plot out some data cleaning stuff and normalize the data, then write out those tiles

## Load in Libraries
packagelist = c("ggplot2", "reshape2", "DESeq2", "grid", "gridExtra", "scales", "ggrepel", "tools")
junk <- lapply(packagelist, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

# path to all of the clinical phenocyte csv files
#filepath = "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia/data/quant-htseq_new/"

# Out path for plotting
#outfilepath = "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia/output/rna_processing3/"
#dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)

## Read in the list of files and combine to a raw count table
## Combine the files into a counttable - then remove columns that are complete duplicates of each other, and sort the table by the PATNUM
## Also remove the rows that are not genes!
#####
create_count_table <- function(countfilepath, outfilepath){
    filelist1 = paste(countfilepath, list.files(path = countfilepath,  pattern = "*.txt"), sep="")
    # Filter out the files that have no information (size = 0)
    filelist = filelist1[!file.size(filelist1) == 0]
      
    # Read in an combine the count files
    tablist = list()
    for (filenum in 1:length(filelist)) {
        if (!file.size(filelist[filenum]) == 0) {
            temp = read.table(filelist[filenum], sep="\t", row.names = 1, stringsAsFactors = FALSE)
            colnames(temp) = sub(pattern = "\\..*","", basename(filelist[filenum]))
            tablist[[filenum]] = temp[,1,drop=FALSE]
        }
    }
    counttab = do.call(cbind, tablist)

    # counttab = counttab[,order(gsub(pattern = ".*_","",colnames(counttab)))]
    write.table(counttab, file = paste(outfilepath, "concat_count_files.txt", sep=""), sep="\t", col.names=NA, row.names=TRUE,quote = FALSE)
  
    rownameomitlist = c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
    counttab = counttab[!rownames(counttab) %in% rownameomitlist,]
  
    return(counttab)
}


#### SELECT FOR COLUMNS THAT FALL INTO OUR SPECIFIED LIST OF FILES
# customsamplelist = read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia/data/new_data_sheets/out/sample_list_DATE_PATNUM_custom.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)[,1]
# 
# write.table(counttab, file = paste(outfilepath, "concat_count_files_rmdup.txt", sep=""), sep="\t", col.names=NA, row.names=TRUE,quote=FALSE)

#####

## Density Curves with each sample labeled - to scroll through and find outliers - DONT WANT TO RUN THIS RIGHT NOW
#####
plot_read_density_curves <- function(counttab, outfilepath) {
    dir.create(paste(outfilepath, "outlier_detection_density_curves/", sep=""), recursive = TRUE, showWarnings = FALSE)

    # NEW STUFF
    ## Determine the average density function
    avg_dens <- density(log2(rowMeans(counttab)))
    # Perform ks test for each curve against the average https://eranraviv.com/test-of-equality-between-two-densities/
    ksoutlist <- list()
    for (kstestnum in seq_len(ncol(counttab))){
        ksout <- suppressWarnings(ks.test(log2(counttab[,kstestnum]), log2(rowMeans(counttab))))
        ksoutlist[[kstestnum]] <- ksout
        names(ksoutlist)[kstestnum] <- colnames(counttab[,kstestnum,drop=FALSE])
    }
    ksouttab <- do.call(rbind, ksoutlist)
    rownames(ksouttab) <- gsub("\\.D","",rownames(ksouttab))
    ksstatscaled <- scale(unlist(ksouttab[,1]))
    rownames(ksstatscaled) <- gsub("\\.D","",rownames(ksstatscaled))
    colnames(ksstatscaled) <- "SD_kstat_diff"
    ## This kind of works again - amazing how this keeps working - check to see which D scores are 3sd from the average
    # densityoutliers <- ksstatscaled[abs(ksstatscaled) > 3,1,drop=FALSE]
    densityoutliers <- ksstatscaled[abs(ksstatscaled) > 2,1,drop=FALSE]

    ## First - plot all of the density curves without emphasizing an outlier
    plotdatatab = melt(log2(counttab), id.vars = NULL)
    
    pout_all <- ggplot(plotdatatab, aes(x=value, group=variable))
    pout_all <- pout_all + geom_density()
    pout_all <- pout_all + labs(title = "Counts per log2(Genecount) for all samples", x = "log2(genecount)", y = "density")
    pout_all <- pout_all + coord_cartesian(ylim = c(0,0.5), xlim=c(0,15))
    
    pdf(paste(outfilepath, "outlier_detection_density_curves/all_samples_outlier_density_curves.pdf", sep=""))
    suppressWarnings(print(pout_all))
    dev.off()
    
    ## Now let's only plot each line if its an outlier? That's what we really care about anyways
    for (outliernum in seq_len(nrow(densityoutliers))) {
        plottab = melt(log2(counttab), id.vars = NULL)
        plotdatatab$outlier = "no"
        plotdatatab[plotdatatab$variable %in% rownames(densityoutliers)[outliernum],"outlier"] <- "yes"
        plotdatatab = plotdatatab[order(match(plotdatatab[,3], c("yes","no"))),]
        
        pout <- ggplot(plotdatatab, aes(x=value, group=variable, color=outlier, size=outlier))
        pout <- pout + geom_density()
        pout <- pout + geom_density(data = subset(plotdatatab, variable == rownames(densityoutliers)[outliernum]), 
                                    aes(x = value, color = outlier, size = outlier))
        pout <- pout + scale_color_manual(limits = c("yes", "no"), values = c("red", "black"))
        pout <- pout + scale_size_manual(limits = c("yes", "no"), values = c(1, 0.2))
        pout <- pout + labs(title = paste0("Counts per log2(Genecount) for all samples with highlighted ", 
                                          rownames(densityoutliers)[outliernum]), x = "log2(genecount)", y = "density")
        pout <- pout + coord_cartesian(ylim = c(0,0.5), xlim=c(0,15))
        
        pdf(paste(outfilepath, "outlier_detection_density_curves/", rownames(densityoutliers)[outliernum], "_outlier_density_curves.pdf", sep=""))
        suppressWarnings(print(pout))
        dev.off()
        
    }

    return(list(density_plot = pout_all, density_outliers = densityoutliers))
    
}

#####

## Define the histogram function that we are using


## Perform filtering by number of reads
#####
QC_filter_readcount <- function(counttab, minreadcutoff) {
    p1datatab = data.frame(sampreadsum = colSums(counttab))
    pout1 <- plot_histogram(data = p1datatab, 
                   fitcurve=FALSE, binparam = 50, limitx = c(floor(min(p1datatab)),ceiling(max(p1datatab))), 
                   labsparam = list(title = "Counts for total reads per Sample prefilter", x = "Number of reads in Sample", y = "Count"))
    pout1 <- pout1 + geom_histogram(color="black", fill="navy", bins = 50)
    
    
    ## Select for less than 2,000,000 reads (empirically determined)
    # NOTE - saving out the readcount to add in as metadata for further QC later on
    counttabfilt = counttab[,colSums(counttab) > minreadcutoff]
    readcounttab = p2datatab = data.frame(sampreadsum = colSums(counttabfilt))
    
    pout2 <- plot_histogram(data = p2datatab, 
                   fitcurve=FALSE, binparam = 50, limitx = c(floor(min(p1datatab)),ceiling(max(p1datatab))), 
                   labsparam = list(title = "Counts for total reads per Sample postfilter", x = "Number of reads in Sample", y = "Count"))
    pout2 <- pout2 + geom_histogram(color="black", fill="navy", bins = 50)
    
    print(dim(counttabfilt))
    omittedsamples = colnames(counttab[,colSums(counttab) <= minreadcutoff,drop=FALSE])
    return(list(prefilterhist = pout1, postfilterhist = pout2, 
                counttabfilt = counttabfilt, outliercounts = data.frame(readcount = colSums(counttab[,omittedsamples,drop=FALSE]))))
}

#####

## Perform filtering by min samples with reads in genes
# mincountcutoff = 4
# minsamplescutoff = round(ncol(counttab)/2)
#####
## Process data for ggplotting
QC_filter_genecount <- function(counttabfilt2, mincountcutoff, minsamplescutoff, outfilepath) {
    p3datatab = data.frame(genemeans = log2(rowMeans(counttabfilt2)))
    p3datatab = p3datatab[is.finite(rowSums(p3datatab)),,drop=FALSE]
    
    pout1 <- plot_histogram(data = p3datatab, 
                   fitcurve=FALSE, binparam = 50, limitx = c(floor(min(p3datatab)),ceiling(max(p3datatab))), 
                   labsparam = list(title = "Counts for rowMeans per Gene prefilter", x = "log2(RowMean)", y = "Count"))
    pout1 <- pout1 + geom_histogram(color="black", fill="navy", bins = 50)
    
    ## Applying some prefiltering - lets say for now we are going to do a minimum value of [4] in at least [50%] of samples
    counttabfilt <- counttabfilt2[rowSums(counttabfilt2 >= mincountcutoff ) >= minsamplescutoff,]
    
    ## Do processing a plotting again
    p4datatab = data.frame(genemeans = log2(rowMeans(counttabfilt)))
    p4datatab = p4datatab[is.finite(rowSums(p4datatab)),,drop=FALSE]
    
    pout2 <- plot_histogram(data = p4datatab, 
                   fitcurve=FALSE, binparam = 50, limitx = c(floor(min(p3datatab)),ceiling(max(p3datatab))), 
                   labsparam = list(title = "Counts for rowMeans per Gene postfilter", x = "log2(RowMean)", y = "Count"))
    pout2 <- pout2 + geom_histogram(color="black", fill="navy", bins = 50)
    
    print(dim(counttabfilt))
    return(list(prefilterhist = pout1, postfilterhist = pout2, counttabfilt = counttabfilt))
}
#####



## Plot the top genes that are proportionately representative in our dataset
numgenesplot = 50
#####
# Get legend function to print on a separate page
# g_legend <- function(a.gplot){ 
#   tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
#   legend <- tmp$grobs[[leg]] 
#   legend
# } 
# Define a color list - extracts 433 colors from base R colors (no grays) and saves it as a list to sample from
plot_gene_prop_histogram <- function(counttabfilt4, numgenesplot) {
    plotcolors = rep(c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00",
                      "black", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
                      "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
                      "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown"),
                    ceiling(numgenesplot/25))
    
    # Turn the table we have into a table of gene proportions
    proptab <- apply(counttabfilt4, MARGIN=2, function(x) x/sum(x, na.rm=TRUE))
    subproptab = proptab[names(sort(rowMeans(counttabfilt4), decreasing = TRUE))[1:numgenesplot],]
    
    # proptaboutlierlist <- list()
    # for (proptabnum in seq_len(nrow(subproptab))){
    #     scaledproptab <- scale(subproptab[proptabnum,])
    #     proptaboutlierlist[[proptabnum]] <- data.frame(rownames(scaledproptab[abs(scaledproptab)>3,1,drop=FALSE]))
    #     names(proptaboutlierlist)[proptabnum] <- rownames(proptab[proptabnum,,drop=FALSE])
    # }
    # proptaboutliertab <- do.call(rbind, proptaboutlierlist)
    
    propplottab = melt(subproptab)
    propplottab[,2] = as.character(propplottab[,2])
    colnames(propplottab) = c("gene", "sample", "value")
    
    g2 <- ggplot()
    g2 <- g2 + geom_bar(aes(y = value, x = sample, fill = gene), data = propplottab, stat="identity", size=0.5, color="black")
    g2 <- g2 + scale_fill_manual(values = plotcolors)
    legend <- g_legend(g2)
    g2 <- g2 + theme(legend.position = "none", axis.text.x = element_text(angle = 90))
    
    # pdf(paste(outfilepath, "read_distribution_bar_chart.pdf", sep=""), height = 15, width = max(10,round(ncol(proptab)/10)))
    # print(g2)
    # grid.newpage()
    # grid.draw(legend) 
    # junk <- dev.off()
    return(list(geneproptab = proptab, geneprop_plot = g2, geneprop_plotlegend = legend))
}
#####


## PCA Plotting
#####

## Define the PCA plotting function
# pca_plotter <- function(pcadata, colorvar, scalecolor=FALSE, labsparam, outfile) {
# pca_plotter <- function(pcadata, colorvar=NULL, shapevar = NULL, scalecolor=FALSE, labelpoints = FALSE, labsparam, outfile) {
#   pcdata = prcomp(pcadata)
#   pc1 = pcdata$x[,1]
#   pc2 = pcdata$x[,2]
#   
#   pcloadings = as.data.frame(pcdata$rotation)
#   pc1genes = rownames(pcloadings[order(pcloadings[,1], decreasing=TRUE),])[1:5]
#   pc2genes = rownames(pcloadings[order(pcloadings[,2], decreasing=TRUE),])[1:5]
#   pcproportions = summary(pcdata)$importance[2,]
#   
#   xlabel = paste(labsparam$x, "_", percent(pcproportions[1]), "_", paste(pc1genes, collapse = "_"), sep="")
#   ylabel = paste(labsparam$y, "_", percent(pcproportions[2]), "_", paste(pc2genes, collapse = "_"), sep="")
#   
#   pcaplotdata = cbind.data.frame(PC1 = pc1, PC2 = pc2)
#   if (!is.null(colorvar)) {pcaplotdata$colorvar = colorvar[,1]}
#   if (!is.null(shapevar)) {pcaplotdata$shapevar = shapevar[,1]}
#   
#   if (labelpoints == TRUE) {pout <- ggplot(data = pcaplotdata, aes(x=pcaplotdata[,1], y = pcaplotdata[,2], label=rownames(pcaplotdata)))} else {pout <- ggplot(data = pcaplotdata, aes(x=pcaplotdata[,1], y = pcaplotdata[,2]))}
#   if (!is.null(colorvar)) {pout <- pout + aes(color=pcaplotdata$colorvar)}
#   if (!is.null(shapevar)) {pout <- pout + aes(shape=pcaplotdata$shapevar)}
#   pout <- pout + geom_point(size = 2.5)
#   if (labelpoints == TRUE) {pout <- pout + geom_text_repel(box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"))}
#   pout <- pout + labs(title = labsparam$title, x = xlabel, y = ylabel, color=labsparam$color, shape=labsparam$shape)
#   if (scalecolor == TRUE) {pout <- pout + scale_color_gradient(low = "yellow", high = "blue")}
#   legend <- g_legend(pout)
#   pout <- pout + theme(legend.position = "none")
#   
#   pdf(file = outfile)
#   print(pout)
#   grid.newpage()
#   grid.draw(legend) 
#   junk <- dev.off()
#   
# }

#####

## DESeq normalization
#####
# print("Checks for DESeq to be run correctly - both should be true")
# all(rownames(metatablefilt3) %in% colnames(counttabfilt4))
# all(rownames(metatablefilt3) == colnames(counttabfilt4))

## Create the DEseq object and extract our normalization factors
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat
# https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
DEseq_Normalization <- function(counttable, metatable, outfilepath, label_extreme_changes = FALSE) {
    dds <- DESeqDataSetFromMatrix(countData = counttable,
                                  colData = metatable,
                                  design = ~ 1 )
    dds <- estimateSizeFactors(dds)
    normcounttab <- counts(dds, normalized=TRUE)
    
    ## Couple of useful visualizations
    rawvnormtab = cbind(rawcount = colSums(counttable), normcount = colSums(normcounttab), 
                        normtoraw_ratio = colSums(normcounttab)/colSums(counttable))
    rawvnormtab = rawvnormtab[order(rawvnormtab[,2], decreasing=TRUE), ]
    
    gridplot1 <- plot_histogram(data = as.data.frame(rawvnormtab[,1]), 
                                fitcurve=FALSE, binparam = NULL, limitx = NULL, 
                                labsparam = list(title = "Raw counts per Sample", x = "reads in sample", y = "Count"))
    
    
    gridplot2 <- plot_histogram(data = as.data.frame(rawvnormtab[,2]), 
                                fitcurve=FALSE, binparam = NULL, limitx = NULL, 
                                labsparam = list(title = "Normalized counts per Sample", x = "reads in sample", y = "Count"))
    
    pdf(paste0(outfilepath, "raw_count_per_sample_grid1.pdf"))
    print(gridplot1)
    junk <- dev.off()
    
    pdf(paste0(outfilepath, "norm_count_per_sample_grid2.pdf"))
    print(gridplot2)
    junk <- dev.off()

    
    indata = as.data.frame(rawvnormtab[,c(1,2)])
    if (nrow(indata) < 20) {datalabels = rownames(rawvnormtab)} else {datalabels = NULL}
    # if (label_extreme_changes == TRUE) {datalabels = rownames(rawvnormtab[1:(min(nrow(rawvnormtab),30)),])}
    if (label_extreme_changes == TRUE) {datalabels = rownames(rawvnormtab[abs(scale(rawvnormtab[,2])) > 1,])}
    labsparam = list(title = "raw to norm count comparison", x = "raw", y = "norm")
    outfile = paste(outfilepath, "norm_ratio_per_sample_grid3.pdf", sep="")
    comp_scat_plot = scatter_plotter(indata = indata, datalabels = datalabels, 
                                     labsparam = labsparam, plotstats = FALSE)
    pdf(file = outfile)
    suppressWarnings(print(comp_scat_plot))
    junk <- dev.off()
    
    readcounttab <- data.frame(colSums(normcounttab))
    scalereadcounttab <- scale(readcounttab)
    blowupsamples <- scalereadcounttab[abs(scalereadcounttab) > 2,1,drop=FALSE]
    colnames(blowupsamples) <- "norm_diff"
    
    write.table(normcounttab, file = paste(outfilepath, "normcounttab.txt", sep=""),  sep= "\t", col.names=NA, row.names = TRUE, quote=FALSE)
    
    return(list(normcounttab = normcounttab, deseqnorm_outliers = blowupsamples))
}
#####

## PCA Plotting - postNORM
#####

## Define the PCA plotting function
# badgenelist = c()
# pcadatanorm = t(normcounttab[!rownames(normcounttab) %in% badgenelist,])
# dir.create(paste(outfilepath, "pca_plots_normalized/", sep=""), recursive = TRUE, showWarnings = FALSE)
# for (desccol in 1:ncol(metatablefilt3)) {
#   colorvar = metatablefilt3[,desccol,drop=FALSE]
#   labsparam = c(list(title = "PCA", x = "PC1", y = "PC2", color=colnames(metatablefilt3)[desccol]))
#   outfile = paste(outfilepath, "pca_plots_normalized/", colnames(metatablefilt3)[desccol], "_pca_plot.pdf", sep="")
#   pca_plotter(pcadata = pcadatanorm, colorvar = colorvar, scalecolor=FALSE, labelpoints = FALSE, labsparam = labsparam, outfile = outfile)
# }

#####

## Apply Upper Quartile Normalization
#####
UpperQ_Normalization <- function(counttabfilt4, metatablefilt3, outfilepath, label_extreme_changes = FALSE) {
    upperquart = apply(counttabfilt4, 2, function(x){quantile(x, 0.75)})
    normcounttabUQ = t(t(counttabfilt4) / upperquart)
    
    ## Couple of useful visualizations
    rawvnormtabUQ = cbind(rawcount = colSums(counttabfilt4), normcount = colSums(normcounttabUQ), normtoraw_ratio = colSums(normcounttabUQ)/colSums(counttabfilt4))
    rawvnormtabUQ = rawvnormtabUQ[order(rawvnormtabUQ[,2], decreasing=TRUE), ]
    
    gridplot1 <- plot_histogram(data = as.data.frame(rawvnormtabUQ[,1]),
                                fitcurve=FALSE, binparam = NULL, limitx = NULL,
                                labsparam = list(title = "Raw counts per Sample", x = "reads in sample", y = "Count"))
    gridplot2 <- plot_histogram(data = as.data.frame(rawvnormtabUQ[,2]),
                                fitcurve=FALSE, binparam = NULL, limitx = NULL,
                                labsparam = list(title = "Normalized counts per Sample", x = "reads in sample", y = "Count"))
    
    pdf(paste0(outfilepath, "UQ_raw_count_per_sample_grid1.pdf"))
    print(gridplot1)
    junk <- dev.off()
    
    pdf(paste0(outfilepath, "UQ_norm_count_per_sample_grid2.pdf"))
    print(gridplot2)
    junk <- dev.off()
    
    indata = as.data.frame(rawvnormtabUQ[,c(1,2)])
    if (nrow(indata) < 20) {datalabels = rownames(rawvnormtabUQ)} else {datalabels = NULL}
    if (label_extreme_changes == TRUE) {datalabels = rownames(rawvnormtabUQ[1:(min(nrow(rawvnormtabUQ),40)),])}
    labsparam = list(title = "raw to UQ norm count comparison", x = "raw", y = "norm")
    outfile = paste(outfilepath, "UQ_norm_ratio_per_sample_grid3.pdf", sep="")
    comp_scat_plot = scatter_plotter(indata = indata, datalabels = datalabels, 
                                     labsparam = labsparam, plotstats = FALSE)
    pdf(file = outfile)
    suppressWarnings(print(comp_scat_plot))
    junk <- dev.off()
    
    write.table(normcounttabUQ, file = paste(outfilepath, "UQ_normcounttab.txt", sep=""),  sep= "\t", col.names=NA, row.names = TRUE, quote=FALSE)
    
    return(normcounttabUQ)
}
#####




## JUNK CODE
#####
