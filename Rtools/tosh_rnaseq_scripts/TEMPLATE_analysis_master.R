################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Master Script TEMPLATE for rnaseq analysis


## Load in Libraries
packagelist = c("ggplot2", "reshape2", "DESeq2", "grid", "gridExtra", "scales", "ggrepel", "tools")
junk <- lapply(packagelist, function(xxx) suppressMessages(
    require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))
ht_opt$message = FALSE

source("/Users/tosh/Desktop/Ruggles_Lab/code/rnaseq_scripts/rnaseq_processing_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/rnaseq_scripts/deseq_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/rnaseq_scripts/geneset_analysis_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/symbol_species_conversion_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/rnaseq_scripts/ssGSEA_custom.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/mgc_plotting_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/mgc_file_formatting.R")


## PROCESSING
countfilepath = ""
outfilepathmaster = ""
outfilepathprocess = paste0(outfilepathmaster, "rna_processing/")
dir.create(outfilepathprocess, recursive = TRUE, showWarnings = FALSE)

## Save and load data
# save.image(file = paste0(outfilepathmaster, "/PROJECT.RData"))
# load(file = paste0(outfilepathmaster, "/PROJECT.RData"))

## Process count files into count table
counttab = create_count_table(countfilepath = countfilepath, outfilepath = outfilepathprocess)
#######OR#######
## Read in Count table
counttab = read.table(countfilepath, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE , comment = "")

## Remove "_S" string if necessary
colnames(counttab) = gsub("_[^_]+$", "", colnames(counttab))
#colnames(counttab) = sapply(strsplit(colnames(counttab), split = "_S"), function(x) (x[1]))

## Reading in Metadata and filter by the metadata
metafile = ""
metatable = read.table(metafile, sep = ",", header=TRUE, stringsAsFactors = FALSE, na.strings = c(NA, "NA", ""))

metatable[,1] = as.character(metatable[,1])
colnames(metatable) = c("sample", colnames(metatable)[2:ncol(metatable)])

## Find the overlapping SAMPLES between the metatable and countdata - then select for that
sampleids = sort(intersect(metatable[,1], colnames(counttab)))
sampleidsordered = sampleids[match(metatable[,1][metatable[,1] %in% colnames(counttab)], sampleids)]

## Pull out the samples that only fell in one group or the other - for my reference
samplesONLYinmeta <- setdiff(metatable[,1], colnames(counttab))
samplesONLYincounts <- setdiff(colnames(counttab), metatable[,1])
if (length(samplesONLYinmeta) > 0) {print(
    paste0("Warning: there were ", length(samplesONLYinmeta), " samples in the meta NOT in the counttable"))}
if (length(samplesONLYincounts) > 0) {print(
    paste0("Warning: there were ", length(samplesONLYincounts), " samples in the counttable NOT in the meta"))}

metatablefilt1 = metatable[match(sampleidsordered, metatable[,1]),]
counttabfilt1 = counttab[,sampleidsordered]
## Check to make sure data lines up
identical(colnames(counttabfilt1), metatablefilt1[,1])
dim(metatablefilt1)
dim(counttabfilt1)
write.table(counttabfilt1, file = paste0(outfilepathprocess, "rawcounttab.txt"), 
            sep= "\t", col.names=NA, row.names = TRUE, quote=FALSE)


## Add in a new output - where we will flag samples throughout the plotting and QC steps
## Then output the flag table, and essentially say that these are the samples which should be double checked
## Tests so far include
# 5 - the density curves == densityoutliers
# 1 - the min number of reads == QC_readcount_outliers
# 3 - PCA outliers (3SD away from the mean) == pca_outliers
# 2 - The blow up after normalization (3SD above the mean norm count) == deseqnorm_outliers
# 4 - SS correlation outliers (3SD from the mean correlation sum) == correlation_outliers

## Remove bad samples from previous run - doesnt matter the file as long as the first column are IDs that match
# badsampfile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2020/data/metadata_sheets/possible_outlier_tab_1_20200108.csv"
# badsamptab <- read.table(badsampfile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
# badsamps <- badsamptab[,1]
# 
# ## filter the initial counttab to exclude those samples
# counttabfilt1 <- counttabfilt1[,!colnames(counttabfilt1) %in% badsamps]



##### Plot Density Curves
density_curve_out <- plot_read_density_curves(counttabfilt1, outfilepathprocess)
densityplot <- density_curve_out$density_plot
densityoutliers <- density_curve_out$density_outliers

##### Filter for Read Count
minreadcutoff <- 1000000
QC_filter_readcount_out <- QC_filter_readcount(counttabfilt1, minreadcutoff)
counttabfilt2 <- QC_filter_readcount_out$counttabfilt
QC_readcount_outliers <- QC_filter_readcount_out$outliercounts

pdf(paste0(outfilepathprocess, "prefiltered_samp_readcount_hist.pdf"))
print(QC_filter_readcount_out$prefilterhist)
junk <- dev.off()

pdf(paste0(outfilepathprocess, "postfiltered_samp_readcount_hist.pdf"))
print(QC_filter_readcount_out$postfilterhist)
junk <- dev.off()

##### Filter for Gene Count
mincountcutoff <- 4
minsamplescutoff <- round(ncol(counttabfilt2)/2)
QC_filter_genecount_out <- QC_filter_genecount(counttabfilt2, mincountcutoff, minsamplescutoff)
counttabfilt3 <- QC_filter_genecount_out$counttabfilt

pdf(paste0(outfilepathprocess, "prefiltered_gene_rowmeans_hist.pdf"))
print(QC_filter_genecount_out$prefilterhist)
junk <- dev.off()

pdf(paste0(outfilepathprocess,  "postfiltered_gene_rowmeans_hist.pdf"))
print(QC_filter_genecount_out$postfilterhist)
junk <- dev.off()

## Adding in additional metadata
# Add in the readcount as a metadata column
# readcounttab = data.frame(sampreadsum = colSums(counttabfilt2))
# readcounttab2 = cbind.data.frame(sample = gsub(pattern = ".*_","",
#                                                rownames(readcounttab)), readcount = readcounttab[,1,drop=TRUE])

# Compile all new metadata columns:
newmetalist <- list()

## Combine new data onto metatable
metatablefilt2 = Reduce(function(dtf1, dtf2)
    merge(dtf1, dtf2, by = "sample", all.x = TRUE, sort = FALSE), c(list(metatablefilt1), newmetalist))
rownames(metatablefilt2) = metatablefilt2[,1]
metatablefilt3 = metatablefilt2[metatablefilt2[,1] %in% colnames(counttabfilt3),2:ncol(metatablefilt2)]
write.table(metatablefilt3, file = paste0(outfilepathprocess, "metatable_in.csv"), 
            sep = ",", col.names = NA, row.names = TRUE, quote = FALSE)


## THESE ARE OUR FINAL TABLES
## Extract our comparison columns for now
compcols = metatablefilt3[,c(grepl("comp_", colnames(metatablefilt3))), drop=FALSE]
metatablefilt4 = metatablefilt3[!grepl("comp_", colnames(metatablefilt3))]

## Change non-character columns to characters for plotting characterization
columns_to_characterify = c()
metatablefilt4[,columns_to_characterify] = apply(metatablefilt4[,columns_to_characterify,drop=FALSE], 2, as.character)

counttabfilt4 = counttabfilt3[,match(colnames(counttabfilt3), rownames(metatablefilt4))]


## Write out the filtered count table and metadata for further process steps (other scripts)
write.table(counttabfilt4, file = paste0(outfilepathprocess, "filtrawcounttab.txt"), 
            sep= "\t", col.names=NA, row.names = TRUE, quote=FALSE)
write.table(metatablefilt4, file = paste0(outfilepathprocess, "metatable_filt.txt"), 
            sep= "\t", col.names=NA, row.names = TRUE, quote=FALSE)

print(dim(counttabfilt4))
print(dim(metatablefilt4))


##### Plot the top genes that are proportionately representative in our dataset
numgenesplot = 50
gene_prop_out <- plot_gene_prop_histogram(counttabfilt4 = counttabfilt4, numgenesplot = numgenesplot)

out_geneprop_file <- paste0(outfilepathprocess, "read_distribution_bar_chart.pdf")
pdf(out_geneprop_file, height = 15, width = max(10,round(ncol(counttabfilt4)/10)))
print(gene_prop_out$geneprop_plot)
grid.newpage()
grid.draw(gene_prop_out$geneprop_plotlegend)
junk <- dev.off()


## Create a table of top varied genes for plotting on global heatmaps and PCAs
upperpercentile = 0.9999
lowerpercentile = 0.95
counttabfilt4_var <- sort(apply(counttabfilt4,1,var), decreasing=TRUE)
## What if I take the 99th - 90th percentile, therefore omitting the crazy drivers at the top, and leaving the rest...
vargeneselect <- names(counttabfilt4_var)[
    round((1-upperpercentile)*length(counttabfilt4_var)):round((1-lowerpercentile)*length(counttabfilt4_var))]
topvartab <- counttabfilt4[vargeneselect,]


##### PCA PLOTTING
outfilepathplot = paste0(outfilepathmaster, "plotting/")
dir.create(outfilepathplot, recursive = TRUE, showWarnings = FALSE)

# pcadata = t(counttabfilt4)
pcadata = t(topvartab)

dir.create(paste0(outfilepathplot, "pca_plots/"), recursive = TRUE, showWarnings = FALSE)
for (desccol in 1:ncol(metatablefilt4)) {
    colorvar = metatablefilt4[,desccol,drop=FALSE]
    labsparam = c(list(title = "PCA", x = "PC1", y = "PC2", color=colnames(metatablefilt4)[desccol]))
    outfile = paste(outfilepathplot, "pca_plots/", colnames(metatablefilt4)[desccol], "_pca_plot.pdf", sep="")
    pcaplotout <- pca_plotter(pcadata = pcadata, colorvar = colorvar, scalecolor=FALSE, separatelegend = FALSE,
                              labelpoints = FALSE, labsparam = labsparam)
    
    pdf(file = outfile)
    print(pcaplotout$pca_out)
    junk <- dev.off()
}
pca_outliers = pcaplotout$pca_outliers

##### DESEQ NORMALIZATION
print("Checks for DESeq to be run correctly - both should be true")
all(rownames(metatablefilt4) %in% colnames(counttabfilt4))
all(rownames(metatablefilt4) == colnames(counttabfilt4))

DEseq_Normalization_out = DEseq_Normalization(counttable = counttabfilt4, metatable = metatablefilt4, 
                                              outfilepath = outfilepathprocess, label_extreme_changes = TRUE)
normcounttab = DEseq_Normalization_out$normcounttab
deseqnorm_outliers = DEseq_Normalization_out$deseqnorm_outliers


## Create a table of top varied genes for plotting on global heatmaps and PCAs
normcounttabfilt4_var <- sort(apply(normcounttab,1,var), decreasing=TRUE)
## What if I take the 99th - 90th percentile, therefore omitting the crazy drivers at the top, and leaving the rest...
normvargeneselect <- names(normcounttabfilt4_var)[
    round((1-upperpercentile)*length(normcounttabfilt4_var)):round((1-lowerpercentile)*length(normcounttabfilt4_var))]
normtopvartab <- normcounttab[normvargeneselect,]


##### POST NORMALIZATION PCA
# pcadatanorm = t(normcounttab)
pcadatanorm = t(normtopvartab)

dir.create(paste0(outfilepathplot, "pca_plots_normalized/"), recursive = TRUE, showWarnings = FALSE)
for (desccol in 1:ncol(metatablefilt4)) {
    colorvar = metatablefilt4[,desccol,drop=FALSE]
    labsparam = c(list(title = "PCA", x = "PC1", y = "PC2", color=colnames(metatablefilt4)[desccol]))
    outfile = paste(outfilepathplot, "pca_plots_normalized/", 
                    colnames(metatablefilt4)[desccol], "_pca_plot.pdf", sep="")
    pcaplotout <- pca_plotter(pcadata = pcadatanorm, colorvar = colorvar, scalecolor=FALSE, separatelegend = FALSE,
                              labelpoints = FALSE, labsparam = labsparam)
    
    pdf(file = outfile)
    print(pcaplotout$pca_out)
    junk <- dev.off()
}

##### UQ Normalization
normcounttabUQ = UpperQ_Normalization(counttabfilt4, metatablefilt4, outfilepathprocess)


##### HEATMAP PLOTTING
dir.create(paste0(outfilepathplot, "heatmaps/"), recursive = TRUE, showWarnings = FALSE)
custcolorlist = NULL

hmpdfoutfile = paste0(outfilepathplot, "heatmaps/hm1_heatmap_allmetadata.pdf")
annotationlist1 = annotationlist_builder(metatablefilt4, customcolorlist = custcolorlist)
outhm1 <- create_heatmap(counttab = normtopvartab, subsetnum = FALSE, colmetatable = metatablefilt4, colannotationlist = annotationlist1, 
                         colclusterparam = TRUE, rowclusterparam = TRUE, separate_legend = TRUE)
pdf(hmpdfoutfile, height = 11.5, width = 10)
draw(outhm1[[1]])
grid.newpage()
draw(outhm1[[2]])
junk <- dev.off()


## Iteratively sort by each metadata entry - and then display by that
for (colnum in 1:ncol(metatablefilt4)) {
    metatablesort = metatablefilt4[order(metatablefilt4[,colnum]),,drop=FALSE]
    metatablesort <- metatablesort[rownames(na.omit(metatablesort[,colnum,drop=FALSE])),,drop=FALSE]
    # normcounttabsort = normcounttab[,rownames(metatablesort)]
    normcounttabsort = normtopvartab[,rownames(metatablesort)]
    
    annotationlist = annotationlist_builder(metatablesort, customcolorlist = custcolorlist)
    outheatmapfile = paste0("hm", colnum, "_heatmap_sorted_", colnames(metatablesort[,colnum,drop=FALSE]),".pdf")
    outhm1 <- create_heatmap(counttab = normcounttabsort, subsetnum = FALSE, 
                             colmetatable = metatablesort, colannotationlist = annotationlist, 
                             colclusterparam = FALSE, rowclusterparam = TRUE, separate_legend = TRUE)
    
    pdfoutfile = paste0(outfilepathplot, "heatmaps/", outheatmapfile)
    pdf(pdfoutfile, height = 11.5, width = 10)
    draw(outhm1[[1]])
    grid.newpage()
    draw(outhm1[[2]])
    junk <- dev.off()
    
}


## Output a sample-sample heatmap for global correlation and outlier detection
SS_outheatmapfile <- paste0(outfilepathplot, "heatmaps/", "hmSS_heatmap_sample_sample.pdf")
outSShm <- create_SS_heatmap(counttab = normtopvartab, metatable = metatablefilt4, annotationlist = annotationlist1, 
                             rowclusterparam=TRUE, colclusterparam=TRUE, separate_legend = TRUE)
# outSShm <- create_SS_heatmap(counttab = normcounttab, metatable = metatablefilt4, annotationlist = annotationlist1, 
#                              rowclusterparam=TRUE, colclusterparam=TRUE, separate_legend = TRUE)
pdfoutfile = paste0(outfilepathplot, outheatmapfile)
pdf(SS_outheatmapfile, height = 11.5, width = 10)
draw(outSShm[[1]])
grid.newpage()
draw(outSShm[[2]])
junk <- dev.off()



#If they are more off then the rest (again, lets say 3 SDs) then they should be flagged)
cortab <- cor(normcounttab, method = "spearman")
scalecortab <- scale(colSums(cortab))
correlation_outliers <- scalecortab[abs(scalecortab) > 3, 1, drop = FALSE]
colnames(correlation_outliers) <- "corr_diff"


### OUTPUT SUGGESTED OUTLIERS
outlierlist <- list(densityoutliers, QC_readcount_outliers, pca_outliers, deseqnorm_outliers, correlation_outliers)
outlierlist2 <- lapply(outlierlist, function(x) {
    temp = cbind(samples = rownames(x), metric = x[,1])
    colnames(temp) = c("samples", colnames(x[,1,drop=FALSE]))
    temp})
## Combine new data onto metatable
suggested_outliertab = Reduce(function(dtf1, dtf2)
    merge(dtf1, dtf2, by = "samples", all = TRUE, sort = FALSE), c(outlierlist2))
outliertab_outfile <- paste0(outfilepathprocess, "possible_outlier_tab.csv")

## Summarize results by saying how many tests the sample failed
suggested_outliertab[,"cumulative_outlier_score"] <- rowSums(!is.na(suggested_outliertab[,2:ncol(suggested_outliertab)]))

write.table(suggested_outliertab, outliertab_outfile, sep = ",", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)


##### DESEQ ANALYSIS
outfilepathDEseq = paste0(outfilepathmaster, "deseq/")
dir.create(outfilepathDEseq, recursive = TRUE, showWarnings = FALSE)

#controlcols <- metatablefilt4[,controlcols,drop=FALSE]
# controlcols <- apply(controlcols, 2, function(x) gsub("-", "_", x))
controlcols <- NULL

DEseqreslist = DESeq_analysis(compcols = compcols, controlcols = controlcols, rawcounttab = counttabfilt4)

##### RUN PAIRED ANALYSIS
# The best way to do this is run an additional analysis and just put in what columns we want to pair off of
# Then append that to the DEseqreslist and everything else should work smoothly

# # Followup vs baseline
# paircompname <- "comp_followups__followup_vs_baseline"
# 
# ## How pairs are assigned has to be highly manual and vary situation to situation!!!
# paircontrolcols <- na.omit(compcols[order(rownames(metatablefilt4)),paircompname,drop=FALSE])
# pairedsamps <- gsub("-1|-2", "", rownames(paircontrolcols))
# paircheck <- table(pairedsamps)
# ## Need to remove those who dont actually have a mate (QC removed or whatever)
# paircheckpass <- paircheck[paircheck == 2]
# paircheckfail <- paircheck[paircheck != 2]
# 
# ## Now select the rows that pass check
# paircompcols <- compcols[grepl(paste(names(paircheckpass), collapse = "|"), rownames(compcols)),
#                          paircompname,drop=FALSE]
# colnames(paircompcols) <- c("PAIRcomp_followups__followup_vs_baseline")
# ## Define the pairs to control against
# paircontrolcols <- data.frame(PAIRparam = paste0("pair", rep(1:length(paircheckpass), each = 2)), row.names = rownames(paircompcols))
# 
# PAIRDEseqreslist = DESeq_analysis(compcols = paircompcols, controlcols = paircontrolcols, rawcounttab = counttabfilt4)
# 
# ## I need to add in the new paired analysis column. But I'll add a failsafe to make sure it only happens once
# if (!colnames(paircompcols) %in% colnames(compcols)) { 
#     temp <- merge(compcols, paircompcols, all.x = TRUE, by = "row.names")
#     rownames(temp) <- temp[,1]
#     compcols <- temp[,-1]
# }
# 
# ## Now just add on the paired analysis and finish the pipeline
######## MUST REORDER AGAIN USING ORIGINAL COMPCOLS ORIENTATION!!!!
# DEseqreslist <- c(DEseqreslist, PAIRDEseqreslist)[colnames(compcols)]

## Write out DEseq analysis
for (DEseqrestab in seq_len(length(DEseqreslist))) {
    deseqanalysisoutfolder = paste0(outfilepathDEseq, "/", names(DEseqreslist)[DEseqrestab], "/")
    dir.create(deseqanalysisoutfolder, recursive = TRUE, showWarnings = FALSE)
    deseqanalysisoutfile = paste0(deseqanalysisoutfolder, "deseq_results_", names(DEseqreslist)[DEseqrestab], ".csv")
    write.table(DEseqreslist[[DEseqrestab]], deseqanalysisoutfile, quote = FALSE, sep = ",", row.names = TRUE, col.names=NA)
}

## Summary Table and Genelists
#####
# stattype = {pvalue | padj}
summaryparams = list(run1 = c("stattype" = "pvalue", "statcutoff" = 0.01, "log2fccutoff" = 1),
                     run2 = c("stattype" = "pvalue", "statcutoff" = 0.01, "log2fccutoff" = 0),
                     run3 = c("stattype" = "pvalue", "statcutoff" = 0.05, "log2fccutoff" = 1),
                     run4 = c("stattype" = "pvalue", "statcutoff" = 0.05, "log2fccutoff" = 0),
                     run5 = c("stattype" = "padj", "statcutoff" = 0.05, "log2fccutoff" = 1),
                     run6 = c("stattype" = "padj", "statcutoff" = 0.05, "log2fccutoff" = 0),
                     run7 = c("stattype" = "padj", "statcutoff" = 0.1, "log2fccutoff" = 1),
                     run8 = c("stattype" = "padj", "statcutoff" = 0.1, "log2fccutoff" = 0))
DESeq_summary(DEseqreslist, summaryparams, outfilepathDEseq)

## Volcano plot and GOI heatmap
# pvalcutoffparam = 0.1
# log2fccutoffparam = 1
# stattype = "padj"
for (deseqanalysisnum in seq_len(length(DEseqreslist))) {
    
    ## Draw a plot to show number of genes available at various cutoffs
    # Run analysis for pvalue
    subtab1 = DEseqreslist[[deseqanalysisnum]][,c("log2FoldChange", "pvalue")]
    deseqanalysislabel = names(DEseqreslist[deseqanalysisnum])
    destatplot_out_pval <- destatplot(subtab = subtab1)
    pdf(paste0(outfilepathDEseq, deseqanalysislabel, "/dge_pvalue_statplot.pdf"))
    print(destatplot_out_pval[["statplot"]])
    junk <- dev.off()
    write.table(destatplot_out_pval[["stattab"]], paste0(outfilepathDEseq, deseqanalysislabel, "/dge_pvalue_statplot_tab.txt"), 
                col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")
    
    # Run analysis for padj
    subtab2 = DEseqreslist[[deseqanalysisnum]][,c("log2FoldChange", "padj")]
    destatplot_out_padj <- destatplot(subtab = subtab2)
    pdf(paste0(outfilepathDEseq, deseqanalysislabel, "/dge_padj_statplot.pdf"))
    print(destatplot_out_padj[["statplot"]])
    junk <- dev.off()
    write.table(destatplot_out_padj[["stattab"]], paste0(outfilepathDEseq, deseqanalysislabel, "/dge_padj_statplot_tab.txt"), 
                col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")  
    
    ## Create summary figures for all of the requested summary params
    for (runnum in seq_len(length(summaryparams))) {
        ## Define the stats for the run and create the outfolder
        stattype = unname(summaryparams[[runnum]]["stattype"])
        pvalcutoffparam = as.numeric(unname(summaryparams[[runnum]]["statcutoff"]))
        log2fccutoffparam = as.numeric(unname(summaryparams[[runnum]]["log2fccutoff"]))
        
        subtab = DEseqreslist[[deseqanalysisnum]]
        deseqanalysislabel = names(DEseqreslist[deseqanalysisnum])
        intable = subtab[,c("log2FoldChange", stattype)]
        
        pout1 = create_volcano_plot(intable, pvalcutoff = pvalcutoffparam, 
                                    log2fccutoff = log2fccutoffparam, labeledgenes = TRUE, 
                                    nameparam = paste0(deseqanalysislabel, "_volcano_plot"))
        pdf(paste0(outfilepathDEseq, "/", deseqanalysislabel, "/", stattype, pvalcutoffparam, "_lfc", log2fccutoffparam, "_", "volcano_plot.pdf"))
        print(pout1)
        dev.off()
        
        metatableGOI = compcols[,deseqanalysisnum,drop=FALSE]
        metatableGOI = na.omit(metatableGOI[order(metatableGOI[,1]),,drop=FALSE])
        metatableGOI[,1] = as.character(metatableGOI[,1])
        
        annotationlist = annotationlist_builder(metatableGOI)
        
        GOIhmtab = normcounttab[rownames(na.omit(intable[intable[,2] < pvalcutoffparam & abs(intable[,1]) > log2fccutoffparam,,drop=FALSE])),,drop=FALSE]
        GOIhmtab = GOIhmtab[order(na.omit(intable[intable[,2] < pvalcutoffparam & abs(intable[,1]) > log2fccutoffparam,1,drop=FALSE])), 
                            rownames(metatableGOI),drop=FALSE]
        
        ## Draw GOI heatmaps with and without column clustering
        if (nrow(GOIhmtab) > 0) {
            GOIhm1 <- create_heatmap(counttab = GOIhmtab, colmetatable = metatableGOI, colannotationlist = annotationlist,
                                     colclusterparam = FALSE, rowclusterparam = TRUE)
            pdfoutfile = paste0(outfilepathDEseq, "/", deseqanalysislabel, "/", stattype, pvalcutoffparam, "_lfc", log2fccutoffparam, "_", "GOI_heatmap_1.pdf")
            pdf(pdfoutfile, height = 11.5, width = 10)
            draw(GOIhm1[[1]])
            junk <- dev.off()
            
            GOIhm2 <-create_heatmap(counttab = GOIhmtab, colmetatable = metatableGOI, colannotationlist = annotationlist,
                                    colclusterparam = TRUE, rowclusterparam = TRUE)
            pdfoutfile = paste0(outfilepathDEseq, "/", deseqanalysislabel, "/", stattype, pvalcutoffparam, "_lfc", log2fccutoffparam, "_", "GOI_heatmap_2.pdf")
            pdf(pdfoutfile, height = 11.5, width = 10)
            draw(GOIhm2$heatmap)
            junk <- dev.off()
            
            ## Add the metadata back in to make sure that our GOI clustering isnt due to some other confounding variable
            metatableGOI2 <- cbind(metatablefilt4[rownames(metatableGOI),], metatableGOI)
            ## HOT FIX TO PROPERLY CODE LATER - IF THERE ARE PURE NA COLUMNS, THEN DONT PLOT THEM
            metatableGOI2<- metatableGOI2[colSums(!is.na(metatableGOI2)) > 0]
            annotationlist2 = annotationlist_builder(metatableGOI2)
            
            GOIhm3 <-create_heatmap(counttab = GOIhmtab, colmetatable = metatableGOI2, colannotationlist = annotationlist2,
                                    colclusterparam = TRUE, rowclusterparam = TRUE, separate_legend = TRUE)
            pdfoutfile = paste0(outfilepathDEseq, "/", deseqanalysislabel, "/", stattype, pvalcutoffparam, "_lfc", log2fccutoffparam, "_", "GOI_heatmap_3.pdf")
            pdf(pdfoutfile, height = 11.5, width = 10)
            draw(GOIhm3$heatmap)
            junk <- dev.off()
            
            ## Special Paired analysis GOI heatmaps
            if (grepl("PAIR", deseqanalysislabel)) {
                ## Need to have a metatable with the comparison and the pairs, then create an annot list off that
                metatableGOI3 <- merge(metatableGOI, paircontrolcols, by = "row.names")
                rownames(metatableGOI3) <- metatableGOI3[,1]
                metatableGOI3 <- metatableGOI3[,-1]
                annotationlist3 = annotationlist_builder(metatableGOI3)
                GOIhmtab2 <- GOIhmtab[,rownames(metatableGOI3)]
                
                GOIhm4 <- create_heatmap(counttab = GOIhmtab2, colmetatable = metatableGOI3, colannotationlist = annotationlist3,
                                         colclusterparam = FALSE, rowclusterparam = TRUE)
                pdfoutfile = paste0(outfilepathDEseq, "/", deseqanalysislabel, "/", 
                                    stattype, pvalcutoffparam, "_lfc", log2fccutoffparam, "_", "GOI_heatmap_4.pdf")
                pdf(pdfoutfile, height = 11.5, width = 10)
                draw(GOIhm4[[1]])
                junk <- dev.off()
                
                GOIhm5 <- create_heatmap(counttab = GOIhmtab2, colmetatable = metatableGOI3, colannotationlist = annotationlist3,
                                         colclusterparam = TRUE, rowclusterparam = TRUE)
                pdfoutfile = paste0(outfilepathDEseq, "/", deseqanalysislabel, "/", 
                                    stattype, pvalcutoffparam, "_lfc", log2fccutoffparam, "_", "GOI_heatmap_5.pdf")
                pdf(pdfoutfile, height = 11.5, width = 10)
                draw(GOIhm5[[1]])
                junk <- dev.off()
                
                GOIpcadata = t(GOIhmtab)
                colorvar = metatableGOI3[rownames(GOIpcadata),2,drop=FALSE]
                shapevar = metatableGOI3[rownames(GOIpcadata),1,drop=FALSE]
                labsparam = c(list(title = "PCA", x = "PC1", y = "PC2", color=colnames(metatableGOI)))
                outfile = paste0(outfilepathDEseq, "/", deseqanalysislabel, "/", stattype, pvalcutoffparam, "_lfc", log2fccutoffparam, "_", "GOI_pca_paired.pdf")
                pcaplotout <- pca_plotter(pcadata = GOIpcadata, colorvar = colorvar, shapevar = shapevar, scalecolor=FALSE, separatelegend = FALSE,
                                          labelpoints = TRUE, labsparam = labsparam, returnoutliers = FALSE)
                
                pdf(file = outfile)
                print(pcaplotout$pca_out)
                junk <- dev.off()
            }
            
            if (nrow(GOIhmtab) > 1) {
                GOIpcadata = t(GOIhmtab)
                # for (desccol in 1:ncol(metatablefilt4)) {
                colorvar = metatableGOI
                labsparam = c(list(title = "PCA", x = "PC1", y = "PC2", color=colnames(metatableGOI)))
                outfile = paste0(outfilepathDEseq, "/", deseqanalysislabel, "/", stattype, pvalcutoffparam, "_lfc", log2fccutoffparam, "_", "GOI_pca.pdf")
                pcaplotout <- pca_plotter(pcadata = GOIpcadata, colorvar = colorvar, scalecolor=FALSE, separatelegend = FALSE,
                                          labelpoints = FALSE, labsparam = labsparam, returnoutliers = FALSE)
                
                pdf(file = outfile)
                print(pcaplotout$pca_out)
                junk <- dev.off()
            }
        }
    }
    print(paste0("Finished Summary for ", deseqanalysislabel))
}




## Gene set analysis
outfilepathGSEA = paste0(outfilepathmaster, "gsea/")
dir.create(outfilepathGSEA, recursive = TRUE, showWarnings = FALSE)
## {"Homo sapiens", "Mus musculus"}
speciesparam = "Homo sapiens"
pstatparam = "pvalue"
numpathways_plotparam <- 10

## Seeding to get rid of the change in results between runs
seedparam = NULL
#sample(1:2147483647, 1)
#seedparam = 

for (deseqanalysisnum in seq_len(length(DEseqreslist))) {
    ## Create folder for each comparison
    outfilepathGSEAcomp <- paste0(outfilepathGSEA, names(DEseqreslist[deseqanalysisnum]), "/")
    dir.create(outfilepathGSEAcomp, recursive = TRUE, showWarnings = FALSE)
    
    ## GSEA run through for KEGG terms
    if (!is.null(seedparam)) {set.seed(seedparam)}
    geneset_analysis_out_KEGG = geneset_analysis(DEseqtable = DEseqreslist[[deseqanalysisnum]], 
                                                 pvalcutoffparam = 1, genesetparam = c("CP:KEGG"), speciesparam = speciesparam,
                                                 seedparam = seedparam)
    write.table(geneset_analysis_out_KEGG, 
                file = paste0(outfilepathGSEAcomp, names(DEseqreslist)[[deseqanalysisnum]], "_gsea_KEGG.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    KEGG_plot_out = gsea_barplot(gseaout = geneset_analysis_out_KEGG, 
                                 pstatparam = pstatparam, numterms = numpathways_plotparam, 
                                 titleparam = names(DEseqreslist)[deseqanalysisnum])
    pdf(width = 11, height = 8.5, paste0(outfilepathGSEAcomp, names(DEseqreslist)[[deseqanalysisnum]], "_gsea_KEGG_plot.pdf"))
    print(KEGG_plot_out)
    junk <- dev.off()
    
    ## GSEA run through for HALLMARK terms
    if (!is.null(seedparam)) {set.seed(seedparam)}
    geneset_analysis_out_HALL = geneset_analysis(DEseqreslist[[deseqanalysisnum]], 
                                                 pvalcutoffparam = 1, genesetparam = c("H"), speciesparam = speciesparam,
                                                 seedparam = seedparam)
    write.table(geneset_analysis_out_HALL, 
                file = paste0(outfilepathGSEAcomp, names(DEseqreslist)[[deseqanalysisnum]], "_gsea_HALL.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    HALL_plot_out = gsea_barplot(geneset_analysis_out_HALL, 
                                 pstatparam = pstatparam, numterms = numpathways_plotparam, 
                                 titleparam = names(DEseqreslist)[deseqanalysisnum])
    pdf(width = 11, height = 8.5, paste0(outfilepathGSEAcomp, names(DEseqreslist)[[deseqanalysisnum]], "_gsea_HALL_plot.pdf"))
    print(HALL_plot_out)
    junk <- dev.off()
    
    ## GSEA run through for GO terms
    if (!is.null(seedparam)) {set.seed(seedparam)}
    geneset_analysis_out_GO = geneset_analysis(DEseqreslist[[deseqanalysisnum]], 
                                               pvalcutoffparam = 1, genesetparam = c("C5"), speciesparam = speciesparam,
                                               seedparam = seedparam)
    write.table(geneset_analysis_out_GO, 
                file = paste0(outfilepathGSEAcomp, names(DEseqreslist)[[deseqanalysisnum]], "_gsea_GO.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    GO_plot_out = gsea_barplot(geneset_analysis_out_GO, 
                               pstatparam = pstatparam, numterms = numpathways_plotparam, 
                               titleparam = names(DEseqreslist)[deseqanalysisnum])
    pdf(width = 11, height = 8.5, paste0(outfilepathGSEAcomp, names(DEseqreslist)[[deseqanalysisnum]], "_gsea_GO_plot.pdf"))
    print(GO_plot_out)
    junk <- dev.off()
    
}

## GSEA hypergeometric test
for (deseqanalysisnum in seq_len(length(DEseqreslist))) {
    ## Create folder for each comparison
    outfilepathGSEAcomp <- paste0(outfilepathGSEA, names(DEseqreslist[deseqanalysisnum]), "/")
    dir.create(outfilepathGSEAcomp, recursive = TRUE, showWarnings = FALSE)
    
    statcutoffparamlist = c("stattype" = "pvalue", "pstatcutoff" = 0.01, "log2fccutoff" = 0)
    hypergeo_genetest_out_KEGG = hypergeo_genetest(DEseqtable = DEseqreslist[[deseqanalysisnum]],
                                                   statcutoffparam = statcutoffparamlist, 
                                                   genesetparam = c("CP:KEGG"), speciesparam = speciesparam)
    write.table(hypergeo_genetest_out_KEGG$enricherUPout, 
                file = paste0(outfilepathGSEAcomp, names(DEseqreslist)[[deseqanalysisnum]], "_hypergeo_UP_KEGG.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    write.table(hypergeo_genetest_out_KEGG$enricherDOWNout, 
                file = paste0(outfilepathGSEAcomp, names(DEseqreslist)[[deseqanalysisnum]], "_hypergeo_DOWN_KEGG.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    
    hypergeo_genetest_out_HALL = hypergeo_genetest(DEseqreslist[[deseqanalysisnum]],
                                                   statcutoffparam = statcutoffparamlist, 
                                                   genesetparam = c("H"), speciesparam = speciesparam)
    write.table(hypergeo_genetest_out_HALL$enricherUPout, 
                file = paste0(outfilepathGSEAcomp, names(DEseqreslist)[[deseqanalysisnum]], "_hypergeo_UP_HALL.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    write.table(hypergeo_genetest_out_HALL$enricherDOWNout, 
                file = paste0(outfilepathGSEAcomp, names(DEseqreslist)[[deseqanalysisnum]], "_hypergeo_DOWN_HALL.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    
    hypergeo_genetest_out_GO = hypergeo_genetest(DEseqreslist[[deseqanalysisnum]],
                                                 statcutoffparam = statcutoffparamlist, 
                                                 genesetparam = c("C5"), speciesparam = speciesparam)
    write.table(hypergeo_genetest_out_GO$enricherUPout, 
                file = paste0(outfilepathGSEAcomp, names(DEseqreslist)[[deseqanalysisnum]], "_hypergeo_UP_GO.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    write.table(hypergeo_genetest_out_GO$enricherDOWNout, 
                file = paste0(outfilepathGSEAcomp, names(DEseqreslist)[[deseqanalysisnum]], "_hypergeo_DOWN_GO.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    
}




# ## Plot pathways of interest:
# keggpathwaytab = as.data.frame(msigdbr(species = speciesparam, subcategory = c("CP:KEGG"))[,c("gs_name", "gene_symbol")])
# hallmarkpathwaytab = as.data.frame(msigdbr(species = speciesparam, category = c("H"))[,c("gs_name", "gene_symbol")])
# gopathwaytab = as.data.frame(msigdbr(species = speciesparam, category = c("C5"))[,c("gs_name", "gene_symbol")])
# pathwaytab = rbind(keggpathwaytab, hallmarkpathwaytab, gopathwaytab)
# deseqtabname <- ""
# plotmetacols <- c(1:4)
# 
# ## NOTE - GO_RNA_PROCESSING is no longer valid with new msigdb - this has been replaced with specific RNA processing pathways
# go_pathwayparam = c()
# hallmark_pathwayparam = c()
# kegg_pathwayparam = c()
# pathwayparam = c(go_pathwayparam, hallmark_pathwayparam, kegg_pathwayparam)
# 
# for (pathwaynum in seq_len(length(pathwayparam))) {
#     pathwaylabel = pathwayparam[pathwaynum]
#     pathwaygenes = pathwaytab[pathwaytab[,1] == pathwaylabel,2]
#     
#     ##amazingly.. theres a bug where some pathways have duplicate genes.. so I am going to remove those
#     pathwaygenes = unique(pathwaygenes)
#     
#     # HEATMAP PLOTTING
#     outfilepathplot = paste0(outfilepathmaster, "GOI/pathway_plotting/", deseqtabname, "/")
#     dir.create(outfilepathplot, recursive = TRUE, showWarnings = FALSE)
#     
#     plotmeta = metatablefilt4[,plotmetacols,drop=FALSE]
#     plotmeta = plotmeta[order(plotmeta[,1]),,drop=FALSE]
#     plottab = normcounttab[pathwaygenes[pathwaygenes %in% rownames(normcounttab)],rownames(plotmeta), drop = FALSE]
#     
#     # Adding in the DESeq information to sort the genes by the log2fc
#     log2fcval = DEseqreslist[[deseqtabname]][pathwaygenes[pathwaygenes %in% rownames(DEseqreslist[[1]])],c("log2FoldChange", "pvalue"),drop=FALSE]
#     geneorder = log2fcval[order(log2fcval[,2], decreasing = FALSE),,drop=FALSE]
#     plottab = plottab[rownames(geneorder),,drop=FALSE]
#     
#     
#     outheatmapfile = paste0(outfilepathplot, pathwaylabel, ".pdf")
#     annotationlist1 = annotationlist_builder(plotmeta,
#                                              customcolorlist = list(group = c("stress" = "red", "control" = "blue")))
#     rowmeta <- geneorder
#     rowannotationlist1 <- annotationlist_builder(rowmeta, customcolorlist = list(
#         log2FoldChange = colorRamp2(c(ceiling(max(rowmeta[,1])), 0, floor(min(rowmeta[,1]))), brewer.pal(3,"RdBu")),
#         pvalue = colorRamp2(seq(0,1,length = 5), brewer.pal(5,"Reds"))
#     ))
#     
#     outhm1 <- create_heatmap(counttab = plottab, scale_data = TRUE,
#                              colmetatable = plotmeta, colannotationlist = annotationlist1,
#                              rowmetatable = rowmeta, rowannotationlist = rowannotationlist1,
#                              colclusterparam = FALSE, rowclusterparam = FALSE)
#     outheatmapfile1 = paste0(outfilepathplot, pathwaylabel, ".pdf")
#     pdf(outheatmapfile1, height = 11.5, width = 10)
#     draw(outhm1$heatmap)
#     junk <- dev.off()
#     
#     outhm2 <- create_heatmap(counttab = plottab, scale_data = TRUE,
#                              colmetatable = plotmeta, colannotationlist = annotationlist1,
#                              rowmetatable = rowmeta, rowannotationlist = rowannotationlist1,
#                              colclusterparam = TRUE, rowclusterparam = FALSE)
#     outheatmapfile2 = paste0(outfilepathplot, pathwaylabel, "_clustered.pdf")
#     pdf(outheatmapfile2, height = 11.5, width = 10)
#     draw(outhm2$heatmap)
#     junk <- dev.off()
#     
# }


# Write out GCT file:
gctincounts = cbind.data.frame(id = rownames(normcounttab), normcounttab)
colmeta = cbind(sample = rownames(metatablefilt4), metatablefilt4[,1, drop=FALSE])
# rowmeta = convertedgenes[,c(2,1)]
outfilegct13 = paste0(outfilepathprocess, "normcounttab_v13.gct")
gctout13 <- create_gct(counts = gctincounts, rowmeta = NULL, colmeta = colmeta, version = "1.3", outfile = outfilegct13)
outfilegct12 = paste0(outfilepathprocess, "normcounttab_v12.gct")
gctout12 <- create_gct(counts = gctincounts, rowmeta = NULL, version = "1.2", outfile = outfilegct12)
# 
# 
# # Running ssGSEA module - uses the same counttable for each one, so the variable we actually change is the database
# genesetdblist <- list(hallmark = "/Users/tosh/Desktop/Ruggles_Lab/databases/msigdb_20191010/h.all.v7.0.symbols.gmt",
#                       kegg = "/Users/tosh/Desktop/Ruggles_Lab/databases/msigdb_20191010/c2.cp.kegg.v7.0.symbols.gmt",
#                       go = "/Users/tosh/Desktop/Ruggles_Lab/databases/msigdb_20191010/c5.all.v7.0.symbols.gmt")
# 
# ssGSEAoutlist <- list()
# for (genesetnum in seq_len(length(genesetdblist))) {
#     ## Set the gene database we are using (dont use whatever weird default)
#     gene.set.databases <- genesetdblist[[genesetnum]]
#     out.prefix <- names(genesetdblist[genesetnum])
#     
#     ## Define the out.dir and gct.file (gct.file is output from previous step)
#     out.dir <- paste0(outfilepathmaster, "ssgsea/", names(genesetdblist[genesetnum]), "/")
#     dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
#     gct.file <- outfilegct12
#     
#     ssGSEAoutlist[[genesetnum]] <- ssGSEAout <- run_ssGSEA2(gct.file = gct.file, gene.set.databases = gene.set.databases, 
#                                                             out.dir = out.dir)
#     names(ssGSEAoutlist[genesetnum]) <- names(genesetdblist)[genesetnum]
#     
#     ## Select out the data we want (the enrichment scores) and then apply a necessary but silly conversion to make numeric
#     ssGSEA_hmtab <- data.frame(ssGSEAout[,colnames(normcounttab)],row.names = ssGSEAout$id, 
#                                stringsAsFactors = FALSE, check.names = FALSE)
#     ssGSEA_hmtab[] <- lapply(ssGSEA_hmtab, function(x) as.numeric(as.character(x)))
#     
#     ## Create the heatmap
#     outheatmapfile = paste0(out.dir, names(genesetdblist[genesetnum]), "_ssgsea", ".pdf")
#     plotmeta = metatablefilt4[,1,drop=FALSE]
#     plotmeta = plotmeta[order(plotmeta[,1]),,drop=FALSE]
#     annotationlist1 = annotationlist_builder(plotmeta, 
#                                              customcolorlist = list(group = c("stress" = "red", "control" = "blue")))
#     
#     numpathways <- numpathways_plotparam
#     if (names(genesetdblist[genesetnum]) == "hallmark") {POI <- geneset_analysis_out_HALL[order(geneset_analysis_out_HALL$p.adjust),]}
#     if (names(genesetdblist[genesetnum]) == "kegg") {POI <- geneset_analysis_out_KEGG[order(geneset_analysis_out_KEGG$p.adjust),]}
#     if (names(genesetdblist[genesetnum]) == "go") {POI <- geneset_analysis_out_GO[order(geneset_analysis_out_GO$p.adjust),]}
#     POI <- POI[POI[,1] %in% rownames(ssGSEA_hmtab),][1:numpathways,]
#     if (names(genesetdblist[genesetnum]) == "go") {POI <- rbind(POI, geneset_analysis_out_GO["GO_ESTABLISHMENT_OF_CELL_POLARITY",])}
#     
#     pathwayplotparam <- POI[,1]
#     rowmeta <- data.frame(POI[,c("NES", "pvalue")], row.names = POI[,1])
#     rowannotationlist1 <- annotationlist_builder(rowmeta, customcolorlist = list(
#         NES = colorRamp2(seq(2, -2, length =3), brewer.pal(3,"RdBu")),
#         pvalue = colorRamp2(seq(0,0.1,length = 5), brewer.pal(5,"Reds"))
#     ))
#     
#     ssGSEA_outhm <- create_heatmap(counttab = ssGSEA_hmtab[pathwayplotparam,], scale_data = TRUE, 
#                                    colmetatable = plotmeta, colannotationlist = annotationlist1, 
#                                    rowmetatable = rowmeta, rowannotationlist = rowannotationlist1,
#                                    colclusterparam = FALSE, rowclusterparam = TRUE)
#     pdf(outheatmapfile, 11.5, 10)
#     draw(ssGSEA_outhm$heatmap)
#     junk <- dev.off()
#     
# }





# GOI ANALYSIS BELOW:




