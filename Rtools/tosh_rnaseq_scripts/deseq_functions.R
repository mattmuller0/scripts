################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Read in each file from output of htseq count, and concat together naming each column by the filename, then plot out some data cleaning stuff and normalize the data, then write out those tiles

## Load in Libraries
packagelist = c("ggplot2", "reshape2", "DESeq2", "gridExtra", "ggpubr", "stringr")
junk <- lapply(packagelist, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

## Cbind fill function
cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function(x) rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

# ## INFILES
# rawcounttabfile = "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia/output/rna_processing3/filtrawcounttab.txt"
# metatablefile = "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia/output/rna_processing3/metatable_filt_comp.txt"
# 
# ## Outpath
# outfilepath = "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia/output/deseq/"

# ## Read in the data
# rawcounttab = read.table(rawcounttabfile, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
# metatable = read.table(metatablefile, header = TRUE, row.names = 1, sep = ifelse(file_ext(metatablefile)=="txt", "\t", ","), stringsAsFactors = FALSE, check.names = FALSE)
# temp = as.data.frame(sapply(metatable, function(x) gsub(x = x, pattern = "\\-", replacement = "\\_")), row.names = rownames(metatable))
# metatable = temp

## DESeq works with binary comparisons - so per comp, we have to coerce our metadata into a table with the one condition column, and then a "control for" column if desired (like batch)
#####
#compcols = colnames(metatable)[grep("comp_", colnames(metatable))]

# DESeq_analysis <- function(compcols, controlcols = NULL, rawcounttab, outfilepathDEseq){
DESeq_analysis <- function(compcols, controlcols = NULL, rawcounttab){
    DEseqreslist = list()
    for (i in seq_len(length(compcols))) {
        # deseqanalysisoutfolder = paste0(outfilepathDEseq, "/", colnames(compcols)[i], "/")
        # dir.create(deseqanalysisoutfolder, recursive = TRUE, showWarnings = FALSE)

        ## Create the metadata we need by merging the current comparison and any control columns
        metasubcomp = compcols[,i,drop=FALSE]
        ## Failsafe to make sure that our treatment and control are used appropriately
        ## https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
        
        ## Select the categories for my analysis
        categories <- na.omit(unique(metasubcomp[,1]))
        ## Special circumstance for my own style tables - with 1s and 0s
        if (sum(categories %in% c(1,0)) == 2) {
            treatcat <- 1
            ctrlcat <- 0
        } else {
            ## For everything else, just select the first cat as the treatcat
            treatcat <- categories[1]
            ctrlcat <- categories[2]
        }
        
        ## Create our failsafed counttable
        treatsamps = rownames(na.omit(metasubcomp[metasubcomp[,1] == treatcat,,drop=FALSE]))
        ctrlsamps = rownames(na.omit(metasubcomp[metasubcomp[,1] == ctrlcat,,drop=FALSE]))
        treattab = rawcounttab[,colnames(rawcounttab) %in% treatsamps, drop=FALSE]
        ctrltab = rawcounttab[,colnames(rawcounttab) %in% ctrlsamps, drop=FALSE]
        ## Create the subsetted counttab that we need
        counttabDEseq = cbind(treattab, ctrltab)
        
        ## Create the subsetted metadata table
        metaDEseq = metasubcomp
        metaDEseq[metaDEseq == treatcat] = "treat"
        metaDEseq[metaDEseq == ctrlcat] = "control"
        metaDEseq[,1] = as.factor(metaDEseq[,1])
        metaDEseq = metaDEseq[colnames(counttabDEseq),,drop=FALSE]
        colnames(metaDEseq) = "DEseq_label"
        
        ## ADDING CONTROLING PARAMETERS
        if (!is.null(controlcols)) {
            ## Adding failsafe to factorize any character columns for the sake of DEseq formula
            # controlcols[sapply(controlcols, is.character)] <- lapply(controlcols[sapply(controlcols, is.character)], as.factor)
            controltab = data.frame(controlcols[match(rownames(metaDEseq), rownames(controlcols)),,drop=FALSE], stringsAsFactors = TRUE)
            metaDEseq = merge(metaDEseq, controltab, by = "row.names", sort = FALSE)
            metaDEseq = data.frame(metaDEseq[,2:ncol(metaDEseq)], row.names = metaDEseq[,1])
            colnames(metaDEseq)[2:ncol(metaDEseq)] =
                paste0(rep("control_var_", ncol(controltab)), seq(ncol(controltab)))
        }
        
        ## Fail safe checks
        dim(metaDEseq)
        dim(counttabDEseq)
        all(rownames(metaDEseq) %in% colnames(counttabDEseq))
        all(rownames(metaDEseq) == colnames(counttabDEseq))
        
        dds <- DESeqDataSetFromMatrix(countData = counttabDEseq,
                                    colData = metaDEseq,
                                    design= ~ DEseq_label)
        
        ## If we are contorlling for something - include that in the design formula
        if (!is.null(controlcols)) {
            designstring <- paste0("~",paste(colnames(metaDEseq)[2:ncol(metaDEseq)], collapse = " + "), " + DEseq_label")
            dds <- DESeqDataSetFromMatrix(countData = counttabDEseq,
                                          colData = metaDEseq,
                                          design= as.formula(designstring))
        }
        
        ## Run DEseq
        suppressMessages(dds <- DESeq(dds))
        DEseqres <- as.data.frame(results(dds))
        
        ## Save out results
        DEseqreslist[[i]] = DEseqres
        names(DEseqreslist)[[i]] <- colnames(compcols)[i]

        
    }
    return(DEseqreslist)

}

#####

# ## Output Summary Table and Genelists
# #####
# # statetype = {pvalue | padj}
# summaryparams = list(run1 = c("stattype" = "pvalue", "statcutoff" = 0.01, "log2fccutoff" = 1),
#                      run2 = c("stattype" = "padj", "statcutoff" = 0.1, "log2fccutoff" = 0),
#                      run3 = c("stattype" = "padj", "statcutoff" = 0.2, "log2fccutoff" = 0))

DESeq_summary <- function(DEseqreslist, summaryparams, outfilepath) {
    for (runnum in seq_len(length(summaryparams))) {
        ## Define the stats for the run and create the outfolder
        stattype = unname(summaryparams[[runnum]]["stattype"])
        statcutoff = as.numeric(unname(summaryparams[[runnum]]["statcutoff"]))
        log2fccutoff = as.numeric(unname(summaryparams[[runnum]]["log2fccutoff"]))
        
        summaryoutpath = paste(outfilepath, stattype, "_", statcutoff, "_log2fc_", log2fccutoff, "_DEseq_summary/", sep="")
        dir.create(summaryoutpath, recursive = TRUE, showWarnings = FALSE)
        
        #####
        deseqsubsetUP = lapply(DEseqreslist, function(x) {
            if (stattype == "pvalue") {statcol = 5}
            if (stattype == "padj") {statcol = 6}
            temp = na.omit(x[x$log2FoldChange > log2fccutoff & x[,statcol] < statcutoff,])
            temp[order(temp[,statcol], decreasing = FALSE),]})
        deseqsubsetDOWN = lapply(DEseqreslist, function(x) {
            if (stattype == "pvalue") {statcol = 5}
            if (stattype == "padj") {statcol = 6}
            temp = na.omit(x[x$log2FoldChange < -log2fccutoff & x[,statcol] < statcutoff,])
            temp[order(temp[,statcol], decreasing = FALSE),]})
        
        sumtable = cbind.data.frame(genessigUP = sapply(deseqsubsetUP, nrow), genessigDOWN = sapply(deseqsubsetDOWN, nrow), comparisons = names(DEseqreslist))
        sumplotdata = melt(sumtable, id.vars="comparisons")
        sumplotdata[,1] = gsub("_", " ", sumplotdata[,1])
        
        outplot <- ggplot(sumplotdata, aes(x=str_wrap(sumplotdata[,1], 20), y=sumplotdata[,3], fill=sumplotdata[,2])) 
        outplot <- outplot + geom_bar(stat='identity', position='dodge')
        outplot <- outplot + labs(title = paste("DGE Summary for Genes with ", stattype, " < ", statcutoff, " and |log2fc| > ", log2fccutoff, sep=""), x = "comparisons", y = "Number of Genes", fill = "")
        outplot <- outplot + scale_x_discrete(limits = str_wrap(gsub("_", " ", names(DEseqreslist)), 20))
        outplot <- outplot + theme_pubr(base_size = 10, x.text.angle = 60)
        pdf(paste0(summaryoutpath, stattype, "_", statcutoff, "_log2fc_", log2fccutoff, "_DEseq_summary.pdf"))
        print(outplot)
        junk <- dev.off()
        
        ## I also want to write out the table that goes into the summary plot - need these numbers
        comparisonorder <- gsub("_", " ", names(DEseqreslist))
        outsumplotdata <- sumplotdata[order(match(sumplotdata[,1], comparisonorder)),]
        write.table(outsumplotdata, paste0(summaryoutpath, stattype, "_", statcutoff, "_log2fc_", log2fccutoff, "_DEseq_summary_table.txt"), 
                    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
        
        ## Write out the genes that match each of those comparisons
        idx <- order(c(seq_along(deseqsubsetUP), seq_along(deseqsubsetDOWN)))
        deseqsumgenetable = do.call(cbind.fill, lapply(c(deseqsubsetUP, deseqsubsetDOWN)[idx], row.names))
        deseqsumlog2fctable = do.call(cbind.fill, lapply(c(deseqsubsetUP, deseqsubsetDOWN)[idx], function(x) x[,"log2FoldChange"]))
        deseqsumpstattable = do.call(cbind.fill, lapply(c(deseqsubsetUP, deseqsubsetDOWN)[idx], function(x) -log10(x[, stattype])))
        colnames(deseqsumgenetable) = colnames(deseqsumlog2fctable) = colnames(deseqsumpstattable) = 
            c(paste0(names(deseqsubsetUP), "_UP"), paste0(names(deseqsubsetDOWN), "_DOWN"))[idx]
        write.table(deseqsumgenetable, paste0(summaryoutpath, stattype, "_", statcutoff, "_log2fc_", log2fccutoff, "_DEseq_DEgenes_table.txt"), 
                    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na = "")
        

        
        ## Utility coding - outputting summary gene tables - table1 swapping the comparison and gene - for purpose of annotation, table2 and table3 do this swap, but instead of filling in with the comparison, it fills in with the log2fc and the pstat (either adjpval or pval, whatever is decided up front)
        templistgene = templistlog2fc = templistpstat = list()
        for (i in 1:ncol(deseqsumgenetable)) {
        templistgene[[i]] = data.frame(na.omit(deseqsumgenetable[,i]), rep(colnames(deseqsumgenetable)[i], length(na.omit(deseqsumgenetable[,i]))))
        templistlog2fc[[i]] = data.frame(na.omit(deseqsumgenetable[,i]), na.omit(deseqsumlog2fctable[,i]))
        templistpstat[[i]] = data.frame(na.omit(deseqsumgenetable[,i]), na.omit(deseqsumpstattable[,i]))
        colnames(templistgene[[i]]) = colnames(templistlog2fc[[i]]) = colnames(templistpstat[[i]]) = c("entries", colnames(deseqsumgenetable[,i,drop=FALSE]))
        }
        names(templistgene) = names(templistlog2fc) = names(templistpstat) = colnames(deseqsumgenetable)
        genecomptab = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "entries", all = TRUE), templistgene)
        log2fccomptab = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "entries", all = TRUE), templistlog2fc)
        pstatcomptab = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "entries", all = TRUE), templistpstat)
        write.table(genecomptab, paste0(summaryoutpath, stattype, "_", statcutoff, "_log2fc_", log2fccutoff, "_DEseq_comparison_annotation_table.txt"), 
                    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na = "")
        write.table(log2fccomptab, paste0(summaryoutpath, stattype, "_", statcutoff, "_log2fc_", log2fccutoff, "_DEseq_log2fc_annotation_table.txt"), 
                    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na = "")
        write.table(pstatcomptab, paste0(summaryoutpath, stattype, "_", statcutoff, "_log2fc_", log2fccutoff, "_DEseq_pstat_annotation_table.txt"), 
                    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na = "")
        
        #####
        }
}
 #####




## Draw a plot to show number of genes available at various cutoffs
# Lets do lines at diff log2fc, and then go across the pvalue on xaxis, and y-axis is number of genes
# And do it for pvalue and adjpval I guess....?
destatplot <- function(subtab) {
    statcutlevels <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2)
    log2fclevels <- c(0, 0.25, 0.5, 1, 2)
    statplotlist <- list()
    counter <- 1
    for (log2fccutnum in seq_len(length(log2fclevels))) {
        for (statcutnum in seq_len(length(statcutlevels))) {
            outval <- nrow(na.omit(subtab[abs(subtab[,1]) > log2fclevels[log2fccutnum] & subtab[,2] < statcutlevels[statcutnum],]))
            statplotlist[[counter]] <- c(log2fclevels[log2fccutnum], statcutlevels[statcutnum], outval)
            counter <- counter + 1
        }
    }
    statplottab <- do.call(rbind.data.frame, statplotlist)
    colnames(statplottab) <- c("log2fccutoff", "pstat", "numbergenes")
    statplottab[,1] <- apply(statplottab[,1,drop=FALSE], 1, as.character)
    
    pout <- linechart_plotter(indata = statplottab[,c(2,3)], groupvar = statplottab[,1,drop=FALSE], colorvar = statplottab[,1,drop=FALSE], 
                              showpoints = TRUE, labsparam = list(title = "statplot", x = "pvalue", y = "number of genes"))
    ybreaks <- c(10,50,100,200,500,1000)[c(10,50,100,200,500,1000) < max(statplottab[,3])]
    pout <- pout + geom_hline(yintercept = ybreaks, linetype = 2, alpha = 0.5)
    
    ## Reformat the out table
    outstatplottab <- dcast(statplottab, statplottab[,1] ~ statplottab[,2], value.var = "numbergenes")
    rownames(outstatplottab) <- outstatplottab[,1]
    outstatplottab <- outstatplottab[,-1]
    
    return(list(statplot = pout, stattab = outstatplottab))
    
}









# This is used to skip the upstream stuff and just add on single comparisons to the RNAseq script
addon_deseq_runs <- function(metatable_addoncomps, controlcols = NULL, outfilepathmaster, summaryparams) {
    
    # Using the project folder (outfilepathmaster) read in the files we need
    filtrawcounttab_file <- paste0(outfilepathmaster, "/rna_processing/filtrawcounttab.txt")
    filtrawcounttab <- read.table(filtrawcounttab_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    normcounttab_file <- paste0(outfilepathmaster, "rna_processing/normcounttab.txt")
    normcounttab <- read.table(normcounttab_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    
    # Create the new compcols object
    compcols <- metatable_addoncomps[,c(grepl("comp_", colnames(metatable_addoncomps))), drop=FALSE]
    metatablefilt4 <- metatable_addoncomps[,!grepl("comp_", colnames(metatable_addoncomps))]
    counttabfilt4 <- filtrawcounttab
    
    # Check and make sure our files are properly aligned
    if (all(rownames(metatable_addoncomps) %in% colnames(counttabfilt4)) + all(rownames(metatable_addoncomps) == colnames(counttabfilt4)) != 2) {
        print("Theres a problem with the metatable and counttable aligning properly")  
        break
    }
    
    ##### DESEQ ANALYSIS
    outfilepathDEseq = paste0(outfilepathmaster, "deseq/")
    
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
    # DEseqreslist <- c(DEseqreslist1, DEseqreslist2)[colnames(compcols)]
    
    ## Write out DEseq analysis
    for (DEseqrestab in seq_len(length(DEseqreslist))) {
        deseqanalysisoutfolder = paste0(outfilepathDEseq, "/", names(DEseqreslist)[DEseqrestab], "/")
        dir.create(deseqanalysisoutfolder, recursive = TRUE, showWarnings = FALSE)
        deseqanalysisoutfile = paste0(deseqanalysisoutfolder, "deseq_results_", names(DEseqreslist)[DEseqrestab], ".csv")
        write.table(DEseqreslist[[DEseqrestab]], deseqanalysisoutfile, quote = FALSE, sep = ",", row.names = TRUE, col.names=NA)
    }
    
    # Dont want to redo the summary - as that will overwrite whats already there, this is just the addons
    #DESeq_summary(DEseqreslist, summaryparams, outfilepathDEseq)
    
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
                                        colclusterparam = TRUE, rowclusterparam = TRUE)
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
    # seedparam = NULL
    # sample(1:2147483647, 1)
    seedparam = 1369634837
    
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
    
}


