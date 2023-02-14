################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Gene Set Analysis - via clusterprofiler - to be expanded upon as needed: https://yulab-smu.github.io/clusterProfiler-book/chapter1.html

## Load in Libraries
packagelist = c("clusterProfiler", "tools", "org.Hs.eg.db", "msigdbr", "GO.db", "stringdist", "enrichplot")
junk <- lapply(packagelist, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

## Load in data
#indeseqfile = "/Users/tosh/Desktop/Ruggles_Lab/projects/platelet-diabetes/output/deseq/deseq_results_comp_diabetes__DM_v_nonDM.csv"
#indeseqtable = read.table(indeseqfile, header = TRUE, row.names = 1, sep = ",")

# http://software.broadinstitute.org/gsea/msigdb/index.jsp
# H - hallmarks
# C1 - positional along each chromosome
# C2 - curated
#   CGP - chemical and genetic perturbations
#   CP - canonical
#   CP:BIOCARTA - Biocarta
#   CP:REACTOME - reactome
#   CP:KEGG
# C3 - motif
#   MIR: microRNA
#   TFT: TFs and their targets
# C4 - computational - data mined cancer microarrays
#   CGN - cancer gene neighborhoods
#   CM - cancer modules
# C5 - GO gene sets
#   BP - GO biological
#   CC - GO cellular components
#   MF - GO molecular function
# C6 - oncogenic gene sets
# C7 - immunological gene sets
# species for msigdbr: "Bos taurus" "Caenorhabditis elegans" "Canis lupus familiaris" 
# "Danio rerio" "Drosophila melanogaster" "Gallus gallus" "Homo sapiens" 
# "Mus musculus" "Rattus norvegicus" "Saccharomyces cerevisiae" "Sus scrofa"
geneset_analysis <- function(DEseqtable, rankmetric = "log2fc", pvalcutoffparam = 1, 
                             genesetparam = c("CP:KEGG"), speciesparam = "Homo sapiens", seedparam = NULL, customgeneset = NULL) {

    if (is.null(seedparam)){
        seedparamin = FALSE
    } else {
        seedparamin = TRUE
    }
    
    if (rankmetric == "log2fc") {
        gseaintable = DEseqtable[DEseqtable[,5] < pvalcutoffparam ,"log2FoldChange", drop = FALSE]
        gsearanklist = sort(setNames(object = gseaintable[,1], nm = rownames(gseaintable)), decreasing = TRUE)
    }
    if (rankmetric == "adjpval") {
        gseaintable = DEseqtable[DEseqtable[,5] < pvalcutoffparam ,c("log2FoldChange","padj"), drop = FALSE]
        gseaintable$pstat <- -log10(gseaintable$padj) * ifelse(gseaintable$log2FoldChange > 0, 1, -1)
        gsearanklist = sort(setNames(object = gseaintable[,"pstat"], nm = rownames(gseaintable)), decreasing = TRUE)
    }
    
    ## Msigdb table generation
    genesetcategories = c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
    genesetsubcategories = c("CGP", "CP", "CP:BIOCARTA", "CP:REACTOME", "CP:KEGG", 
                             "MIR", "TFT", "CGN", "CM", "BP", "CC", "MF")
    genesetcategories_select = genesetcategories[genesetcategories %in% genesetparam]
    genesetsubcategories_select = genesetsubcategories[genesetsubcategories %in% genesetparam]
    
    if (length(genesetcategories_select) > 0) {
        m_t2g_cat <- as.data.frame(msigdbr(species = speciesparam, 
            category = c(genesetcategories_select))[,c("gs_name", "gene_symbol")])
    } else {m_t2g_cat = NULL}
    if (length(genesetsubcategories_select) > 0) {
        m_t2g_subcat <- as.data.frame(msigdbr(species = speciesparam, 
            subcategory = c(genesetsubcategories_select))[,c("gs_name", "gene_symbol")])
    } else {m_t2g_subcat = NULL}
    m_t2g = rbind(m_t2g_cat, m_t2g_subcat)
    
    ## If you want to test a single custom geneset - input here to override the null argument, and will remake the m_t2g object
    ## Custom object needs to be 2 columns, with the genesets in the first col, and the genes in the second col
    # gs_name gene_symbol
    # 1 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS     ALDH1L1
    # 2 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS     ALDH1L2
    if (!is.null(customgeneset)) {
        colnames(customgeneset) <- c("gs_name", "gene_symbol")
        m_t2g <- rbind(customgeneset, m_t2g)
    }
    
    gseaout = GSEA(geneList = gsearanklist, TERM2GENE = m_t2g, by = "fgsea",
                   minGSSize = 1, maxGSSize = 1000, 
                   # nPerm = 10000, ## Apparently they dont want us to use nPerm anymore, so removing this! 2021-08-18
                   seed = seedparamin,
                   pvalueCutoff = 1.1, verbose = TRUE)
    return(gseaout)
}


hypergeo_genetest <- function(DEseqtable, statcutoffparam = c("stattype" = "pvalue", "pstatcutoff" = 0.01, "log2fccutoff" = 0), 
                              genesetparam = c("C5"), speciesparam = "Homo sapiens",
                              customgeneset = NULL) {
    
    ## This allows you to just input a list of genes and override the input for a DEseq table
    if (length(DEseqtable) == 1){
        geotestUPlist <- geotestDOWNlist <- DEseqtable[,1]
    } else {
        geotestUPlist = rownames(DEseqtable[
            DEseqtable[,statcutoffparam["stattype"]] < as.numeric(statcutoffparam["pstatcutoff"]) &
            !is.na(DEseqtable[,statcutoffparam["stattype"]]) &
            DEseqtable[,"log2FoldChange"] > as.numeric(statcutoffparam["log2fccutoff"]),, drop = FALSE])
        geotestDOWNlist = rownames(DEseqtable[
            DEseqtable[,statcutoffparam["stattype"]] < as.numeric(statcutoffparam["pstatcutoff"]) &
            !is.na(DEseqtable[,statcutoffparam["stattype"]]) &
            DEseqtable[,"log2FoldChange"] < -as.numeric(statcutoffparam["log2fccutoff"]),, drop = FALSE])
    }
    
    ## Msigdb table generation
    genesetcategories = c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
    genesetsubcategories = c("CGP", "CP", "CP:BIOCARTA", "CP:REACTOME", "CP:KEGG", 
                             "MIR", "TFT", "CGN", "CM", "BP", "CC", "MF")
    genesetcategories_select = genesetcategories[genesetcategories %in% genesetparam]
    genesetsubcategories_select = genesetsubcategories[genesetsubcategories %in% genesetparam]
    
    if (length(genesetcategories_select) > 0) {
        m_t2g_cat <- as.data.frame(msigdbr(species = speciesparam, 
                                           category = c(genesetcategories_select))[,c("gs_name", "gene_symbol")])
    } else {m_t2g_cat = NULL}
    if (length(genesetsubcategories_select) > 0) {
        m_t2g_subcat <- as.data.frame(msigdbr(species = speciesparam, 
                                              subcategory = c(genesetsubcategories_select))[,c("gs_name", "gene_symbol")])
    } else {m_t2g_subcat = NULL}
    m_t2g = rbind(m_t2g_cat, m_t2g_subcat)
    
    ## If you want to test a single custom geneset - input here to override the null argument, and will remake the m_t2g object
    ## Custom object needs to be 2 columns, with the genesets in the first col, and the genes in the second col
    # gs_name gene_symbol
    # 1 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS     ALDH1L1
    # 2 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS     ALDH1L2
    if (!is.null(customgeneset)) {
        colnames(customgeneset) <- c("gs_name", "gene_symbol")
        m_t2g <- rbind(customgeneset, m_t2g)
    }
    
    enricherUPout = enricher(gene = geotestUPlist, qvalueCutoff = 2, pvalueCutoff = 2,
                        minGSSize = 0, maxGSSize = 10000, TERM2GENE = m_t2g)
    enricherDOWNout = enricher(gene = geotestDOWNlist, qvalueCutoff = 2, pvalueCutoff = 2,
                        minGSSize = 0, maxGSSize = 10000, TERM2GENE = m_t2g)
    if (is.null(enricherUPout) & is.null(enricherDOWNout)) {
        enricherUPout <- enricherDOWNout <- data.frame(matrix(nrow = 0, ncol = 9, 
                                                              dimnames = list(NULL, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"))))
    }
    return(list(enricherUPout = enricherUPout, enricherDOWNout = enricherDOWNout))
}


gsea_barplot <- function(gseaout, pstatparam, numterms = 10, titleparam, clean_genesetIDs = FALSE){
    if (pstatparam == "pvalue"){pstat = 6}
    if (pstatparam == "p.adjust"){pstat = 7}
    
    gseaplotin = data.frame(gseaout)[order(abs(gseaout[,pstatparam]), decreasing = FALSE)[1:numterms],c(1,5,pstat)]
    gseaplotin = gseaplotin[order(gseaplotin[,"NES"], decreasing = TRUE),]
    if (clean_genesetIDs == TRUE) {
        gseaplotin[,1] <- unlist(lapply(gsub("_", " ", gsub("GOBP_|GOMF_|GOCC_", "", gseaplotin[,1])), 
                                        function(x) simpleCap(tolower(x))))
    } else {
        gseaplotin[,1] = gsub("_", " ", gseaplotin[,1])
    }
    
    #pout <- ggplot(gseaplotin, mapping = aes(x = gseaplotin[,1], y = gseaplotin[,2], fill = gseaplotin[,3]))
    pout <- ggplot(gseaplotin, mapping = aes(x = str_wrap(gseaplotin[,1], 20), y = gseaplotin[,2], fill = gseaplotin[,3]))
    pout <- pout + geom_bar(stat = "identity")
    # pout <- pout + scale_fill_gradient(low = "palegreen", high = "darkgreen", limits=c(0,0.1))
    pout <- pout + scale_fill_viridis(limits = c(0,max(0.1, max(gseaplotin[,3]))), direction = -1, option = "magma")
    pout <- pout + scale_x_discrete(limits=rev(str_wrap(gseaplotin[,1], 20)))
    pout <- pout + coord_flip()
    pout <- pout + labs(x="Geneset", y = "NES", fill = pstatparam, title = titleparam)
    pout
    return(pout)
}


## Takes output from gsea_analysis function above ^^
gsea_custom_randomwalkplot <- function(gseaobject, geneSetID, addNESvalue = TRUE) {
    gseaplot_out <- gseaplot(gseaobject, geneSetID = geneSetID, title = geneSetID)
    if (addNESvalue) {
        # Place text very specifically
        if (sign(gseaobject[geneSetID,"enrichmentScore"]) < 0) {
            xcor <- 0.16
            ycor <- 0.16
        } else {
            xcor <- 0.84
            ycor <- 0.32
        }
        gseaplot_out <- gseaplot_out + annotate("text", x = xcor, y = ycor,
                            label = paste0("NES: ", round(gseaobject[geneSetID,"NES"], 3), "\n",
                                           "p.adjust: ", round(gseaobject[geneSetID,"p.adjust"],4)
                            ))
    }
    return(gseaplot_out)
}



## Takes output from gsea_analysis function above ^^
gsea_custom_randomwalkplot_v2 <- function(gseaobject, geneSetID, addNESvalue = TRUE) {
    gseaplot_out <- gseaplot2(x = gseaobject, geneSetID = geneSetID, title = geneSetID, base_size = 16, )
    if (addNESvalue) {
        # Place text very specifically
        if (sign(gseaobject[geneSetID,"enrichmentScore"]) < 0) {
            xcor <- 0.16
            ycor <- 0.16
        } else {
            xcor <- 0.84
            ycor <- 0.32
        }
        gseaplot_out <- gseaplot_out + annotate("text", x = xcor, y = ycor,
                                        label = paste0("NES: ", round(gseaobject[geneSetID,"NES"], 3), "\n",
                                                       "p.adjust: ", round(gseaobject[geneSetID,"p.adjust"],4)
                                        ))
        gseaplot_out + geom_vline(xintercept = 0.9, linetype = 2, color = "red")
    }
    return(gseaplot_out)
}


gene_in_geneset_heatmap <- function(plotintable, gene_delineator = "/", clean_genesetIDs = FALSE, clusterparam = FALSE) {
    melted_ID_v_gene_table <- do.call(rbind, apply(plotintable, 1, function(x) {
        genes <- unlist(strsplit(x[2], split = "/"))
        cbind(ID = x[1], genes)
    }))
    ID_v_gene_matrix <- dcast(data.frame(melted_ID_v_gene_table), ID ~ genes, fun.aggregate = length)
    rownames(ID_v_gene_matrix) <- ID_v_gene_matrix[,"ID"]
    
    if (clean_genesetIDs == TRUE) {
        rownames(ID_v_gene_matrix) <- unlist(lapply(gsub("_", " ", gsub("GOBP_|GOMF_|GOCC_", "", ID_v_gene_matrix[,"ID"])), 
                                                    function(x) simpleCap(tolower(x))))
    }
    
    create_heatmap(counttab = ID_v_gene_matrix[,!grepl("ID", colnames(ID_v_gene_matrix))], subsetnum = FALSE, scale_data = FALSE,
                   colclusterparam = clusterparam, rowclusterparam = clusterparam, 
                   heatmapcolorparam = colorRamp2(breaks = c(0,1), colors = c("white", "black")), addborders = TRUE)
    
}





## Convert goterms to go ids
goterm_to_goID <- function(ingoterms) {
    ## First read in our GO.db list of all go terms
    allgoterms <- Term(GOTERM)
    cleaned_allgoterms <- tolower(gsub("-|\\/", " ", allgoterms))
    
    ## Then clean our inputted list of goterms
    cleaned_ingoterms <- gsub("_", " ", tolower(gsub("GO_", "", ingoterms)))
    
    ## Then find which we have, which we done, and return all
    # Separated for clarity and failsafing
    p1 <- data.frame(cleaned_GOterm = cleaned_allgoterms[cleaned_allgoterms %in% cleaned_ingoterms],
                     GOID = names(cleaned_allgoterms[cleaned_allgoterms %in% cleaned_ingoterms]), stringsAsFactors = FALSE)
    p2 <- data.frame(cleaned_GOterm = cleaned_ingoterms,
                     Input_GOterm = ingoterms)
    foundIDs <- merge(p1, p2, all.x = TRUE, by = "cleaned_GOterm", sort = FALSE)
    
    ## For the not found - lets find a closest match
    notfoundIDs <- cleaned_ingoterms[!cleaned_ingoterms %in% cleaned_allgoterms]
    
    ## Can we return a closest go term for each one of these? YES
    ClosestMatch2 = function(string, stringVector){stringVector[amatch(string, stringVector, maxDist=Inf)]}
    closematchtab1 <- data.frame(sapply(notfoundIDs, function(x) ClosestMatch2(x, cleaned_allgoterms)))
    closematchtab1[,c(2,3)] <- do.call(rbind, strsplit(rownames(closematchtab1), split = "\\."))
    colnames(closematchtab1) <- c("ClosestMatch_GOterm", "cleaned_GOterm", "GOID")
    closematchtab1 <- closematchtab1[,c("GOID", "ClosestMatch_GOterm", "cleaned_GOterm")]
    rownames(closematchtab1) <- NULL
    # Reattach the original input term
    closematchtab2 <- cbind(closematchtab1, Input_GOterm = ingoterms[!cleaned_ingoterms %in% cleaned_allgoterms])
    
    return(list(converted_goterm_to_goID = foundIDs, missing_goterm_to_id = closematchtab2))
}


## Import a vecotr of goIDS - return a table of the ancestors
find_goterm_ancestors <- function(ingoids) {
    ## Define the function we need
    getAllBPAncestors <- function(goids, GOBPANCESTOR_inlist = NULL){
        if (is.null(GOBPANCESTOR_inlist)) {GOBPANCESTOR_inlist = as.list(GOBPANCESTOR)}
        ## Apparently not all goids HAVE ancestors. So.... we need an if statement to only run this if they have an ancestor, if not, return itself?
        if (goids %in% names(GOBPANCESTOR_inlist)) {
            ans <- unique(unlist(mget(goids, GOBPANCESTOR), use.names=FALSE))
            ans <- c(ans[!is.na(ans)], goids)
        } else {
            ans <- goids
        }
        return(ans)
    }
    
    ## Then apply over our ingoid vector this function to find our ancestors
    # Need to predefine out BPAncestors object because thats the major time crunch to do each time
    GOBPANCESTOR_inlist = as.list(GOBPANCESTOR)
    anc_out1 <- lapply(ingoids, function(x) getAllBPAncestors(x, GOBPANCESTOR_inlist))
    ancestor_count <- table(unlist(anc_out1))
    ancestor_count <- data.frame(ancestor_count[order(ancestor_count, decreasing = TRUE)])
    rownames(ancestor_count) <- ancestor_count[,1]
    
    # with the ancestor count table - see what GO terms we have and attach as annotation
    allgoterms <- data.frame(Term(GOTERM), stringsAsFactors = FALSE)
    ancestor_counttab_annot <- merge(ancestor_count, allgoterms, by = "row.names", all.x = TRUE, sort = FALSE)
    
    outtab <- ancestor_counttab_annot[,c(2:4)]
    colnames(outtab) <- c("GOID", "Number_of_Mentions", "GOTerm")
    
    # Return this table
    return(ancestral_GOterm_counttable = outtab, ancesterlist = anc_out1)
    
}

# getAllBPChildren <- function(goids)
# {
#     ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names=FALSE))
#     ans <- ans[!is.na(ans)]
# }
# getAllBPAncestors <- function(goids)
# {
#     ans <- unique(unlist(mget(goids, GOBPANCESTOR), use.names=FALSE))
#     ans <- ans[!is.na(ans)]
# }


# ## Classic Random Walk Scoring table for GSEA Plot
# library(DOSE)
# data(geneList)
# de <- names(geneList)[abs(geneList) > 2]
# edo2 <- gseDO(geneList)
# p1 <- gseaplot(edo2, geneSetID = 1, by = "runningScore", title = edo2$Description[1])
# p2 <- gseaplot(edo2, geneSetID = 1, by = "preranked", title = edo2$Description[1])
# p3 <- gseaplot(edo2, geneSetID = 1, title = edo2$Description[1])
















# http://software.broadinstitute.org/gsea/msigdb/index.jsp
# H - hallmarks
# C1 - positional along each chromosome
# C2 - curated
#   CGP - chemical and genetic perturbations
#   CP - canonical
#   CP:BIOCARTA - Biocarta
#   CP:REACTOME - reactome
#   CP:KEGG
# C3 - motif
#   MIR: microRNA
#   TFT: TFs and their targets
# C4 - computational - data mined cancer microarrays
#   CGN - cancer gene neighborhoods
#   CM - cancer modules
# C5 - GO gene sets
#   BP - GO biological
#   CC - GO cellular components
#   MF - GO molecular function
# C6 - oncogenic gene sets
# C7 - immunological gene sets
# species for msigdbr: "Bos taurus" "Caenorhabditis elegans" "Canis lupus familiaris" 
# "Danio rerio" "Drosophila melanogaster" "Gallus gallus" "Homo sapiens" 
# "Mus musculus" "Rattus norvegicus" "Saccharomyces cerevisiae" "Sus scrofa"
# geneset_heatmap <- function(DEseqtable, pvalcutoffparam = 1, genesetparam = c("CP:KEGG"), speciesparam = "Homo sapiens", 
#                             pathwayparam = "KEGG_PENTOSE_PHOSPHATE_PATHWAY") {
#     
#     gseaintable = DEseqtable[DEseqtable[,5] < pvalcutoffparam ,c(2), drop = FALSE]
#     gsearanklist = sort(setNames(object = gseaintable[,1], nm = rownames(gseaintable)), decreasing = TRUE)
#     
#     ## Msigdb table generation
#     genesetcategories = c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
#     genesetsubcategories = c("CGP", "CP", "CP:BIOCARTA", "CP:REACTOME", "CP:KEGG", 
#                              "MIR", "TFT", "CGN", "CM", "BP", "CC", "MF")
#     genesetcategories_select = genesetcategories[genesetcategories %in% genesetparam]
#     genesetsubcategories_select = genesetsubcategories[genesetsubcategories %in% genesetparam]
#     
#     if (length(genesetcategories_select) > 0) {
#         m_t2g_cat <- as.data.frame(msigdbr(species = speciesparam, 
#                                            category = c(genesetcategories_select))[,c("gs_name", "gene_symbol")])
#     } else {m_t2g_cat = NULL}
#     if (length(genesetsubcategories_select) > 0) {
#         m_t2g_subcat <- as.data.frame(msigdbr(species = speciesparam, 
#                                               subcategory = c(genesetsubcategories_select))[,c("gs_name", "gene_symbol")])
#     } else {m_t2g_subcat = NULL}
#     m_t2g = rbind(m_t2g_cat, m_t2g_subcat)
# 
#     
#     as.data.frame(msigdbr(species = "Mus musculus", subcategory = c("CP:KEGG"))[,c("gs_name", "gene_symbol")])
#     
# }


















# dotplot(edo, showCategory=10)
# 
# 
# cnetplot(gseaout, showCategory = c("KEGG_CELL_ADHESION_MOLECULES_CAMS",
#                                    "KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS",
#                                    "KEGG_O_GLYCAN_BIOSYNTHESIS"), foldChange = gsearanklist, categorySize="pvalue")



### GOI HEATMAP FOR A SINGLE GENE SET
# dbtab = as.data.frame(msigdbr(species = "Homo sapiens", subcategory = c("CP:KEGG")))[,c("gs_name", "gene_symbol")]
# genesetOIcounttab = data.frame(KEGG_CELL_ADHESION_MOLECULES_CAMS = 
#                     gsearanklist[dbtab[dbtab[,1] == "KEGG_CELL_ADHESION_MOLECULES_CAMS",2]], 
#                     row.names = dbtab[dbtab[,1] == "KEGG_CELL_ADHESION_MOLECULES_CAMS",2])
# genesetOIcounttab = genesetOIcounttab[order(genesetOIcounttab[,1], decreasing = TRUE),,drop=FALSE]
# temppdfoutfile = "/Users/tosh/Desktop/Ruggles_Lab/projects/platelet-cholesterol-human/output/test.pdf"
# genesetOI = create_heatmap(counttab = genesetOIcounttab, pdfoutfile = temppdfoutfile)



# 
# heatplot(geneset_analysis_out_HALL, foldChange = gsearanklist)
# 
# emapplot(geneset_analysis_out_HALL, showCategory = 10)
# 
# gseaplot(geneset_analysis_out_HALL, geneSetID = "HALLMARK_COAGULATION", by="runningScore")
# gseaplot2(geneset_analysis_out_HALL, geneSetID = "HALLMARK_COAGULATION")
# 
# 
# 
# library("pathview")
# pathview(gene.data  = gsearanklist,
#          pathway.id = "hsa04110",
#          species    = "hsa",
#          limit      = list(gene=max(abs(geneList)), cpd=1))

# filename=system.file("extdata/gse16873.demo", package = "pathview")
# gse16873=read.delim(filename, row.names=1)
# gse16873.d=gse16873[,2*(1:6)]-gse16873[,2*(1:6)-1]
# 
# i <- 1
# pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
#                      species = "hsa", out.suffix = "gse16873", kegg.native = T)
# list.files(pattern="hsa04110", full.names=T)
# 



######### POSSIBLE PLOTTING FUNCTIONS
# library(DOSE)
# data(geneList)
# de <- names(geneList)[abs(geneList) > 2]
# 
# edo <- enrichDGN(de)
# library(enrichplot)
# barplot(edo, showCategory=20)
# barplot(geneset_analysis_out_GO, showCategory=5)
# barplot(gseaout)
# dotplot(gseaout, color = "pvalue", showCategory = 10, x = "NES", size = "setSize")
# heatplot(gseaout, foldChange = gsearanklist)
# 
# barplot(edo)



### THIS ALLOWS FOR (NOT) - downloaded in 2018 LIVE CURATION OF MSIGDB GENESETS
#library(msigdbr)
# msigdbsets <- as.data.frame(msigdbr(species = "Homo sapiens"))
# #t1 = msigdbsets[msigdbsets[,3] == "C2",]
# #kegg_genesets = msigdbsets[msigdbsets[,4] == "CP:KEGG",]
# 
# m_t2g <- msigdbr(species = "Homo sapiens", subcategory = "CP:KEGG") %>% 
#     dplyr::select(gs_name, gene_symbol)

#### WE NEED A NAMED DECREASING VECTOR - WITH THE VALUES AS LOG2FC AND NAMES AS GENES (according to vignette)
#### NOTE THAT YOU CAN DEFINE GENES AS DE FIRST THOUGH (vignette usees logfc of 2)

## Create our gene list

#gseagenelist = names(gsearanklist)



# m_t2g <- msigdbr(species = "Homo sapiens", subcategory = "CP:KEGG") %>% 
#     dplyr::select(gs_name, gene_symbol)
# 
# 
# 
# 
# 
# data(geneList, package="DOSE")
# gene <- names(geneList)[abs(geneList) > 2]
# 
# gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
# c5 <- read.gmt(gmtfile)
# 
# egmt <- enricher(gene, TERM2GENE=c5)
# egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)



## ARCHIVE
#####
# genesetfile1 = "/Users/tosh/Desktop/Ruggles_Lab/databases/msigdb_20190711/c2.cp.kegg.v6.2.symbols.gmt"
# intablefile = "/Users/tosh/Desktop/Ruggles_Lab/projects/newman-pace/data/count_normalized_formatted_filt005.txt"
# indeseqfile = "/Users/tosh/Desktop/Ruggles_Lab/projects/newman-pace/data/dge.PACE-vs-Healthy.csv"
# ingenelistfile = "/Users/tosh/Desktop/Ruggles_Lab/projects/newman-pace/data/dge.PACE-vs-Healthy_GOI_qval01.txt"
# 
# intable = read.table(intablefile, header = TRUE, row.names = 1, sep = ifelse(file_ext(intablefile)=="txt", "\t", ","), stringsAsFactors = FALSE, check.names = FALSE)
# 
# indeseq = read.table(indeseqfile, header = TRUE, row.names = 1, sep = ifelse(file_ext(indeseqfile)=="txt", "\t", ","), stringsAsFactors = FALSE, check.names = FALSE)
# deseqGOI = sort(setNames(indeseq$log2FoldChange, as.character(row.names(indeseq))), decreasing=TRUE)
# 
# ingenelist = read.table(ingenelistfile, header = FALSE, sep = ifelse(file_ext(ingenelistfile)=="txt", "\t", ","), stringsAsFactors = FALSE, check.names = FALSE)
# genelistGOI = ingenelist[,1]
# 
# ## Read in GMT file
# gs1 = read.gmt(genesetfile1)
# 
# 
# wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
# wp2gene <- read.gmt(wpgmtfile)
# wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
# wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
# wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
# ## Looks like TERM2GENE - is just 2 columns, melted table with the first column being genesets (terms) and the second column being genes in the geneset
# ## TERM2NAME - is the conversion of "term" to the "name" - where the name is the geneset(redunant if youre geneset is already named)
# 
# test1 = enricher(gene = genelistGOI, TERM2GENE = gs1, qvalueCutoff = 1, pvalueCutoff = 1)
# test2 = GSEA(geneList = deseqGOI, TERM2GENE = gs1, pvalueCutoff = 1, verbose = TRUE)
# 
# #ewp <- setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID") ## Code snipper to convert entrez to genesymbol
# 
# ### CELL MARKERS
# cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
#     tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
#     dplyr::select(cellMarker, geneID) %>%
#     dplyr::mutate(geneID = strsplit(geneID, ', '))
# 
# test1 = enricher(gene = genelistGOI, TERM2GENE = cell_markers, qvalueCutoff = 1, pvalueCutoff = 1)