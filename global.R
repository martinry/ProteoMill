library(igraph)
library(visNetwork)

#setwd("C://Users/martinry/qodb-shiny/")
setwd("~/qodb-shiny/")
# 
# setClass("Experiment", representation(
#     expData         = "data.table",
#     expDataSrc      = "data.table",
#     annotationData  = "data.table",
#     aliases         = "data.table",
#     sampleData      = "data.frame",
#     contrast        = "data.table",
#     pathwayData     = "data.table")
# )

# Upload dataset ----

undup <- function(genes){
    genes[!is.na(genes)]
}

upload_data <- function(path, sep, i){
    
    data_wide <- data.table::fread(
        path,
        sep = sep,
        dec = ".",
        header = T)
    
    data_wide <- data_wide[!duplicated(names(data_wide)[1])]
    #data_wide <- data_wide[apply(data_wide[, 2:ncol(data_wide)], 1, function(x) sum(x, na.rm = T) > 500)]
    
    for(j in seq_along(data_wide)){
        set(data_wide, i = which(data_wide[[j]] == 0 & is.numeric(data_wide[[j]])), j = j, value = NA)
    }
    
    empty_rows <- apply(data_wide[,2:ncol(data_wide)], 1, function(x) all(is.na(x)))
    data_wide <- data_wide[!empty_rows,]
    
    assign('tmpData', data_wide, envir = .GlobalEnv)
    
    data_origin <- data_wide
    
    convertColumns <- c("UNIPROTID", "ENTREZID", "SYMBOL", "GENEID", "PROTEINID")
    
    assign('convertColumns', convertColumns, envir = .GlobalEnv)
    
    keys <- data_wide[, as.character(.SD[[1L]])][1:10]
    
    if(i == "auto") {
        
        tr1 <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = keys, column = "SYMBOL", keytype = "UNIPROTID", multiVals = "first")
        tr2 <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = keys, column = "UNIPROTID", keytype = "ENTREZID", multiVals = "first")
        tr3 <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = keys, column = "UNIPROTID", keytype = "SYMBOL", multiVals = "first")
        tr4 <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = keys, column = "UNIPROTID", keytype = "GENEID", multiVals = "first")
        tr5 <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = keys, column = "UNIPROTID", keytype = "PROTEINID", multiVals = "first")
        
        trs <- list(tr1[!is.na(tr1)],
                    tr2[!is.na(tr2)],
                    tr3[!is.na(tr3)],
                    tr4[!is.na(tr4)],
                    tr5[!is.na(tr5)])
        
#        tr <- as.data.table(trs[which.max(lapply(lapply(trs, lengths), sum))])
        
        i <- convertColumns[which.max(lapply(lapply(trs, lengths), sum))]
        
        if(i == "UNIPROTID"){
            
            tr_all <- data.table(AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = data_wide[, as.character(.SD[[1L]])], columns = convertColumns, keytype = i))
            tr_all <- tr_all[!duplicated(UNIPROTID)]
            
        } else {
            tr <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = data_wide[, as.character(.SD[[1L]])], column = "UNIPROTID", keytype = i, multiVals = "first")
            tr <- tr[!is.na(tr)]
            tr <- tr[!duplicated(tr)]
            
            tr <- data.table(i = names(tr), "UNIPROTID" = unname(tr))
            
            tr_all <- data.table(AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = data_wide[, as.character(.SD[[1L]])], columns = convertColumns, keytype = i))
            tr_all <- tr_all[UNIPROTID %in% tr$UNIPROTID]
            tr_all <- tr_all[!duplicated(UNIPROTID)]
        }
        
        

        
    } else {
        i <- convertColumns[which.max(lapply(lapply(trs, lengths), sum))]
        
        if(i == "UNIPROTID"){
            
            tr_all <- data.table(AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = data_wide[, as.character(.SD[[1L]])], columns = convertColumns, keytype = i))
            tr_all <- tr_all[!duplicated(UNIPROTID)]
            
        } else {
            tr <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = data_wide[, as.character(.SD[[1L]])], column = "UNIPROTID", keytype = i, multiVals = "first")
            tr <- tr[!is.na(tr)]
            tr <- tr[!duplicated(tr)]
            
            tr <- data.table(i = names(tr), "UNIPROTID" = unname(tr))
            
            tr_all <- data.table(AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = data_wide[, as.character(.SD[[1L]])], columns = convertColumns, keytype = i))
            tr_all <- tr_all[UNIPROTID %in% tr$UNIPROTID]
            tr_all <- tr_all[!duplicated(UNIPROTID)]
        }
    }
    
    names(data_wide)[1] <- i
    
    setkeyv(tr_all, i)
    
    setkeyv(data_wide, i)
    
    tr_all <- tr_all[data_wide[, ..i], on = i]
    
    data_wide <- data_wide[tr_all, nomatch = 0]
    
    data_wide$ENTREZID <- as.character(data_wide$ENTREZID)
    
    for(j in seq_along(data_wide)){
        set(data_wide, i = which(duplicated(data_wide[[j]]) & is.character(data_wide[[j]])), j = j, value = NA)
    }
    
    data_wide[is.na(ENTREZID), ENTREZID := paste0("MISSING_", seq(1:length(is.na(ENTREZID))))]
    data_wide[is.na(SYMBOL), SYMBOL := paste0("MISSING_", seq(1:length(is.na(SYMBOL))))]
    data_wide[is.na(UNIPROTID), UNIPROTID := paste0("MISSING_", seq(1:length(is.na(UNIPROTID))))]
    data_wide[is.na(GENEID), GENEID := paste0("MISSING_", seq(1:length(is.na(GENEID))))]
    data_wide[is.na(PROTEINID), PROTEINID := paste0("MISSING_", seq(1:length(is.na(PROTEINID))))]
    
    setcolorder(data_wide, c(convertColumns, names(data_origin[,2:ncol(data_origin)])))
    
    assign('data_wide', data_wide, envir = .GlobalEnv)
    assign('data_origin', data_wide, envir = .GlobalEnv)
    
    sample_data(data_wide)
    
    
    # 
    # 
    # 
    # 
    # tr <- tr[!is.na(UNIPROTID)]
    # tr <- tr[!duplicated(get(i))]
    # 
    # 
    # setkeyv(tr, i)
    # names(data_wide)[1] <- i
    # setkeyv(data_wide, i)
    # data_wide <- data_wide[tr, nomatch = 0]
    # 
    # 
    # dt <- data_wide[, lapply(.SD, sum), by = "UNIPROTID", .SDcols = c(2:(ncol(data_wide)-4))]
    # setkeyv(dt, "UNIPROTID")
    # 
    # setkeyv(tr, "UNIPROTID")
    # data_wide <- dt[tr, nomatch = 0]
    # data_wide <- data_wide[!duplicated(UNIPROTID)]
    # 
    
    # for(j in seq_along(data_wide)){
    #     set(data_wide, i = which(data_wide[[j]] == 0 & is.numeric(data_wide[[j]])), j = j, value = NA)
    # }
    # 
    # empty_rows <- apply(data_wide[, -..convertColumns], 1, function(x) all(is.na(x)))
    # data_wide <- data_wide[!empty_rows,]
    
    

    
    # setkeyv(tr, i)
    # setkeyv(data_wide, names(data_wide)[1])
    # data_wide <- data_wide[tr, nomatch = 0]
    # dt <- data_wide[, lapply(.SD, sum), by="UNIPROTID", .SDcols = c(2:(ncol(data_wide)-5))]
    # setkey(dt, "UNIPROTID")
    # setkeyv(tr, "UNIPROTID")
    # data_wide <- dt[tr]
    # 
    # #n <- names(data_wide[, ncol(data_wide)-4:0, with = F])
    # 
    # #userIDColumn <- convertColumns[!convertColumns %in% n]
    # 
    # colnames(data_wide) <- c(i, colnames(data_wide[,2:ncol(data_wide)]))
    # 
    # setcolorder(data_wide, c(convertColumns, names(data_origin[,2:ncol(data_origin)])))
    # 
    # data_wide <- data_wide[!duplicated(data_wide[, c(3, 7:ncol(data_wide)), with = F])]
    # 
    # found_ids <- data_wide[, as.character(.SD[[i]])]
    # 
    # not_found_ids <- tmpData[!(colnames(data_origin)[1] %in% found_ids)]
    # not_found_ids[,i] <- NA
    # colnames(not_found_ids)[1] <- i
    # setcolorder(not_found_ids, c(convertColumns, names(data_origin[,2:ncol(data_origin)])))
    # 
    # data_wide <- rbindlist(list(data_wide, not_found_ids))
    # 
    # empty_rows <- apply(data_wide[, -..convertColumns], 1, function(x) all(is.na(x)))
    # data_wide <- data_wide[!empty_rows,]
    # 
    # data_wide$ENTREZID <- as.character(data_wide$ENTREZID)
    # 
    # data_wide[is.na(ENTREZID), ENTREZID := paste0("MISSING_", seq(1:length(is.na(ENTREZID))))]
    # data_wide[is.na(SYMBOL), SYMBOL := paste0("MISSING_", seq(1:length(is.na(SYMBOL))))]
    # data_wide[is.na(UNIPROTID), UNIPROTID := paste0("MISSING_", seq(1:length(is.na(UNIPROTID))))]
    # data_wide[is.na(GENEID), GENEID := paste0("MISSING_", seq(1:length(is.na(GENEID))))]
    # data_wide[is.na(TXID), TXID := paste0("MISSING_", seq(1:length(is.na(TXID))))]
    # data_wide[is.na(PROTEINID), PROTEINID := paste0("MISSING_", seq(1:length(is.na(PROTEINID))))]
    # 
    # assign('data_wide', data_wide, envir = .GlobalEnv)
    # assign('data_origin', data_wide, envir = .GlobalEnv)
    
    #sample_data(data_wide)
}


# Interaction data ----

#if(!exists("interactions")){
# 
# interactions <- data.table::fread("C://Users/martinry/interactions6.txt")
# assign("interactions", interactions, envir = .GlobalEnv)
#}

# if(!exists("uniprot_to_string_src")){
#     
#     uniprot_to_string_src <- knee::collect('https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz')
#     uniprot_to_string_src <- data.table::fread(uniprot_to_string_src)
#     
#     uniprot_to_string_src$up <- gsub("\\|.*", "", uniprot_to_string_src$V2)
#     
#     uniprot_to_string_src <- uniprot_to_string_src[, c(3,6)]
#     
#     assign("uniprot_to_string_src", uniprot_to_string_src, envir = .GlobalEnv)
# }
# 
# if(!exists("interactions")){
#     if(!file.exists(file.path(system.file(package = 'knee'), "data", "9606.protein.links.v11.0.txt"))){
#         interactions <- knee:::get_interactions()
#     } else {
#         
#         interactions <- data.table::fread(file.path(system.file(package = 'knee'), "data", "9606.protein.links.v11.0.txt"),
#                                           sep = ' ',
#                                           header = T)
#     }
#     assign("interactions", interactions, envir = .GlobalEnv)
    
    # 
    # setkey(interactions, protein1)
    # setkey(uniprot_to_string_src, V3)
    # 
    # interactions2 <- interactions[uniprot_to_string_src, nomatch = 0]
    # 
    # setkey(interactions2, protein2)
    # 
    # interactions3 <- interactions2[uniprot_to_string_src, nomatch = 0]
    # 
    # keycols = c("protein1","protein2")
    # setkeyv(interactions3, keycols)
    # 
    # 
    # keycols = c("item_id_a","item_id_b")
    # setkeyv(actions, keycols)
    # 
    # interactions4 <- interactions3[actions, nomatch = 0]
    # 
    # fwrite(interactions4)
    
    # interactions5 <- interactions4[, c(4:10)]
    # 
    # colnames(interactions5) <- c("protein1", "protein2", "mode", "action", "is_directional", "a_is_acting", "score")
    # 
    # interactions5$score <- interactions5$score / 100
    # 
    # fwrite(interactions5, "interactions5.txt")
    
#}

# if(!exists("actions")){
#     if(!file.exists(file.path(system.file(package = 'knee'), "data", "9606.protein.actions.v11.0.txt"))){
#         actions <- knee::collect('https://stringdb-static.org/download/protein.actions.v11.0/9606.protein.actions.v11.0.txt.gz')
#     } else {
#         
#         actions <- data.table::fread(file.path(system.file(package = 'knee'), "data", "9606.protein.actions.v11.0.txt"),
#                                      header = T)
#     }
#     assign("actions", actions, envir = .GlobalEnv)
# }

dframe <- function(dt, r){
    dt <- dt[!is.na(get(r))]                    # Remove NA
    rn <- dt[, as.character(.SD[[r]])]          # Rownames to vector
    df <- as.data.frame(dt[,-..convertColumns]) # 
    rownames(df) <- rn
    
    return(df)
}

# Build sample info ----

generate_experiment <- function(){
    
    experiment <- new("Experiment",
                      expData = data_wide,
                      expDataSrc = data_origin,
                      annotationData = data_annotation,
                      aliases = data.table(),
                      sampleData = samples,
                      contrast = contrast,
                      pathwayData = res
    )
    
    assign("experiment", experiment, envir = .GlobalEnv)
    
}

group <- list()
sample_data <- function(data) {
    samples <- names(data_wide[, -..convertColumns])
    condition <- as.factor(gsub('_.*', '', samples))
    replicate <- as.factor(gsub('.*_', '', samples))
    samples <- data.frame(samples, condition, replicate)
    rownames(samples) <- samples$samples

    group <- sapply(levels(condition), function(x) paste("condition", x, sep = ''))
    
    assign("samples", samples, envir = .GlobalEnv)
    assign("condition", condition, envir = .GlobalEnv)
    assign("replicate", replicate, envir = .GlobalEnv)
    assign("group", group, envir = .GlobalEnv)
    
    return (FALSE)
}



filter_na <- function(threshold) {
    
    # Which elements are NA?
    allNA <- is.na(data_origin[, -..convertColumns])
    
    # Summary of how many TRUEs there are in each row
    NA_frequency <- table(rowSums(allNA))
    
    # Subset to NA threshold ----
    
    subset_NA <- function(condition)
    {
        # Subset columns by condition
        condition_subset <- data_origin[,grep(condition,names(data_origin)), with = F]
        
        # Determine if rows pass NA threshold
        rows_to_keep <- rowSums(is.na(condition_subset)) <= threshold # !! set global
        
        # Subset rownames
        keep <- data_origin[rows_to_keep, UNIPROTID]
        
        return(keep)
        
    }
    
    # Apply function to all regions
    condition_sub <- lapply(names(group), subset_NA)
    
    # Reduce to shared proteins
    condition_sub <- Reduce(intersect, condition_sub)
    
    # Subset dataframe
    data_wide <- data_origin[UNIPROTID %in% condition_sub,]
    
    return(data_wide)
    
}



# PCA ----
plotPCA <- function(contribs, ellipse, type) {
    
    dt <- dframe(data_origin, sID)
    
    pca.data <- log2(dt)  # Log2 transform data
    pca.data[is.na(pca.data)] <- 0 # "Impute" missing values as 0
    pca.data <- t(pca.data)        # # Transpose dataset
    p.pca <- prcomp(pca.data, center = TRUE, scale. = TRUE) # Perform principal component analysis
    
    # Biplot extension displaying top contributing proteins currently only available for 2D plot.
    
    if(type == '2d') {
        
        pcaplot <- factoextra::fviz_pca_biplot(p.pca, title = '', label = "var", habillage = condition,
                                               addEllipses = TRUE, ellipse.level = ellipse,
                                               select.var = list(contrib = contribs), repel = TRUE)
        
        return (pcaplot)
        
    } else if (type == '3d') {
        pcaplot <- plotly::plot_ly(x = p.pca$x[,1],
                                   y = p.pca$x[,2],
                                   z = p.pca$x[,3],
                                   color = condition,
                                   colors = c("red","green","blue"),
                                   sizes = c(100, 150)) %>%
            plotly::add_markers() %>%
            plotly::layout(scene = list(xaxis = list(title = 'PC1'),
                                        yaxis = list(title = 'PC2'),
                                        zaxis = list(title = 'PC3')))
        
        return (pcaplot)
    } else if (type == 'UMAP') {
        
        um <- umap::umap(pca.data, n_neighbors = ncol(dt))
        
        df <- data.frame(x = um$layout[,1],
                         y = um$layout[,2],
                         Sample <- condition)
        
        ggplot(df, aes(x, y, colour = condition, shape = condition)) +
            geom_point(size = 4)
        

        
    } else { return (FALSE) }
    
}




# Differential expression ----
diff_exp <- function(coeff, pairing) {
    
    best_fit = 'normal'
    
    if(best_fit == 'nbinom') {
        dds <- DESeqDataSetFromMatrix(countData  = data_wide,
                                      colData    = samples,
                                      design     = ~ condition + replicate)
        dds <- DESeq(dds)
    }
    
    # Create Annotation data and expression set (Biobase)
    phenoData <- new("AnnotatedDataFrame", data = samples)
    exampleSet <- ExpressionSet(assayData = as.matrix(log2(dframe(data_wide, sID))), phenoData = phenoData)
    
    unpaired <- model.matrix( ~ 0 + condition )
    paired <- model.matrix( ~ 0 + condition + replicate )
    
    if(pairing == 1) {
        design <- paired
    } else {
        design <- unpaired
    }
    
    # Fit the linear model
    fit <- lmFit(exampleSet, design)
    
    # Decide possible contrasts
    c <- expand.grid(group, group)
    cc <- factor(ifelse(c$Var1 != c$Var2, paste(c$Var1, c$Var2, sep = '-'), NA ))
    cc <- cc[!is.na(cc)]
    names(cc) <- gsub('-','', gsub('condition','',cc))
    
    cont.matrix <- makeContrasts(contrasts = cc, levels = design) # All possible contrasts
    
    # Contrast groups, run empirical bayes statistics
    fit.cont <- contrasts.fit(fit,cont.matrix)
    fit.cont <- eBayes(fit.cont, robust = T)
    
    # Generate data frame with results from linear model fit, with confidence intervals.
    contrast <- toptable(fit.cont, number = Inf, coef = coeff, confint = TRUE)
    
    # Confidence intervals used for plot, global var
    cint <- contrast
    cint$protein <- rownames(cint)
    cint$protein <- factor(cint$protein, levels = cint$protein[order(cint$logFC)])
    
    assign("cint", cint, envir = .GlobalEnv)
    
    contrast <- contrast[order(contrast$P.Value, decreasing = F),]
    contrast <- data.table::as.data.table(contrast, keep.rownames = T)
    
    setkeyv(contrast, "rn")
    setkeyv(data_wide, sID)
    contrast <- contrast[data_wide[,..convertColumns], nomatch = 0]
    names(contrast) <- c(sID, names(contrast[, 2:ncol(contrast)]))
    setcolorder(contrast, c(convertColumns, "logFC", "CI.L", "CI.R", "t", "P.Value", "adj.P.Val", "B"))

    return( contrast )
    
}










