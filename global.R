library(igraph)
library(visNetwork)


# Gene annotation data ----

# tissues <- data.table::fread('~/Large_files/tissues.uniprot.csv', sep = '\t', header = T, data.table = F, strip.white=TRUE)
# tissue_names <- colnames(tissues)


# Interaction data ----

if(!exists("interactions")){
    
    interactions <- data.table::fread("C://Users/martinry/interactions6.txt")
    assign("interactions", interactions, envir = .GlobalEnv)
}

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


# Read input file ----

read_file <- function(infile, separator, type) {
    if(type == "main") {
        data_wide <- data.table::fread(infile, sep = separator, dec = '.', header = T)
        
        #rownames(data_wide) <- data_wide[,1]
        #data_wide <- data_wide[,2:ncol(data_wide)]
        #data_wide[data_wide=="Filtered"] <- NA

        empty_rows <- apply(data_wide, 1, function(x) all(is.na(x)))
        
        data_wide <- data_wide[!empty_rows,]
        
        # Set factor -> numeric
        data_wide[] <- lapply(data_wide, function(x) {
            if(is.factor(x)) as.numeric(as.character(x)) else x
        })
        
        assign('data_wide', data_wide, envir = .GlobalEnv)
        assign('data_origin', data_wide, envir = .GlobalEnv)
        
    } else if(type == "anno") {
        data_annotation <- read.csv(infile, sep = separator, dec = '.', row.names = 1)
        assign('data_annotation', data_annotation, envir = .GlobalEnv)
    }
    
    
    
}

# Build sample info ----

group <- list()
sample_data <- function(data) {
    samples <- names(data_wide[, -..convertColumns])
    condition <- as.factor(gsub('_.*', '', samples))
    replicate <- as.factor(gsub('.*_', '', samples))

    group <- sapply(levels(condition), function(x) paste("condition", x, sep = ''))
    
    assign("samples", samples, envir = .GlobalEnv)
    assign("condition", condition, envir = .GlobalEnv)
    assign("replicate", replicate, envir = .GlobalEnv)
    assign("group", group, envir = .GlobalEnv)
    
    return (FALSE)
}



filter_na <- function(threshold) {
    
    # Which elements are NA?
    allNA <- is.na(data_wide[, -..convertColumns])
    
    # Summary of how many TRUEs there are in each row
    NA_frequency <- table(rowSums(allNA))
    
    # Subset to NA threshold ----
    
    subset_NA <- function(condition)
    {
        # Subset columns by condition
        condition_subset <- data_wide[,grep(condition,names(data_wide)), with = F]
        
        # Determine if rows pass NA threshold
        rows_to_keep <- rowSums(is.na(condition_subset)) <= threshold # !! set global
        
        # Subset rownames
        keep <- data_wide[rows_to_keep, UNIPROTID]
        
        return(keep)
        
    }
    
    # Apply function to all regions
    condition_sub <- lapply(names(group), subset_NA)
    
    # Reduce to shared proteins
    condition_sub <- Reduce(intersect, condition_sub)
    
    # Subset dataframe
    data_wide <- data_wide[UNIPROTID %in% condition_sub,]
    
    # # Set factor -> numeric
    # data_wide[] <- lapply(data_wide, function(x) {
    #     if(is.character(x) || is.factor(x)) as.numeric(as.character(x)) else x
    # })
    
    data_wide <- log2(data_wide[, -..convertColumns])
    
    return(data_wide)
    
}

# # Filter NA ----
# filter_na <- function(threshold) {
#     
#     d <- data_wide
#     
#     for(i in seq_along(names(group))){
#         
#         g <- names(group[i])
#         
#         g <- samples[startsWith(samples, g)]
#         
#         dt2 <- d[, lapply(.SD, function(x) is.na(x)), .SDcols = g]
#         dt2[, `:=`(SUM = rowSums(.SD)), .SDcols = g]
#         
#         d <- d[dt2[, SUM <= threshold]]
#         
#     }
#     
#     return(d)
# 
#     
# }


# PCA ----
plotPCA <- function(contribs, ellipse, type) {
    
    pca.data <- log2(data_origin)  # Log2 transform data
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
        
        um <- umap::umap(pca.data, n_neighbors = ncol(data_origin))
        
        df <- data.frame(x = um$layout[,1],
                         y = um$layout[,2],
                         Sample <- condition)
        
        ggplot(df, aes(x, y, colour = condition, shape = condition)) +
            geom_point(size = 4)
        

        
    } else { return (FALSE) }
    
}




# Differential expression ----
diff_exp <- function(coeff, pairing) {
    
    # Which model should be used? Currently (May 2019) limma is used as default. Additional testing needed here.
    best_fit = 'normal'
    
    if(best_fit == 'nbinom') {
        dds <- DESeqDataSetFromMatrix(countData  = data_wide,
                                      colData    = samples,
                                      design     = ~ condition + replicate)
        dds <- DESeq(dds)
    }
    
    # Create Annotation data and expression set (Biobase)
    phenoData <- new("AnnotatedDataFrame", data=samples)
    exampleSet <- ExpressionSet(assayData=as.matrix(data_wide), phenoData=phenoData)
    
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
    
    #contrast$FDR <- p.adjust(contrast$P.Value, n = 510)
    
    return( contrast )
    
}










