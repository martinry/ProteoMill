library(igraph)
source("bin/pathway.enrichment.R")

# Background data ----

tissues <- data.table::fread('~/Large_files/tissues.uniprot.csv', sep = '\t', header = T, data.table = F, strip.white=TRUE)
tissue_names <- colnames(tissues)

# Reactome
reactome <- data.table::fread('~/Large_files/reactome/reactome.uniprot.levels.csv', sep = '\t', header = F, data.table = F, strip.white=TRUE)
reactome <- reactome[,2:ncol(reactome)]
colnames(reactome) <- c("UniprotID", "ReactomeID", "URL", "Pathway_name", "Evidence_code", "Species", "Top", "TopPathway")

# Lowest abstraction level
lowest <- data.table::fread('~/Large_files/reactome/reactome.uniprot.lowest.levels.csv', sep = '\t', header = F, data.table = F, strip.white = T)
lowest <- lowest[,2:ncol(lowest)]
colnames(lowest) <- c("UniprotID", "ReactomeID", "URL", "Pathway_name", "Evidence_code", "Species", "Top", "TopPathway")


# Read input file ----

read_file <- function(infile, separator, type) {
    #assign('infile', infile, envir = .GlobalEnv) # Make infile global var
    #assign('separator', separator, envir = .GlobalEnv) # Make infile global var
    
    if(type == "main") {
        data_wide <- data.table::fread(infile, sep = separator, dec = '.', header = T, data.table=FALSE)
        
        rownames(data_wide) <- data_wide[,1]
        data_wide <- data_wide[,2:ncol(data_wide)]
        data_wide[data_wide=="Filtered"] <- NA
        
        empty_rows <- apply(data_wide, 1, function(x) all(is.na(x)))
        
        data_wide <- data_wide[!empty_rows,]
        
        # Set factor -> numeric
        data_wide[] <- lapply(data_wide, function(x) {
            if(is.factor(x)) as.numeric(as.character(x)) else x
        })
        
        assign('data_wide', data_wide, envir = .GlobalEnv)
        assign('data_origin', data_wide, envir = .GlobalEnv)
    } else if(type == "anno") {
        data_annotation <- read.csv(infile, sep = separator, dec = '.')
        assign('data_annotation', data_annotation, envir = .GlobalEnv)
    }
    
    
    
}

# Build sample info ----

groups <- list()
sample_data <- function(data) {
    samples <- data.frame(colnames(data_wide))
    colnames(samples) = "SampleNames"
    rownames(samples) <- samples$SampleNames
    samples$condition <- as.factor(gsub('_.*', '', samples$SampleNames))
    samples$replicate <- as.factor(gsub('.*_', '', samples$SampleNames))
    
    condition <- factor( samples$condition )
    replicate <- factor( samples$replicate )
    
    groups <- sapply(levels(condition), function(x) paste("condition", x, sep = ''))
    
    assign("samples", samples, envir = .GlobalEnv)
    assign("condition", condition, envir = .GlobalEnv)
    assign("replicate", replicate, envir = .GlobalEnv)
    assign("groups", groups, envir = .GlobalEnv)
    
    return (FALSE)
}

# Filter NA ----
filter_na <- function(threshold) {
    
    # Which elements are NA?
    allNA <- is.na(data_wide)
    
    # Summary of how many TRUEs there are in each row
    NA_frequency <- table(rowSums(allNA))
    
    # Subset to NA threshold ----
    
    subset_NA <- function(condition)
    {
        # Subset columns by condition
        condition_subset <- data_wide[,grep(condition,colnames(data_wide))]
        
        # Determine if rows pass NA threshold
        rows_to_keep <- rowSums(is.na(condition_subset)) <= threshold # !! set global
        
        # Subset rownames
        keep <- rownames(condition_subset)[rows_to_keep]
        
        return(keep)
        
    }
    
    # Apply function to all regions
    condition_sub <- lapply(names(groups), subset_NA)
    
    # Reduce to shared proteins
    condition_sub <- Reduce(intersect, condition_sub)
    
    # Subset dataframe
    data_wide <- data_wide[condition_sub,]
    
    # Set factor -> numeric
    data_wide[] <- lapply(data_wide, function(x) {
        if(is.character(x) || is.factor(x)) as.numeric(as.character(x)) else x
    })
    
    data_wide <- log2(data_wide)
    
    return(data_wide)
    
}


# PCA ----
plotPCA <- function(contribs, ellipse, type) {
    
    pca.data <- log2(data_origin)  # Log2 transform data
    pca.data[is.na(pca.data)] <- 0 # "Impute" missing values as 0
    pca.data <- t(pca.data)        # # Transpose dataset
    p.pca <- prcomp(pca.data, center = TRUE, scale. = TRUE) # Perform principal component analysis
    
    # Biplot extension displaying top contributing proteins currently only available for 2D plot.
    
    if(type == '2d') {
        
        pcaplot <- factoextra::fviz_pca_biplot(p.pca, title = '', label = "ind", habillage = condition,
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
    c <- expand.grid(groups, groups)
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
    
    #contrast$FDR <- p.adjust(contrast$P.Value, n = 510)
    
    return( contrast )
    
}










