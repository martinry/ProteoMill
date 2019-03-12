library(igraph)
source("bin/pathway.enrichment.R")

# Background data ----

tissues <- data.table::fread('~/Large_files/tissues.uniprot.csv', sep = '\t', header = T, data.table = F, strip.white=TRUE)
tissue_names <- colnames(tissues)

# Reactome
reactome <- data.table::fread('~/Large_files/reactome/reactome.uniprot.levels.csv', sep = '\t', header = F, data.table = F, strip.white=TRUE)
reactome <- reactome[,2:ncol(reactome)]
colnames(reactome) <- c("UniprotID", "ReactomeID", "URL", "Pathway_name", "Evidence_code", "Species", "Top", "TopPathway")
# reactome <- data.table::fread('~/Large_files/reactome/uniprot.reactome.hsa.TAS.txt', sep = '\t', header = F, data.table = F, strip.white=TRUE)
# colnames(reactome) <- c("UniprotID", "ReactomeID", "URL", "Pathway_name", "Evidence_code", "Species")
# 
# # Reactome hierarchy
# hier <- data.table::fread('~/Large_files/reactome/reactome.relations.hsa.txt', sep = '\t', header = F, data.table = F, strip.white=TRUE)
# 
# G <- graph.data.frame(d = data.frame(hier$V2, hier$V1), directed = T)
# #plot(df.g, vertex.label = V(df.g)$name)
# 
# reactome$Top <- ''
# 
# top_level <- function(id) {
#     SUB = induced_subgraph(G, subcomponent(G, id, mode="out"))
#     top = farthest.nodes(SUB)$vertices[2]$name
#     return(top)
# }
# 
# seen <- c()
# for(i in 1:nrow(reactome)) {
#     id <- reactome[i,"ReactomeID"]
#     if(!id %in% seen) {
#         print(paste(i, "out of", nrow(reactome), sep = ' '))
#         top <- top_level(reactome[i,"ReactomeID"])
#         reactome[reactome$ReactomeID == id, "Top"] <- top
#         seen <- c(seen, id)
#     }
#     
# }
# 
# toppaths <- qob::switch.items(reactome$Top, reactome, 2, 4, no_targets = "NA", multiple_targets = "first")
# reactome$TopPathways <- toppaths
# 
# write.table(file = "reactome.uniprot.levels.csv", reactome, sep = '\t', quote = F)

# repl <- function() {
#     tissues2 <- data.frame()
#     for(i in 1:ncol(tissues)) {
#         print(paste("Running ", i, " out of ", ncol(tissues), sep = ''))
#         tmpvec <- qob::mapify(tissues[,i], source_id = "Gene_Name", target_id = "UniprotAC")
#         length(tmpvec) <- 17000
#         tissues2[1:length(tmpvec),i] <- tmpvec
#         
#     }
#     
#     return (tissues2)
# }
# 
# temp <- repl()
# 
# colnames(temp) <- colnames(tissues)
# 
# write.table(file = "tissues.uniprot2.csv", temp, sep = '\t', quote = F, row.names = F)

# Read input file ----

read_file <- function(infile, separator) {
    assign('infile', infile, envir = .GlobalEnv) # Make infile global var
    assign('separator', separator, envir = .GlobalEnv) # Make infile global var
    
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

}

# Build sample info ----
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


# Differential expression ----
diff_exp <- function(coeff, pairing) {
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
    
    c <- expand.grid(groups, groups)
    cc <- factor(ifelse(c$Var1 != c$Var2, paste(c$Var1, c$Var2, sep = '-'), NA ))
    cc <- cc[!is.na(cc)]
    names(cc) <- gsub('-','', gsub('condition','',cc))
    
    
    cont.matrix <- makeContrasts(contrasts = cc, levels = design) # All possible contrasts
    
    
    fit.cont <- contrasts.fit(fit,cont.matrix)
    fit.cont <- eBayes(fit.cont, robust = T)
    
    res.all <- topTable(fit.cont, n = Inf)
    res.all.sign <- res.all[res.all$adj.P.Val < 0.01,]
    
    contrast <- topTable(fit.cont, number = Inf, coef = coeff)
    
    contrast <- contrast[order(contrast$P.Value, decreasing = F),]
    
    contrast$FDR <- p.adjust(contrast$P.Value, n = 510)
    
    return( contrast )
    
}