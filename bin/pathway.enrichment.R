
run_pathway_enrichment <- function(db, background_df) {

  contrast.sign <- contrast[contrast$adj.P.Val < 0.05,]
  contrast.sign <- rownames(contrast.sign)

  specify_decimal <- function(x, k) as.numeric( trimws(format(round(x, k), nsmall=k, scientific = F)) )

  unique_pathways <- unique(reactome$ReactomeID)

  M_annotated_background <- vector()
  N_num_background <- length(background_df)
  n_num_de <- length(contrast.sign)
  k_num_annotated <- vector()

  enriched_pathways <- data.frame(matrix(0, nrow = length(unique_pathways), ncol = 8))
  colnames(enriched_pathways) <- c("Pathway_name", "M", "n", "Pvalue", "AdjP", "iScore", "Pathway_topname", "Genes")

  for(i in 1:length(unique_pathways)){
    pathway1 <- reactome[reactome$ReactomeID == unique_pathways[i],]
    M_annotated_background <- length( pathway1[pathway1$UniprotID %in% background_df, 'ReactomeID'] )
    
    k_annotated_de <- length( pathway1[pathway1$UniprotID %in% contrast.sign, 'ReactomeID'] )
    
    # dhyper(x, m, n, k) gives the probability of drawing exactly x
    # phyper(x, m, n, k) gives the probability of getting x or fewer,
    
    # n in sample
    # n in pop
    # pop size
    # sample size
    
    #pval <- phyper(k_annotated_de, n_num_de, N_num_background, M_annotated_background, lower.tail=FALSE) +
    #  dhyper(k_annotated_de, n_num_de, N_num_background, M_annotated_background)
    
    pval <- phyper(k_annotated_de, M_annotated_background, N_num_background-M_annotated_background, n_num_de, lower.tail = FALSE) +
        dhyper(k_annotated_de, M_annotated_background, N_num_background-M_annotated_background, n_num_de)

    sampled <- pathway1[pathway1$UniprotID %in% contrast.sign, 'UniprotID']

    lfc <- contrast[rownames(contrast) %in% sampled, 'logFC']
    
    lfc.mean <- mean( abs(lfc) )
    
    iscore <- lfc.mean * -log10( pval )
    
    enriched_pathways[i, "Pathway_name"] <- reactome[reactome$ReactomeID == unique_pathways[i], 4][1]
    enriched_pathways[i, "Pathway_topname"] <- reactome[reactome$ReactomeID == unique_pathways[i], 8][1]

    #enriched_pathways[i, "Pathway"] <- ifelse(nchar(unique_pathways[i]) > 50, paste(substr(unique_pathways[i], 1,50), '...', sep = ''), unique_pathways[i])
    enriched_pathways[i, "M"] <- M_annotated_background
    enriched_pathways[i, "n"] <- k_annotated_de
    enriched_pathways[i, "Pvalue"] <- ifelse(is.na(pval), NA, specify_decimal(pval, 9))
    enriched_pathways[i, "iScore"] <- ifelse(is.na(iscore), NA, specify_decimal(iscore, 9))
    enriched_pathways$Genes[i] <- list(sampled)
  }
  
  enriched_pathways$AdjP <- p.adjust(enriched_pathways$Pvalue, method = "BH")
  enriched_pathways <- enriched_pathways[order(enriched_pathways$iScore, decreasing = T),]

  interesting_pathways <- enriched_pathways[order(enriched_pathways$iScore, decreasing = T),][1:150,]
  interesting_pathways <- interesting_pathways[interesting_pathways$Pvalue < 0.01,]

  return( list(enriched_pathways, interesting_pathways) )

}

run_similarity_plot <- function(interesting_pathways) {
    
    interesting_pathways <- interesting_pathways[1:40,]

    # Plot similarity matrix
    mat <- matrix(nrow = nrow(interesting_pathways), ncol = nrow(interesting_pathways))
    rownames(mat) <- interesting_pathways$Pathway_name
    colnames(mat) <- interesting_pathways$Pathway_name

  for(i in 1:nrow(interesting_pathways)) {
    for(j in 1:nrow(interesting_pathways)) {
      p.i <- interesting_pathways[["Genes",8]]
      p.j <- interesting_pathways[["Genes",8]]
      n.i <- length(p.i)

      x <- intersect(p.i, p.j)

      fraction <- length(x)/n.i

      mat[i,j] <- fraction
    }
  }

  return( mat )

}

run_volcano_plot <- function(contrast, results) {

  # Volcano plot


   qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
   col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

   col_vector[1] = "grey"
   col_vector[2] = "black"

   vol <- ggplot(results, aes(x = logFC, y = -log10(P.Value))) +
     geom_point(shape = 21, size = 2.5, aes(fill=Rep.Path.Top)) +
     scale_fill_manual(values=col_vector) +
     theme_bw(base_size = 14) +
     theme(legend.position = "right") +
     geom_hline(yintercept = 1.3,
                colour = "#990000",
                linetype = "dashed") +
     geom_vline(xintercept = 1,
                colour = "#990000",
                linetype = "dashed") +
     geom_vline(xintercept = -1,
                colour = "#990000",
                linetype = "dashed") +
     scale_y_continuous(trans = "log1p")

   volcano_plot <- vol + geom_text_repel(
     data = filter(results, logFC > 3.5),
     aes(label = Gene),
     box.padding = unit(0.5, 'lines'),
     point.padding = unit(1, 'lines'),
     segment.color = '#cccccc',
     segment.size = 0.5,
     arrow = arrow(length = unit(0.01, 'npc'))
   )

  return ( volcano_plot )

}

run_sankey_diagram <- function(results) {
    require(networkD3)
    
    results$regulation <- ifelse(results$logFC > 0, "Up-regulated", "Down-regulated")
    
    results.sign <- results[results$Rep.Path.Score > 0,]
    
    paths <-
        data.frame(
            "Gene" = results.sign$Gene,
            "Regulation" = results.sign$regulation,
            "TopPath" = results.sign$Rep.Path.Top,
            "PathName" = results.sign$Rep.Path.Name
        )
    
    paths2 <- paths[,c(1,3,4)]
    fr1 <- as.data.frame(table(paths[,2:3]))
    fr2 <- as.data.frame(table(paths2[,2:3]))
    
    nodes <- data.frame("name" = c("Up-regulated", "Down-regulated", levels(paths$TopPath), levels(paths$PathName)))
    links <- data.frame("source" = c(), "target" = c(), "value" = c())
    
    ntop <- 3+length(levels(paths$TopPath))
    nname <- ntop + length(levels(paths$PathName)) - 1
    
    for(i in 3:ntop) {
        n <- as.character(nodes[i,])
        k <- i-1
        
        if(n %in% fr1$TopPath) {
            up <- fr1[fr1$Regulation == "Up-regulated",]
            if(up[up$TopPath == n,3] > 0){
                links <- rbind(links, c(0, k, up[up$TopPath == n,3]))
            }
            down <- fr1[fr1$Regulation == "Down-regulated",]
            if(down[down$TopPath == n,3] > 0){
                links <- rbind(links, c(1, k, down[down$TopPath == n,3]))
            }
        }
        
        if(n %in% fr2$TopPath) {
            for(j in (ntop):nname){
                m <- as.character(nodes[j,])
                l <- j-1
                
                tp <- fr2[fr2$TopPath == n & fr2$PathName == m & fr2$Freq > 0,]
                
                if(nrow(tp) > 0) {
                    links <- rbind(links, c(k, l, tp[,3]))
                }
            }
        }
    }
    
    colnames(links) <- c("source", "target", "value")
    
    P <- list("nodes" = nodes, "links" = links)
    
    return (P)
    
}

