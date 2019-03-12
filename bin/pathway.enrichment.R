
run_pathway_enrichment <- function(db, background_df) {

  contrast.sign <- contrast[contrast$adj.P.Val < 0.01,]
  contrast.sign <- rownames(contrast.sign)

  specify_decimal <- function(x, k) as.numeric( trimws(format(round(x, k), nsmall=k, scientific = F)) )

  unique_pathways <- unique(reactome$ReactomeID)

  M_annotated_background <- vector()
  N_num_background <- length(background_df)
  n_num_de <- length(contrast.sign)
  k_num_annotated <- vector()

  enriched_pathways <- data.frame(matrix(0, nrow = length(unique_pathways), ncol = 7))
  colnames(enriched_pathways) <- c("Pathway_name", "M", "n", "Pvalue", "iScore", "Pathway_topname", "Genes")

  for(i in 1:length(unique_pathways)){
    pathway1 <- reactome[reactome$ReactomeID == unique_pathways[i],]
    M_annotated_background <- length( pathway1[pathway1$UniprotID %in% background_df, 'ReactomeID'] )
    
    k_annotated_de <- length( pathway1[pathway1$UniprotID %in% contrast.sign, 'ReactomeID'] )
    
    
    pval <- phyper(k_annotated_de, n_num_de, N_num_background, M_annotated_background, lower.tail=FALSE) +
      dhyper(k_annotated_de, n_num_de, N_num_background, M_annotated_background)

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

  enriched_pathways <- enriched_pathways[order(enriched_pathways$iScore, decreasing = T),]

  interesting_pathways <- enriched_pathways[order(enriched_pathways$iScore, decreasing = T),][1:150,]

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
      p.i <- interesting_pathways[[i,1]]
      p.j <- interesting_pathways[[j,1]]
      n.i <- length(p.i)

      x <- intersect(p.i, p.j)

      fraction <- length(x)/n.i

      mat[i,j] <- fraction
    }
  }

  pathway_similarity_hm <- pheatmap::pheatmap(mat, cluster_rows = F, cluster_cols = F, fontsize = 7.5)

  return( pathway_similarity_hm )

}

run_volcano_plot <- function(interesting_pathways) {

  # Volcano plot
  contrast$Rep.Path <- '* NOT SIGNFICANT'
  contrast$Rep.Path.Score <- 0

   for (i in 1:nrow(interesting_pathways)) {
     path_name <- interesting_pathways[i,'Pathway_name']
     path_topname <- interesting_pathways[i,'Pathway_topname']
     score <- interesting_pathways[i,'iScore']

     prots <- unlist(interesting_pathways[i,'Genes'])

     for (p in 1:length(prots)) {
       protein <- prots[p]

       if (score > contrast[rownames(contrast) == protein,'Rep.Path.Score']) {
         contrast[rownames(contrast) == protein,'Rep.Path.Name'] <- path_name
         contrast[rownames(contrast) == protein,'Rep.Path.Top'] <- path_topname
         contrast[rownames(contrast) == protein,'Rep.Path.Score'] <- score
       }
     }
   }

   contrast$Gene <- rownames(contrast)
   results = dplyr::mutate(contrast, sig=ifelse(contrast$adj.P.Val<0.01, "FDR<0.01", "Not Sig"))

   assign("results", results, envir = .GlobalEnv)

   #coul = brewer.pal(11, "Spectral")

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


