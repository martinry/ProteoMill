
get_top_pathways <- function(pathways){

  # pathways | "ReactomeID", "genes", "p.adj"

  pathways$p.adj <- as.numeric(pathways$p.adj)

  pathways <- pathways[, unlist(genes),
                 by = .(ReactomeID, Pathway_name, TopReactomeName, p.adj)]


  pathways <- pathways[, .SD[which.min(p.adj)], by = V1]

  colnames(pathways) <- c("Gene", "ReactomeID", "Pathway_name", "TopReactomeName", "P.Adj")

  return(pathways)

}


enrichment_results <- function(UPREGULATED_pathways, DOWNREGULATED_pathways, contrast){

  dt_up <- UPREGULATED_pathways[ lengths(genes) > 0L, c("Pathway_name", "TopReactomeName", "ReactomeID", "genes", "p.adj")]
  dt_up <- get_top_pathways(dt_up)

  dt_down <- DOWNREGULATED_pathways[ lengths(genes) > 0L, c("Pathway_name", "TopReactomeName", "ReactomeID", "genes", "p.adj")]
  dt_down <- get_top_pathways(dt_down)

  dt_all <- rbindlist(list(dt_up, dt_down))

  # dt_all <- merge.data.table(dt_all, contrast[,c("UNIPROTID", "logFC", "CI.L", "CI.R", "P.Value")], by.x = "Gene", by.y = "UNIPROTID")
  dt_all <- merge.data.table(dt_all, contrast[,c("UNIPROTID", "logFC", "P.Value")], by.x = "Gene", by.y = "UNIPROTID")


  return(dt_all)

}


volcano <- function(res, abstraction = "Global", sID){

  res[P.Adj >= 0.05, "TopReactomeName"] <- "[No significant over-representation]"
  res[P.Adj >= 0.05, "Pathway_name"] <- "[No significant over-representation]"

  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- c("grey", "black", col_vector)

  if(abstraction == "Global"){
    abstraction <- "TopReactomeName"
  } else {
    abstraction <- "Pathway_name"
  }


  sorted_paths <- sort(unique(res[,get(abstraction)]))
  sorted_path_names <- col_vector[seq_along(sorted_paths)]
  names(sorted_path_names) <- sorted_paths
  res$col <- sorted_path_names[res[,get(abstraction)]]

  key <- res[, as.character(.SD[[..sID]])]

  volcano_plot <- ggplot(res, aes(x = logFC, y = -log10(P.Value), key = key)) +
    geom_point(shape = 21, size = 2.5, stroke = .2, aes_string(fill = abstraction)) +
    #geom_errorbarh(aes(xmin = CI.L, xmax = CI.R, color = TopReactomeName), height = .02) +
    scale_fill_manual(values = sorted_path_names) +
    theme_bw(base_size = 10) +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 7.5)) +
    geom_hline(yintercept = 1.3,
               colour = "#990000",
               linetype = "dashed") +
    geom_vline(xintercept = 1,
               colour = "#990000",
               linetype = "dashed") +
    geom_vline(xintercept = -1,
               colour = "#990000",
               linetype = "dashed") +
    scale_y_continuous(trans = "log1p") +
    guides(fill = guide_legend(ncol = 1))

  return(list("volcano_plot" = volcano_plot, "res" = res))


}





