

ora <- function(gene,
                taxid = 9606,
                database = NULL,
                pAdjMethod = "BH") {


  database <- fread(paste0("lib/", taxid, "/", taxid, "_REACTOME_low.tsv.gz"), sep = '\t', header = T)
  
  gene <- gene[gene %in% database$UniprotID]

  # Genes and their counts in input and background, for each pathway
  dt <- database[, list(genes      = list(UniprotID[UniprotID %in% gene]),
                        background = list(UniprotID),
                        q          = length(UniprotID[UniprotID %in% gene]),
                        m          = length(UniprotID)),
                 by = .(ReactomeID, Pathway_name, TopReactomeName)]

  dt <- dt[q > 0]

  # Total number of genes
  dtN = database[!duplicated(UniprotID), .N]

  # Calculate and format p-values
  dt[ , p := mapply(phyper, q-1, m, n = (dtN - m), k = length(gene), lower.tail = F)]
  dt$p.adj  <- p.adjust(dt$p, method = pAdjMethod)
  dt$p      <- formatC(dt$p, format = "e", digits = 2)
  dt$p.adj  <- formatC(dt$p.adj, format = "e", digits = 2)
  dt        <- dt[order(as.numeric(p.adj))]

  res <- new("ora_result",
             gene           = gene,
             database       = database,
             pAdjustMethod  = pAdjMethod,
             output         = dt
  )
  return (res)

}
