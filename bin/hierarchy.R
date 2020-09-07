

hierarchy <- function(){
# 
#   require(data.table)
# 
#   # Pathways hierarchy relationship
#   REACTOME_hierarchy <- data.table::fread(knee:::collect("https://reactome.org/download/current/ReactomePathwaysRelation.txt"),
#                                           header = F)
# 
#   # All levels of the pathway hierarchy
#   # https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt
#   #
#   REACTOME_all <- data.table::fread(knee:::collect("https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt"))
#   colnames(REACTOME_all) <- c("UniprotID",
#                               "ReactomeID",
#                               "URL",
#                               "Pathway_name",
#                               "Evidence_code",
#                               "Species")
#   REACTOME_all <- REACTOME_all[Species == "Homo sapiens"]
#   REACTOME_all <- REACTOME_all[Evidence_code == "TAS"]
# 
#   # Lowest level pathway diagram / Subset of the pathway
# 
# 
#   REACTOME_low <- data.table::fread(knee:::collect("https://reactome.org/download/current/UniProt2Reactome.txt"))
# 
#   colnames(REACTOME_low) <- c("UniprotID",
#                               "ReactomeID",
#                               "URL",
#                               "Pathway_name",
#                               "Evidence_code",
#                               "Species")
# 
#   REACTOME_low <- REACTOME_low[Species == "Homo sapiens"]
#   REACTOME_low <- REACTOME_low[Evidence_code == "TAS"]
#   REACTOME_low <- REACTOME_low[ReactomeID %in% REACTOME_hierarchy$V2]
# 
# 
#   require(igraph)
#   G <- graph.data.frame(d = data.frame(REACTOME_hierarchy$V2, REACTOME_hierarchy$V1), directed = T)
# 
#   top_level <- function(id) {
#     SUB = induced_subgraph(G, subcomponent(G, id, mode="out"))
#     top = farthest.nodes(SUB)$vertices[2]$name
#     return(top)
#   }
# 
# 
#   REACTOME_uniq_low <- unique(REACTOME_low$ReactomeID)
#   REACTOME_uniq_high <- sapply(REACTOME_uniq_low, FUN = top_level)
# 
#   REACTOME_low$TopReactomeID <- REACTOME_uniq_high[REACTOME_low$ReactomeID]
# 
#   REACTOME_ref <- REACTOME_all[,c("ReactomeID", "Pathway_name")]
#   REACTOME_ref <- REACTOME_ref[!duplicated(ReactomeID)]
# 
#   setkey(REACTOME_low, TopReactomeID)
#   setkey(REACTOME_ref, ReactomeID)
# 
#   REACTOME_low$TopReactomeName <- REACTOME_low[REACTOME_ref, nomatch = 0, "i.Pathway_name"]
# 
#   outpath <- file.path(system.file(package = 'knee'), "data")
# 
#   fwrite(REACTOME_all, file = file.path(outpath, "REACTOME_all.tsv.gz"), sep = '\t', compress = "gzip")
#   fwrite(REACTOME_low, file = file.path(outpath, "REACTOME_low.tsv.gz"), sep = '\t', compress = "gzip")

}

