
collect <- function (source_fp, target_fp = "") {
	target_fp <- file.path(basename(source_fp))
	out_name <- gsub(".gz", "", target_fp)
	if (file.exists(out_name)) {
		return(out_name)
	}
	else {
		download.file(url = source_fp, destfile = target_fp)
		if (grepl(".gz", target_fp, fixed = T)) {
			R.utils::gunzip(target_fp, overwrite = F)
		}
		return(out_name)
	}
}

hierarchy <- function(){

  require(data.table)

  # Pathways hierarchy relationship
  REACTOME_hierarchy <- data.table::fread(collect("https://reactome.org/download/current/ReactomePathwaysRelation.txt"),
                                          header = F)

  # All levels of the pathway hierarchy
  # https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt
  #
  REACTOME_all <- data.table::fread(collect("https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt"))
  colnames(REACTOME_all) <- c("UniprotID",
                              "ReactomeID",
                              "URL",
                              "Pathway_name",
                              "Evidence_code",
                              "Species")
  REACTOME_all <- REACTOME_all[Species == "Homo sapiens"]
  REACTOME_all <- REACTOME_all[Evidence_code == "TAS"]

  # Lowest level pathway diagram / Subset of the pathway


  REACTOME_low <- data.table::fread(collect("https://reactome.org/download/current/UniProt2Reactome.txt"))

  colnames(REACTOME_low) <- c("UniprotID",
                              "ReactomeID",
                              "URL",
                              "Pathway_name",
                              "Evidence_code",
                              "Species")

  REACTOME_low <- REACTOME_low[Species == "Homo sapiens"]
  REACTOME_low <- REACTOME_low[Evidence_code == "TAS"]
  REACTOME_low <- REACTOME_low[ReactomeID %in% REACTOME_hierarchy$V2]
  
  protein_descriptions <- data.table::fread(collect("https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz"))

  require(igraph)
  G <- graph.data.frame(d = data.frame(REACTOME_hierarchy$V2, REACTOME_hierarchy$V1), directed = T)

  top_level <- function(id) {
    SUB = induced_subgraph(G, subcomponent(G, id, mode="out"))
    top = farthest.nodes(SUB)$vertices[2]$name
    return(top)
  }


  REACTOME_uniq_low <- unique(REACTOME_low$ReactomeID)
  REACTOME_uniq_high <- sapply(REACTOME_uniq_low, FUN = top_level)

  REACTOME_low$TopReactomeID <- REACTOME_uniq_high[REACTOME_low$ReactomeID]

  REACTOME_ref <- REACTOME_all[,c("ReactomeID", "Pathway_name")]
  REACTOME_ref <- REACTOME_ref[!duplicated(ReactomeID)]

  setkey(REACTOME_low, TopReactomeID)
  setkey(REACTOME_ref, ReactomeID)

  REACTOME_low$TopReactomeName <- REACTOME_low[REACTOME_ref, nomatch = 0, "i.Pathway_name"]

  outpath <- file.path("lib/")

  fwrite(REACTOME_all, file = file.path(outpath, "REACTOME_all.tsv.gz"), sep = '\t', compress = "gzip")
  fwrite(REACTOME_low, file = file.path(outpath, "REACTOME_low.tsv.gz"), sep = '\t', compress = "gzip")
  fwrite(protein_descriptions, file = file.path(outpath, "protein_descriptions.txt.gz"), sep = '\t', compress = "gzip")

}

get_interactions <- function(){
	
	actions <- collect('https://stringdb-static.org/download/protein.actions.v11.0/9606.protein.actions.v11.0.txt.gz')
	actions <- fread(actions)
	
	interactions <- data.table::fread(collect("https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz"), sep = " ", 
									  header = T)
	
	uniprot_to_string_src <- fread("human.uniprot_2_string.2018.tsv")
	
	uniprot_to_string_src$up <- gsub("\\|.*", "", uniprot_to_string_src$V2)
	
	uniprot_to_string_src <- uniprot_to_string_src[, c(3,6)]
	
	setkey(interactions, protein1)
	setkey(uniprot_to_string_src, V3)
	
	interactions2 <- interactions[uniprot_to_string_src, nomatch = 0]
	
	setkey(interactions2, protein2)
	
	interactions3 <- interactions2[uniprot_to_string_src, nomatch = 0]
	
	keycols = c("protein1","protein2")
	setkeyv(interactions3, keycols)
	
	keycols = c("item_id_a","item_id_b")
	setkeyv(actions, keycols)
	
	interactions4 <- interactions3[actions, nomatch = 0]
	
	fwrite(interactions4)
	
	interactions5 <- interactions4[, c(4:10)]
	
	colnames(interactions5) <- c("protein1", "protein2", "mode", "action", "is_directional", "a_is_acting", "score")
	
	interactions5$score <- interactions5$score / 100
	
	fwrite(interactions5, "interactions5.txt", compress = "gzip")
	fwrite(actions, "9606.protein.actions.v11.0.txt.gz", compress = "gzip")
	
}