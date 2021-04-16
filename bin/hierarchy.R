
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
	
	#actions <- collect('https://stringdb-static.org/download/protein.actions.v11.0/9606.protein.actions.v11.0.txt.gz')
	actions <- fread("lib/9606/9606.protein.actions.v11.0.txt.gz")
	
	# interactions <- data.table::fread(collect("https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz"), sep = " ", 
	# 								  header = T)
	
	interactions <- fread("lib/9606/9606.protein.links.v11.0.txt.gz")
	
	#protein_descriptions <- data.table::fread(collect("https://string-db.org/mapping_files/uniprot/all_organisms.uniprot_2_string.2018.tsv.gz"))
	
	uniprot_to_string_src <- fread("all_organisms.uniprot_2_string.2018.tsv.gz", skip = 1)
	
	colnames(uniprot_to_string_src) <- c("Species", "UNIPROTID", "STRINGID", "IDENTITY", "BITSCORE")
	
	uniprot_to_string_src_human <- uniprot_to_string_src[Species == 9606]
	
	uniprot_to_string_src_human$UNIPROTID <- sub("\\|.*", "", uniprot_to_string_src_human$UNIPROTID)
	
	uniprot_to_string_src_human <- uniprot_to_string_src_human[, c("UNIPROTID", "STRINGID")]
	
	setkey(interactions, protein1)
	setkey(uniprot_to_string_src_human, STRINGID)
	
	interactions2 <- interactions[uniprot_to_string_src_human]
	
	setkey(interactions2, protein2)
	
	interactions3 <- interactions2[uniprot_to_string_src_human]
	
	keycols = c("protein1","protein2")
	setkeyv(interactions3, keycols)
	
	keycols = c("item_id_a","item_id_b")
	setkeyv(actions, keycols)
	
	interactions4 <- interactions3[actions]
	
	interactions5 <- interactions4[, c(4, 5, 10)]
	
	colnames(interactions5) <- c("Interactor1", "Interactor2", "Score")
	
	interactions5$Score <- interactions5$Score / 100
	
	fwrite(interactions5, "9606/9606.string.interactions.txt.gz", compress = "gzip")
	#fwrite(actions, "9606.protein.actions.v11.0.txt.gz", compress = "gzip")
	
	
	protein_info <- fread("9606/9606.protein.info.v11.0.txt.gz")
	
	setkey(protein_info, protein_external_id)
	setkey(uniprot_to_string_src_human, STRINGID)
	
	protein_info <- protein_info[uniprot_to_string_src_human]
	
	protein_info <- protein_info[, c(4,5)]
	
	
	fwrite(protein_info, "lib/9606/9606.protein.info.txt.gz", compress = "gzip")

	
}
