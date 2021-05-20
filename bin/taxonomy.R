
#	Readme ----
#	
#	Thu Apr 15, 2021
#	
#	Currently supported species:
#	Homo sapiens | HSA | 9606
#	Mus musculus | MMU | 10090
#	Rattus norvegicus | RNO | 10116
#	
#	Adding more is dependent on obtaining .sqlite files via ensembldb package and to have UNIPROTID support in each file.
#	
#	This file will collect and prepare annotation data for all currently supported species.
#	Annotation data is collected from reactome and stringdb.

#	File structure ----
#
#	In lib/ there is a corresponding subfolder for each species (named by its taxonomical id, for example human dir is named 9606).
#	Each folder contains the following files (human as example):
#
#			============ Name ===============			=============== Description ===============
#		-	[taxid].protein.info.txt.gz				|	STRINGDB processed file used directly in ProteoMill for protein descriptions (used when hovering a node in network analysis)
#		-	[taxid].protein.info.v11.0.txt.gz		|	STRINGDB raw file used to create above mentioned file.
#		-	[taxid].protein.actions.v11.0.txt.gz	|	STRINGDB raw interaction file
#		-	[taxid].reactome.interactions.txt.gz	|	Reactome interaction file. Currently unused.
#		-	[taxid].string.interactions.txt.gz		|	STRINGDB processed interaction file. Used in network analysis.
#		-	[taxid]_REACTOME_all.tsv.gz				|	Reactome processed pathway file.
#		-	[taxid]_REACTOME_low.tsv.gz				|	Reactome processed pathway file.
#	
#	Additionally, in the top-level lib/ dir these files contain information for all species
#
#			============ Name ===============			=============== Description ===============
#		-	ReactomePathwaysRelation.txt.gz			|	Reactome raw file that contains the hierarchical information for each sets of pathways. Used to determine top-lowest level pathways.
#		-	UniProt2Reactome_All_Levels.txt.gz		|	All levels of the pathway hierarchy.
#		-	UniProt2Reactome.txt.gz					|	Lowest levels of the pathway hierarchy.
#
#	Only some of these files (the ones described as processed) are necessary for running the ProteoMill app, but the raw files are necessary for updates.
#

#	! NB make sure your work directory is set appropriately !

if(basename(getwd()) != "ProteoMill") stop('Please set working directory to "ProteoMill"')

typewrite <- function(m, pace = .03, wait = 1) {
	m <- unlist(strsplit(m, ""))
	for(l in m){cat(l); Sys.sleep(pace)}
	Sys.sleep(wait)
	cat("\n")
}

# Some additional information about this script ----
typewrite("This script will download all necessary external data (interaction data, pathway data, etc.) to your local machine.")
typewrite("The total disk usage in lib/ will be about 430MB.")
typewrite("=================================================", pace = 0.01)
typewrite("This takes a while. When the script has finished and all files are downloaded, the following sound will play:", wait = 2)
beepr::beep()
typewrite("Ping!", wait = 2)
typewrite("=================================================", pace = 0.01)
typewrite("Starting in 3...")
typewrite("2...")
typewrite("1...")
Sys.sleep(1)
typewrite("=================================================", pace = 0.01)
typewrite("Collecting data...")

#	Now let's begin downloading and processing the files.

# 	Function to download and gzip files ----

#	If the file exists but not as .gz, original file will be removed and gzipped version returned.
#	If the file exists as .gz it will be returned as is.
#	If the file does not exist it will be downloaded, gzipped (if not already) and returned.

library(data.table)
library(R.utils)
library(beepr)

Bundle <- function (source_fp, subdir = "") {
	
	target_fp <- file.path(subdir, basename(source_fp))
	
	if(file.exists(target_fp)) {
		if(endsWith(target_fp, ".gz")) {
			return(target_fp)
		} else {
			gzip(target_fp, skip = T, remove = T)
			Bundle(paste0(target_fp, ".gz"), subdir)
		}
	} else if(file.exists(paste0(target_fp, ".gz"))) {
		Bundle(paste0(target_fp, ".gz"), subdir)
	} else {
		print(paste0("File ", basename(source_fp), " not found in ", subdir))
		print(paste0("Downloading from ", source_fp, "..."))
		download.file(url = source_fp, destfile = target_fp, method = "auto", mode = "wb")
		Bundle(source_fp, subdir)
	}
}




#	Step 1 Reactome data ----

#	Download (if not already present locally) and process Reactome data by species
hierarchy <- function(species = "Homo sapiens"){
	
	require(data.table)
	
	# Pathways hierarchy relationship
	REACTOME_hierarchy <- data.table::fread(Bundle("https://reactome.org/download/current/ReactomePathwaysRelation.txt", subdir = "lib"),
											header = F)
	
	# All levels of the pathway hierarchy
	REACTOME_all <- data.table::fread(Bundle("https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt", subdir = "lib"))
	
	colnames(REACTOME_all) <- c("UniprotID",
								"ReactomeID",
								"URL",
								"Pathway_name",
								"Evidence_code",
								"Species")
	
	REACTOME_all <- REACTOME_all[Species == species]
	REACTOME_all_TAS <- REACTOME_all[Evidence_code == "TAS"]
	REACTOME_all_IEA <- REACTOME_all[Evidence_code == "IEA"]
	
	if(REACTOME_all_TAS[, .N] >= REACTOME_all_IEA[, .N]) {
		REACTOME_all <- REACTOME_all_TAS
	} else {
		REACTOME_all <- REACTOME_all_IEA
	}
	
	# Lowest level pathway diagram / Subset of the pathway
	
	
	REACTOME_low <- data.table::fread(Bundle("https://reactome.org/download/current/UniProt2Reactome.txt", subdir = "lib"))
	
	colnames(REACTOME_low) <- c("UniprotID",
								"ReactomeID",
								"URL",
								"Pathway_name",
								"Evidence_code",
								"Species")
	
	REACTOME_low <- REACTOME_low[Species == species]
	
	
	REACTOME_low_TAS <- REACTOME_low[Evidence_code == "TAS"]
	REACTOME_low_IEA <- REACTOME_low[Evidence_code == "IEA"]
	
	if(REACTOME_low_TAS[, .N] >= REACTOME_low_IEA[, .N]) {
		REACTOME_low <- REACTOME_low_TAS
	} else {
		REACTOME_low <- REACTOME_low_IEA
	}
	
	
	REACTOME_low <- REACTOME_low[ReactomeID %in% REACTOME_hierarchy$V2]
	
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
	
	return(list("REACTOME_all" = REACTOME_all, "REACTOME_low" = REACTOME_low))
	
}

#	Make a dt of species of interest
species_ <- data.table(sp = c("Homo sapiens", "Mus musculus", "Rattus norvegicus"),
					   txid = c(9606, 10090, 10116))

#	Iterate over our species dt and collect data with above function
for(i in 1:species_[, .N]) {
	
	print(paste0("Preparing pathway files for ", species_[i][[1]], ". Please wait..."))
	
	hres <- hierarchy(species_[i][[1]])
	
	taxid <- species_[i][[2]]
	
	dir.create(file.path("lib", as.character(taxid)), showWarnings = F)
	
	top_name <- paste0(taxid, "_REACTOME_all.tsv.gz")
	low_name <- paste0(taxid, "_REACTOME_low.tsv.gz")
	
	fwrite(hres$REACTOME_all, file = file.path("lib", taxid, top_name), sep = '\t', compress = "gzip")
	fwrite(hres$REACTOME_low, file = file.path("lib", taxid, low_name), sep = '\t', compress = "gzip")
	
}



#	Step 2 STRINGdb data ----

get_interactions <- function(taxid){
	
	subd <- file.path("lib", taxid)
	
	interactions <- fread(Bundle(paste0("https://stringdb-static.org/download/protein.actions.v11.0/", taxid, ".protein.actions.v11.0.txt.gz"), subdir = subd))
	all_organisms <- fread(Bundle("https://string-db.org/mapping_files/uniprot/all_organisms.uniprot_2_string.2018.tsv.gz", subdir = "lib"), skip = 1, header = F)
	
	colnames(all_organisms) <- c("Species", "UNIPROTID", "STRINGID", "IDENTITY", "BITSCORE")
	
	all_organisms <- all_organisms[Species == taxid]
	
	all_organisms$UNIPROTID <- sub("\\|.*", "", all_organisms$UNIPROTID)
	
	all_organisms <- all_organisms[, c("UNIPROTID", "STRINGID")]
	
	setkey(interactions, item_id_a)
	setkey(all_organisms, STRINGID)
	
	interactions <- interactions[all_organisms]
	
	setkey(interactions,  item_id_b)
	
	interactions <- interactions[all_organisms]
	
	interactions <- interactions[, c(8, 9, 7)]
	interactions <- interactions[!duplicated(interactions)]
	
	#interactions <- interactions[, c(4, 5, 3)]
	
	colnames(interactions) <- c("Interactor1", "Interactor2", "Score")
	
	interactions$Score <- interactions$Score / 100
	
	return(interactions)
	
	
}

get_protein_descriptions <- function(taxid){
	
	subd <- file.path("lib", taxid)
	
	protein_info <- fread(Bundle(paste0("https://stringdb-static.org/download/protein.info.v11.0/", taxid, ".protein.info.v11.0.txt.gz"), subdir = subd))
	
	all_organisms <- fread(Bundle("https://string-db.org/mapping_files/uniprot/all_organisms.uniprot_2_string.2018.tsv.gz", subdir = "lib"), skip = 1, header = F)
	
	colnames(all_organisms) <- c("Species", "UNIPROTID", "STRINGID", "IDENTITY", "BITSCORE")
	
	all_organisms <- all_organisms[Species == taxid]
	
	all_organisms$UNIPROTID <- sub("\\|.*", "", all_organisms$UNIPROTID)
	
	all_organisms <- all_organisms[, c("UNIPROTID", "STRINGID")]
	
	setkey(protein_info, protein_external_id)
	setkey(all_organisms, STRINGID)
	
	protein_info <- protein_info[all_organisms]
	
	protein_info <- protein_info[, c(4,5)]
	
	return(protein_info)
	
}

#	In case we didn't run part one, here is the species list again:

#	Make a dt of species of interest
species_ <- data.table(sp = c("Homo sapiens", "Mus musculus", "Rattus norvegicus"),
					   txid = c(9606, 10090, 10116))


# Interaction file. TBD merge these functions/loops.
for(i in 1:species_[, .N]) {
	
	print(paste0("Preparing interaction files for ", species_[i][[1]], ". Please wait..."))
	
	taxid <- species_[i][[2]]
	
	interactions <- get_interactions(taxid = taxid)
	
	dir.create(file.path("lib", as.character(taxid)), showWarnings = F)
	
	fname <- paste0(taxid, ".string.interactions.txt.gz")
	
	fwrite(interactions, file = file.path("lib", taxid, fname), compress = "gzip")
	
}

# Description file. TBD merge these functions/loops.
for(i in 1:species_[, .N]) {
	
	print(paste0("Preparing description files for ", species_[i][[1]], ". Please wait..."))
	
	taxid <- species_[i][[2]]
	
	protein_info <- get_protein_descriptions(taxid)
	
	dir.create(file.path("lib", as.character(taxid)), showWarnings = F)
	
	fname <- paste0(taxid, ".protein.info.txt.gz")
	
	fwrite(protein_info, file = file.path("lib", taxid, fname), compress = "gzip")
}

beepr::beep()
print("All done!")
