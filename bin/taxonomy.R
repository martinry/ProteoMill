
#	Readme ----
#	
#	Thu May 25, 2021
#	
#	Currently supported species:
#							 Name Abbreviation   DbID TaxID
#	1:                  Bos taurus          BTA  48898  9913
#	2:      Caenorhabditis elegans          CEL  68320  6239
#	3:            Canis familiaris          CFA  49646  9615
#	4:                 Danio rerio          DRE  68323  7955
#	5:    Dictyostelium discoideum          DDI 170941 44689	IDENTIFIER DATA NOT AVAILABLE
#	6:     Drosophila melanogaster          DME  56210  7227
#	7:               Gallus gallus          GGA  49591  9031
#	8:                Homo sapiens          HSA  48887  9606
#	9:                Mus musculus          MMU  48892 10090
#	10: Mycobacterium tuberculosis          MTU 176806  1773	IDENTIFIER DATA NOT AVAILABLE
#	11:      Plasmodium falciparum          PFA 170928  5833	IDENTIFIER DATA NOT AVAILABLE
#	12:          Rattus norvegicus          RNO  48895 10116
#	13:   Saccharomyces cerevisiae          SCE  68322  4932
#	14:  Schizosaccharomyces pombe          SPO  68324  4896	IDENTIFIER DATA NOT AVAILABLE
#	15:                 Sus scrofa          SSC  49633  9823
#	16:         Xenopus tropicalis          XTR 205621  8364
#	
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

if(!basename(getwd()) %in% c("ProteoMill", "proteomill-dev")) stop('Please set working directory to "ProteoMill"')

library(data.table)
library(R.utils)
library(beepr)
library(dplyr)
library(rvest)
library(xml2)
library(RSQLite)
library(AnnotationHub)

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



#	Step 0 Species data ----

# Get a list of species from REACTOME
url1 <- "https://reactome.org/content/schema/objects/Species?page=1"
url2 <- "https://reactome.org/content/schema/objects/Species?page=2"

wnodes <- function(url) {
	url %>%
		read_html %>%
		html_nodes("table") %>%
		html_table() %>%
		as.data.table()
}

w1 <- wnodes(url1) 
w2 <- wnodes(url2)

# This table contains REACTOME ID and species name
species_table <- rbindlist(list(w1, w2))

# Get species info (taxonomical ID, etc) by based on above table
rspecies_lookup <- function(species) {
	
	wurl <- paste0('https://reactome.org/content/schema/instance/browser/', species_table[Name == species, Identifier])
	
	wnodes <- wurl %>%
		read_html %>%
		html_nodes("tr > td") %>% html_nodes("span") %>%
		xml_contents() %>% as.character() %>% trimws()
	
	abbr <- wnodes[1]
	dbid <- wnodes[2]
	name <- wnodes[3]
	taxid <- wnodes[length(wnodes)]
	
	data.table("Name" = name, "Abbreviation" = abbr, "DbID" = dbid, "TaxID" = taxid)
}

# Species of interest
sp <- c("Drosophila melanogaster",
		"Rattus norvegicus",
		"Danio rerio",
		"Caenorhabditis elegans",
		"Canis familiaris",
		"Mus musculus",
		"Homo sapiens",
		"Sus scrofa",
		"Bos taurus",
		"Gallus gallus",
		"Plasmodium falciparum",
		"Xenopus tropicalis",
		"Schizosaccharomyces pombe",
		"Dictyostelium discoideum",
		"Saccharomyces cerevisiae",
		"Mycobacterium tuberculosis")

# Subset to matching
species_table <- species_table[Name %in% sp]

# Initialize species_ dt
species_ <- data.table()

for(i in seq(species_table$Name)) {
	print(paste0("Fetch details for ", species_table$Name[i]))
	sp <- rspecies_lookup( species_table$Name[i] )
	species_ <- rbindlist(list(species_, sp))
}

# Exception: Canis familiaris, 9615 -> 9612 (according to STRINGdb)
species_[Name == "Canis familiaris", "TaxID"] <- 9612


#	Step 1 Reactome data ----

#	Download (if not already present locally) and process Reactome data by species
hierarchy <- function(species = "Homo sapiens"){
	
	require(data.table)
	
	# Pathways hierarchy relationship
	if(!exists("xREACTOME_hierarchy")){
		xREACTOME_hierarchy <<- data.table::fread(Bundle("https://reactome.org/download/current/ReactomePathwaysRelation.txt", subdir = "lib"),
												  header = F)
	}
	
	REACTOME_hierarchy <- copy(xREACTOME_hierarchy)
	
	# All levels of the pathway hierarchy
	if(!exists("xREACTOME_all")){
		xREACTOME_all <<- data.table::fread(Bundle("https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt", subdir = "lib"), header = F)
	}
	
	REACTOME_all <- copy(xREACTOME_all)
	
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
	
	if(!exists("xREACTOME_low")){
		xREACTOME_low <<- data.table::fread(Bundle("https://reactome.org/download/current/UniProt2Reactome.txt", subdir = "lib"))	
	}
	
	REACTOME_low <- copy(xREACTOME_low)
	
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

#	Iterate over our species dt and collect data with above function
for(i in 1:species_[, .N]) {
	
	print(paste0("Preparing pathway files for ", species_[i, "Name"][[1]], ". Please wait..."))
	
	hres <- hierarchy(species_[i, "Name"][[1]])
	
	taxid <- species_[i, "TaxID"][[1]]
	
	dir.create(file.path("lib", as.character(taxid)), showWarnings = F)
	
	top_name <- paste0(taxid, "_REACTOME_all.tsv.gz")
	low_name <- paste0(taxid, "_REACTOME_low.tsv.gz")
	
	fwrite(hres$REACTOME_all, file = file.path("lib", taxid, top_name), sep = '\t', compress = "gzip")
	fwrite(hres$REACTOME_low, file = file.path("lib", taxid, low_name), sep = '\t', compress = "gzip")
	
}



#	Step 2 STRINGdb data ----

get_interactions <- function(taxid){
	
	subd <- file.path("lib", taxid)
	
	if(!exists("xinteractions")){
		xinteractions <<- fread(Bundle(paste0("https://stringdb-static.org/download/protein.actions.v11.0/", taxid, ".protein.actions.v11.0.txt.gz"), subdir = subd))
	}
	
	if(!exists("xall_organisms")){
		xall_organisms <<- fread(Bundle("https://string-db.org/mapping_files/uniprot/all_organisms.uniprot_2_string.2018.tsv.gz", subdir = "lib"), skip = 1, header = F)
	}
	
	interactions <- copy(xinteractions)
	all_organisms <- copy(xall_organisms)
	
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
	
	if(!exists("xprotein_info")){
		xprotein_info <<- fread(Bundle(paste0("https://stringdb-static.org/download/protein.info.v11.0/", taxid, ".protein.info.v11.0.txt.gz"), subdir = subd))
	}
	
	if(!exists("xall_organisms")){
		xall_organisms <<- fread(Bundle("https://string-db.org/mapping_files/uniprot/all_organisms.uniprot_2_string.2018.tsv.gz", subdir = "lib"), skip = 1, header = F)
	}
	
	protein_info <- copy(xprotein_info)
	all_organisms <- copy(xall_organisms)
	
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

# Interaction file. TBD merge these functions/loops.
for(i in 1:species_[, .N]) {
	
	print(paste0("Preparing interaction files for ", species_[i, "Name"][[1]], ". Please wait..."))
	
	taxid <- species_[i, "TaxID"][[1]]
	
	interactions <- get_interactions(taxid = taxid)
	
	dir.create(file.path("lib", as.character(taxid)), showWarnings = F)
	
	fname <- paste0(taxid, ".string.interactions.txt.gz")
	
	fwrite(interactions, file = file.path("lib", taxid, fname), compress = "gzip")
	
}

# Description file. TBD merge these functions/loops.
for(i in 1:species_[, .N]) {
	
	print(paste0("Preparing description files for ", species_[i, "Name"][[1]], ". Please wait..."))
	
	taxid <- species_[i, "TaxID"][[1]]
	
	protein_info <- get_protein_descriptions(taxid)
	
	dir.create(file.path("lib", as.character(taxid)), showWarnings = F)
	
	fname <- paste0(taxid, ".protein.info.txt.gz")
	
	fwrite(protein_info, file = file.path("lib", taxid, fname), compress = "gzip")
}



# Get identifier info
ah <- AnnotationHub()

for(i in 1:species_[, .N]) {
	name <- species_[i, "Name"][[1]]
	taxid <- species_[i, "TaxID"][[1]]
	
	print(paste0("Preparing identifier files for ", name, ". Please wait..."))

	edb <- query(ah, c("EnsDb", name))
	
	if(length(edb) > 0){
		ah_id <- edb$ah_id[length(edb$ah_id)]
		ver <- gsub("[A-z ]", "", edb$title[length(edb$ah_id)])
		edb <- ah[[ah_id]]
		
		fname <- paste0(taxid, ".EnsDb.v", ver)
		
		if(!file.exists(file.path("lib", taxid, fname))){
			file.copy(from = file.path(edb@ensdb@dbname), to = file.path("lib", taxid, fname))
		}
	} else {
		print(paste0("No records found for ", name, "..."))
	}
	

	
}



beepr::beep()
print("All done!")
