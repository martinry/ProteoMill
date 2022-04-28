require(data.table)

source("bin/All_Classes.R")
source("bin/obsolete.R")
source("bin/ora.R")
source("bin/heatmap.R")
source("bin/pca.R")
source("bin/tokenization.R")
#source("bin/differential_expression.R")
#source("auth/auth.R")

# Upload dataset ----

convertColumns <- c("UNIPROTID", "ENTREZID", "SYMBOL", "GENEID", "PROTEINID")
assign("convertColumns", convertColumns, envir = .GlobalEnv)

organisms <- list(
	"Homo sapiens [9606]",
	"Bos taurus [9913]",
	"Caenorhabditis elegans [6239]",
	"Canis familiaris [9612]",
	"Danio rerio [7955]",
	"Drosophila melanogaster [7227]",
	"Gallus gallus [9031]",
	"Mus musculus [10090]",
	"Rattus norvegicus [10116]",
	"Saccharomyces cerevisiae [4932]",
	"Sus scrofa [9823]",
	"Xenopus tropicalis [8364]",
	"Other"
)


# User timeout ----

timeoutMinutes <- 20

inactivity <- sprintf("function idleTimer() {
var t = setTimeout(logout, %s);
window.onmousemove = resetTimer; // catches mouse movements
window.onmousedown = resetTimer; // catches mouse movements
window.onclick = resetTimer;     // catches mouse clicks
window.onscroll = resetTimer;    // catches scrolling
window.onkeypress = resetTimer;  //catches keyboard actions

function logout() {
Shiny.setInputValue('timeOut', '%s minutes of')
}

function resetTimer() {
clearTimeout(t);
t = setTimeout(logout, %s);  // time is in milliseconds (1000 is 1 second)
}
}
idleTimer();", timeoutMinutes*6e4, timeoutMinutes, timeoutMinutes*6e4)


# Global functions ----

separator <- function(s) {
	switch(s,
		   "1" = "auto",
		   "2" = ",",
		   "3" = ";",
		   "4" = "\t")
}

identifier <- function(i) {
	switch(i,
		   "1" = "auto",
		   "2" = "UNIPROTID",
		   "3" = "ENTREZID",
		   "4" = "SYMBOL",
		   "5" = "GENEID",
		   "6" = "PROTEINID")
}

get_delim <- function(c){
	if(grepl("\n", c)){
		return("\n")
	} else if(grepl(",", c)){
		return(",")
	} else if(grepl(";", c)){
		return(";")
	} else {
		return ("None")
	}
}

npals <- function(s) {
	switch(s,
		   "Accent"	= 8,
		   "Dark2" = 8,
		   "Paired" = 12,
		   "Pastel1" = 9,
		   "Pastel2" = 8,
		   "Set1" =	9,
		   "Set2" =	8,
		   "Set3" =	12
	)
}

getFilemd5sum <- function(fp) {
	fp <- file.path(fp)
	return( digest::digest(object = fp, algo = "md5", file = T) )
}

getlibfPaths <- function(taxid, lib, version = "current", token = NULL) {
	if(is.null(token)) {
		if(version == "1.0.4") {
			switch(lib,
				   "STRINGDB" = paste0("lib/", taxid, "/", taxid, ".string.interactions.txt.gz"),
				   "REACTOMEDB" = paste0("lib/", taxid, "/", taxid, "_REACTOME_low.tsv.gz"),
				   "ORGDB" = if(taxid == 9606) {
				   	"EnsDb.Hsapiens.v86"
				   } else if(taxid == 10090) {
				   	"EnsDb.Mmusculus.v79"
				   } else if(taxid == 10116) { 
				   	"EnsDb.Rnorvegicus.v79"
				   })
		} else {
			switch(lib,
				   "STRINGDB" = paste0("lib/", taxid, "/", taxid, ".string.interactions.txt.gz"),
				   "REACTOMEDB" = paste0("lib/", taxid, "/", taxid, "_REACTOME_low.tsv.gz"),
				   "ORGDB" = if(taxid == 9606) {
				   	"EnsDb.Hsapiens.v86"
				   } else if(taxid == 10090) {
				   	"EnsDb.Mmusculus.v79"
				   } else if(taxid == 10116) { 
				   	"EnsDb.Rnorvegicus.v79"
				   })
		}
	}
	
}

fpaths <- function(s) {
	switch(s,
		   "1" = "data/donors.uniprot.csv",
		   "2" = "data/E-GEOD-1675-A-AFFY-18-normalized-expressions.txt",
		   "3" = "data/E-GEOD-14226-A-AFFY-23-normalized-expressions.txt"
	)
}