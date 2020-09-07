library(data.table)

source("bin/All_Classes.R")
source("bin/hierarchy.R")
source("bin/obsolete.R")
source("bin/ora.R")
source("bin/volcano.R")

# Upload dataset ----

undup <- function(genes){
    genes[!is.na(genes)]
}




# Interaction data ----

if(!exists("interactions")){
 
 interactions <- data.table::fread("lib/interactions5.txt.gz")
 assign("interactions", interactions, envir = .GlobalEnv)
 
 pdesc <- data.table::fread("lib/protein_descriptions.txt.gz")
 assign("pdesc", pdesc, envir = .GlobalEnv)
 
}
