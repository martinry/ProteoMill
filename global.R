library(igraph)
library(visNetwork)

#setwd("C://Users/martinry/qodb-shiny/")


# Upload dataset ----

undup <- function(genes){
    genes[!is.na(genes)]
}




# Interaction data ----

if(!exists("interactions")){
 
 interactions <- data.table::fread("C://Users/martinry/interactions5.txt")
 assign("interactions", interactions, envir = .GlobalEnv)
 
 pdesc <- data.table::fread("C://Users/martinry/protein_descriptions.txt")
 assign("pdesc", pdesc, envir = .GlobalEnv)
 
}
