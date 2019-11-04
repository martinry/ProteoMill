# Reactome hierarchy
hier <- data.table::fread('~/Large_files/reactome/reactome.relations.hsa.txt', sep = '\t', header = F, data.table = F, strip.white = TRUE)

lowest <- data.table::fread('~/Large_files/reactome/lowest.level.paths.HS.TAS.txt', sep = '\t', header = F, data.table = F, strip.white=TRUE)
colnames(lowest) <- c("UniprotID", "ReactomeID", "URL", "Pathway_name", "Evidence_code", "Species")
lowest <- lowest[lowest$ReactomeID %in% hier$V2,]

G <- graph.data.frame(d = data.frame(hier$V2, hier$V1), directed = T)
#plot(df.g, vertex.label = V(df.g)$name)

lowest$Top <- ''

top_level <- function(id) {
    SUB = induced_subgraph(G, subcomponent(G, id, mode="out"))
    top = farthest.nodes(SUB)$vertices[2]$name
    return(top)
}

seen <- c()
for(i in 1:nrow(lowest)) {
    id <- lowest[i,"ReactomeID"]
    if(!id %in% seen) {
        print(paste(i, "out of", nrow(lowest), sep = ' '))
        top <- top_level(lowest[i,"ReactomeID"])
        lowest[lowest$ReactomeID == id, "Top"] <- top
        seen <- c(seen, id)
    }
    
}

toppaths <- qob::switch.items(lowest$Top, reactome, 2, 4, no_targets = "NA", multiple_targets = "first")
lowest$TopPathways <- toppaths
# 
write.table(file = "reactome.uniprot.lowest.levels.csv", lowest, sep = '\t', quote = F)
