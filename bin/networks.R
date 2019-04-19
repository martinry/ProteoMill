
require(networkD3)
require(igraph)

# STRING interactions file ----
interactions <- data.table::fread('~/Large_files/9606.protein.links.v11.0.txt', sep = ' ', header = T, data.table = F)
actions <- data.table::fread('~/Large_files/9606.protein.actions.v11.0.txt', sep = '\t', header = T, data.table = F)
mapping <- data.table::fread('~/Large_files/HUMAN_9606_idmapping.dat', sep = '\t', header = F, data.table = F)
colnames(mapping) <- c('Accession', 'ID', 'Name')
uniprot_to_string <- mapping[mapping$ID == 'STRING', c('Accession', 'Name')]

uniprot_to_string <- uniprot_to_string[uniprot_to_string$Accession %in% rownames(data_wide),]

interactions2 <- interactions[(interactions$protein1 %in% uniprot_to_string$Name) & (interactions$protein2 %in% uniprot_to_string$Name),]

interactions3 <- interactions2

interactions3$protein1 <- qob::switch.items(interactions2$protein1, uniprot_to_string, 2, 1)
interactions3$protein2 <- qob::switch.items(interactions2$protein2, uniprot_to_string, 2, 1)
