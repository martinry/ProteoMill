
require(networkD3)
require(igraph)


# STRING interactions file ----
interactions <- data.table::fread('~/Large_files/9606.protein.links.v11.0.txt', sep = ' ', header = T, data.table = F)

uniprot_to_string <- mapping[mapping$ID == 'STRING', c('Accession', 'Name')]

uniprot_to_string <- uniprot_to_string[uniprot_to_string$Accession %in% rownames(data_wide_pre),]

interactions2 <- interactions[(interactions$protein1 %in% uniprot_to_string$Name) & (interactions$protein2 %in% uniprot_to_string$Name),]



merged <- merge(uniprot_to_string, interactions2, by.x = 'Name', by.y = 'protein1')

merged2 <- merge(uniprot_to_string, merged, by.x = 'Name', by.y = 'protein2')

merged2 <- merged2[,c(2,4,5)]



merged2[,1] <- qob::switch.items(merged2[,1], mapped_genes, 1, 3)
merged2[,2] <- qob::switch.items(merged2[,2], mapped_genes, 1, 3)


