
# Create mock data

set.seed(100)

# Major disease path 3 - 3, 5, 150
# Unrelated path     3 - 3, 15, 500
# Minor disease path 4 - 5, 12, 800
# Major disease path 3 - 5, 7, 60

mock = as.data.frame(matrix(runif(600, min=4e2, max=6e2), ncol=6))
colnames(mock) <- c("Control_1", "Control_2", "Control_3", "GS_1", "GS_2", "GS_3")

reactome <- data.table::fread('~/Large_files/reactome/reactome.uniprot.levels.csv', sep = '\t', header = F, data.table = F, strip.white=TRUE)
reactome <- reactome[,2:ncol(reactome)]
colnames(reactome) <- c("UniprotID", "ReactomeID", "URL", "Pathway_name", "Evidence_code", "Species", "Top", "TopPathway")

tissues <- data.table::fread('~/Large_files/tissues.uniprot.csv', sep = '\t', header = T, data.table = F, strip.white=TRUE)

random_proteins <- unique(reactome$UniprotID[!reactome$UniprotID %in% tissues$epidermis])
rownames(mock) <- sample(random_proteins, 100, replace = F)

other <- reactome$UniprotID[!reactome$UniprotID %in% tissues$epidermis]
other <- other[!other %in% rownames(mock)]

# Sample 3 majdp proteins from own dataset
majp1 <- sample(rownames(mock), 3, replace = F)
majp1tmp <- majp1
majp1 <- c(majp1, sample(reactome$UniprotID[reactome$UniprotID %in% tissues$epidermis], 2, replace = F))
majp1 <- c(majp1, sample(other, 200, replace = F))

unrel <- sample(rownames(mock), 3, replace = F)
unreltmp <- unrel
unrel <- c(unrel, sample(reactome$UniprotID[reactome$UniprotID %in% tissues$epidermis], 12, replace = F))
unrel <- c(unrel, sample(other, 500, replace = F))

minp1 <- sample(rownames(mock), 4, replace = F)
minp1tmp <- minp1
minp1 <- c(minp1, sample(reactome$UniprotID[reactome$UniprotID %in% tissues$epidermis], 8, replace = F))
minp1 <- c(minp1, sample(other, 785, replace = F))

majp2 <- sample(rownames(mock), 3, replace = F)
majp2tmp <- majp2
majp2 <- c(majp2, sample(reactome$UniprotID[reactome$UniprotID %in% tissues$epidermis], 4, replace = F))
majp2 <- c(majp2, sample(other, 53, replace = F))

#mock[1:7,4:6] <- mock[1:7,4:6] * 4

mock[rownames(mock) %in% majp1tmp,4:6] <- mock[rownames(mock) %in% majp1tmp,4:6] * 10
mock[rownames(mock) %in% minp1tmp,4:6] <- mock[rownames(mock) %in% minp1tmp,4:6] * 6
mock[rownames(mock) %in% unreltmp,4:6] <- mock[rownames(mock) %in% unreltmp,4:6] * 12
mock[rownames(mock) %in% majp2tmp,4:6] <- mock[rownames(mock) %in% majp2tmp,4:6] * 16

# ABC -> P1, not in epi
# 2 more -> P1, in epi
# 145 more, not in epi, not in df

#mydata <- read.csv("mock.data.csv", row.names = 1)

write.csv(mock, file = "mock.data.csv", quote = F, row.names=T)
write.csv(majp1, file = "majp1.data.csv", quote = F, row.names=F)
write.csv(unrel, file = "unrel.data.csv", quote = F, row.names=F)
write.csv(minp1, file = "minp1.data.csv", quote = F, row.names=F)
write.csv(majp2, file = "majp2.data.csv", quote = F, row.names=F)



#View(reactome)
#p1 <- reactome[reactome$Pathway_name == "Activated NTRK2 signals through FYN",]
#p2 <- reactome[reactome$Pathway_name == "Activation of C3 and C5",]
#p3 <- reactome[reactome$Pathway_name == "Signaling by Activin",]


# We have 600 randomly selected proteins
# We want to assign ~ 30 of these, selected from skin related 
# There are roughly 2000 skin related proteins in total
#length(unique(tissues$skin[!is.na(tissues$skin)]))

disease <- sample(tissues$epidermis[!is.na(tissues$epidermis)], 10, replace = F)
rn <- rownames(mock)
rn[1:5] <- sample(disease, 5, replace = F)
rownames(mock) <- rn


disease2 <- sample(unique(reactome$UniprotID), 1000, replace = F)
disease2 <- disease2[1:900]

# proteins <- gsub("-.*", "", c(p1$UniprotID, p2$UniprotID, p3$UniprotID))
# rn <- rownames(mock)
# rn[1:25] <- proteins
# rownames(mock) <- rn
# 
# rownames(mock) %in% proteins

# Increase expression in disease
mock[1:4,4:6] <- mock[1:4,4:6] * 50
mock[5:8,4:6] <- mock[5:8,4:6] * 4

write.csv(mock, file = "mock.data.csv", quote = F, row.names=T)
write.csv(disease, file = "disease.csv", quote = F, row.names=F)
write.csv(disease2, file = "additional.csv", quote = F, row.names=F)
