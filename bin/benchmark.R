
# Create mock data

set.seed(100)

mock = as.data.frame(matrix(runif(600, min=4e2, max=6e2), ncol=6))
colnames(mock) <- c("Control_1", "Control_2", "Control_3", "GS_1", "GS_2", "GS_3")
random_proteins <- unique(reactome$UniprotID[!reactome$UniprotID %in% tissues$epidermis])
rownames(mock) <- gsub("-.*", "", sample(random_proteins, 100, replace = F))

#View(reactome)
#p1 <- reactome[reactome$Pathway_name == "Activated NTRK2 signals through FYN",]
#p2 <- reactome[reactome$Pathway_name == "Activation of C3 and C5",]
#p3 <- reactome[reactome$Pathway_name == "Signaling by Activin",]


# We have 600 randomly selected proteins
# We want to assign ~ 30 of these, selected from skin related 
# There are roughly 2000 skin related proteins in total
#length(unique(tissues$skin[!is.na(tissues$skin)]))

disease <- sample(tissues$epidermis[!is.na(tissues$epidermis)], 100, replace = F)
rn <- rownames(mock)
rn[1:5] <- sample(disease, 5, replace = F)
rownames(mock) <- rn

# proteins <- gsub("-.*", "", c(p1$UniprotID, p2$UniprotID, p3$UniprotID))
# rn <- rownames(mock)
# rn[1:25] <- proteins
# rownames(mock) <- rn
# 
# rownames(mock) %in% proteins

# Increase expression in disease
mock[1:5,4:6] <- mock[1:5,4:6] * 50
mock[6:20,4:6] <- mock[6:20,4:6] * 4

write.csv(mock, file = "mock.data.csv", quote = F, row.names=T)
write.csv(disease, file = "disease.csv", quote = F, row.names=F)
