

require(networkD3)


results$regulation <- ifelse(results$logFC > 0, "Up-regulated", "Down-regulated")

results.sign <- results[results$Rep.Path.Score > 0,]



paths <-
    data.frame(
        "Gene" = results.sign$Gene,
        "Regulation" = results.sign$regulation,
        "TopPath" = results.sign$Rep.Path.Top,
        "PathName" = results.sign$Rep.Path.Name
    )

paths2 <- paths[,c(1,3,4)]
fr1 <- as.data.frame(table(paths[,2:3]))
fr2 <- as.data.frame(table(paths2[,2:3]))

nodes <- data.frame("name" = c("Up-regulated", "Down-regulated", levels(paths$TopPath), levels(paths$PathName)))
links <- data.frame("source" = c(), "target" = c(), "value" = c())

ntop <- 3+length(levels(paths$TopPath))
nname <- ntop + length(levels(paths$PathName)) - 1

for(i in 3:ntop) {
    n <- as.character(nodes[i,])
    k <- i-1
    
    if(n %in% fr1$TopPath) {
        up <- fr1[fr1$Regulation == "Up-regulated",]
        if(up[up$TopPath == n,3] > 0){
            links <- rbind(links, c(0, k, up[up$TopPath == n,3]))
        }
        down <- fr1[fr1$Regulation == "Down-regulated",]
        if(down[down$TopPath == n,3] > 0){
            links <- rbind(links, c(1, k, down[down$TopPath == n,3]))
        }
    }

    if(n %in% fr2$TopPath) {
        for(j in (ntop):nname){
            m <- as.character(nodes[j,])
            l <- j-1
            
            tp <- fr2[fr2$TopPath == n & fr2$PathName == m & fr2$Freq > 0,]
            
            if(nrow(tp) > 0) {
                links <- rbind(links, c(k, l, tp[,3]))
            }
        }
    }
}

colnames(links) <- c("source", "target", "value")

P <- list("nodes" = nodes, "links" = links)






