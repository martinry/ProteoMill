
history <- function(accession, index = 1){

  res <- data.frame()

  for(i in 1:length(accession)){

    print(paste('Processing', accession[i], '-', i, 'out of', length(accession)))

    xml.url <- paste('https://www.uniprot.org/uniprot/', accession[i], '.rss?version=*', sep = '')

    x <- readChar(xml.url, nchars = 600)
    y <- htmlTreeParse(x, useInternalNodes = T)

    descriptions <- XML::xpathSApply(y,'//item/description',xmlValue)
    pubdates <- XML::xpathSApply(y,'//item/pubdate',xmlValue)

    res[i, "Accession"] <- accession[i]
    res[i, "Event"] <- gsub("</a>", "", gsub("<([[:alpha:]][[:alnum:]]*)(.[^>]*)>([.^<]*)", "", descriptions[1]))
    res[i, "Date"] <- readr::parse_date(pubdates[1], "%a, %d %B %Y")


  }

  return(res)

}
