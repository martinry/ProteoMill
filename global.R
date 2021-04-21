require(data.table)

source("bin/All_Classes.R")
source("bin/obsolete.R")
source("bin/ora.R")

# Upload dataset ----

convertColumns <- c("UNIPROTID", "ENTREZID", "SYMBOL", "GENEID", "PROTEINID")
assign("convertColumns", convertColumns, envir = .GlobalEnv)

