# Load packages ----

# Annotation
library(parallel)
library(BiocGenerics)
library(Biobase)
library(stats4)
library(IRanges)
library(S4Vectors)
library(AnnotationDbi)
library(GenomeInfoDb)
library(GenomicRanges)
# library(GenomicFeatures)
# library(AnnotationFilter)
# library(ensembldb)

# Plotting/UI
library(heatmaply)
require(ggplot2)
require(ggrepel)
require(RColorBrewer)
require(factoextra)
# require(fitdistrplus)
require(plotly)
require(networkD3)
require(igraph)
require(shinycssloaders)
require(shinyjs)
require(shinyBS)

# Statistics
require(limma)
require(DESeq2)
require(mixOmics)

# Parsing and reshaping
require(stringr)
require(dplyr)
require(data.table)
require(XML)
require(rmarkdown)
require(R.utils)
require(knitr)



# Generic functions ----

separator <- function(s) {
    switch(s,
           "1" = "auto",
           "2" = ",",
           "3" = ";",
           "4" = "\t")
}

identifier <- function(i) {
    switch(i,
           "1" = "auto",
           "2" = "UNIPROTID",
           "3" = "ENTREZID",
           "4" = "SYMBOL",
           "5" = "GENEID",
           "6" = "PROTEINID")
}

get_delim <- function(c){
    if(grepl("\n", c)){
        return("\n")
    } else if(grepl(",", c)){
        return(",")
    } else if(grepl(";", c)){
        return(";")
    } else {
        return ("None")
    }
}

# org_lib <- function(s) {
#     switch(s,
#            "Homo sapiens|HSA|9606" = EnsDb("lib/9606/Homo_sapiens.GRCh38.103.sqlite"),
#            "Bos taurus|BTA|9913" = EnsDb("lib/9913/Bos_taurus_hybrid.UOA_Angus_1.103.sqlite"),
#            "Caenorhabditis elegans|CEL|6239" = EnsDb("lib/6239/Caenorhabditis_elegans.WBcel235.103.sqlite"),
#            "Danio rerio|DRE|7955" = EnsDb("lib/7955/Danio_rerio.GRCz11.103.sqlite"),
#            "Drosophila melanogaster|DME|7227" = EnsDb("lib/7227/Drosophila_melanogaster.BDGP6.32.103.sqlite"),
#            "Gallus gallus|GGA|9031" = EnsDb("lib/9031/Gallus_gallus.GRCg6a.103.sqlite"),
#            "Mus musculus|MMU|10090" = EnsDb("lib/10090/Mus_musculus_wsbeij.WSB_EiJ_v1.103.sqlite"),
#            "Rattus norvegicus|RNO|10116" = EnsDb("lib/10116/Rattus_norvegicus.Rnor_6.0.103.sqlite"),
#            "Saccharomyces cerevisiae|SCE|4932" = EnsDb("lib/4932/Saccharomyces_cerevisiae.R64-1-1.103.sqlite"),
#            "Sus scrofa|SSC|9823" = EnsDb("lib/9823/Sus_scrofa_wuzhishan.minipig_v1.0.103.sqlite"),
#            "Xenopus tropicalis|XTR|8364" = EnsDb("lib/8364/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.103.sqlite")
#     )
# }

# Keep only in global.R?
convertColumns <- c("UNIPROTID", "ENTREZID", "SYMBOL", "GENEID", "PROTEINID")
assign("convertColumns", convertColumns, envir = .GlobalEnv)


# Server ----
server <- function(session, input, output) {
    
    options(shiny.maxRequestSize=30*1024^2) 
    
    # Remove text "Loading packages, please wait..."
    removeUI(selector = "#notifications-wrapper > span")
    
    # Define reactive variables ----
    notifications <- reactiveValues()
    tasks <- reactiveValues()
    sampleinfo <- reactiveValues()
    maindata <- reactiveValues()
    app_meta <- reactiveValues(palette = "Accent")
    rcont <- reactiveValues()
    pathways <- reactiveValues()
    plots <- reactiveValues()
    
    sampleinfo$sID <- identifier(2) # Set UNIPROTID as default
    
    # Notifications ----
    
    # Initialize as empty list
    # Populate with updateNotifications() function
    notification_list <- list() 
    
    tasksdf <- data.table(
        text = c("Upload a dataset",
                 "Upload annotation data",
                 "Set a filter",
                 "Inspect data",
                 "Set contrast",
                 "Run pathway enrichment",
                 "Run network analysis"),
        value = rep(0, 7),
        color = rep("green", 7),
        id = sprintf("%04d", seq(1:7))
    )
    
    tasks$tasks <- tasksdf
    
    # Make notification_list global var
    notifications$notification_list <- notification_list
    
    
    # Initialize notification menu
    output$notifMenu <- renderMenu({
        dropdownMenu(type = "notifications", .list = notifications$notification_list)
    })
    
    # Initialize notification menu
    output$helpMenu <- renderMenu({
        dropdownMenu(type = "tasks", .list = list(taskItem(text = tasks$tasks[1, text],
                                                           value = tasks$tasks[1, value],
                                                           color = tasks$tasks[1, color])))
    })
    # A function to append task items to the menu
    updateTasks <- function(text, value, color, i) {
        
        i = sprintf("%04d", i)
        
        value <- ifelse(value > 99, 100, round(value, digits = 0))
        
        l = list(text = text, value = value, color = color, id = i)
        
        if(i %in% tasks$tasks$id){
            tasks$tasks[id == i, names(tasks$tasks) := l][]
        } else {
            # Create a new item
            
            tasks$tasks <- rbindlist(list(tasks$tasks, l))
            
        }
        
        tasks$task_list <- list()
        
        tmp <- tasks$tasks[value == 0]
        tmp <- tmp[order(id)][1:2]
        tmp <- tmp[!is.na(id)]
        tmp <- rbindlist(list(tasks$tasks[value > 0], tmp))
        
        
        for(x in 1:nrow(tmp)){
            item <- taskItem(text = tasks$tasks[x, text],
                             value = tasks$tasks[x, value],
                             color = tasks$tasks[x, color])
            
            # A hack to make the notification item clickable
            # Onclick opens a modal dialog
            item$children[[1]] <- a(
                "onclick" = paste0("clickFunction('", paste0(tasks$tasks[x, id]), "'); return false;"),
                list(item$children[[1]]$children))
            
            
            item <- list(item)
            
            tasks$task_list <- append(item, tasks$task_list)
        }
        
        # Render the menu with appended notification list
        output$helpMenu <- renderMenu({
            
            dropdownMenu(type = "tasks", .list = tasks$task_list)
        })
        
    }
    
    observeEvent("", {
        updateTasks(text = "Upload a dataset", value = 0, color = "green", i = 0001)
    }, ignoreNULL = T, ignoreInit = F, once = T)
    
    
    # A function to append notification items to the menu
    updateNotifications <- function(notif, icon, status, clickable = F) {
        
        if (clickable) {
            
            item <- notificationItem(
                text = notif,
                icon = shiny::icon(icon),
                status = status,
                paste0("noti")
            )
            
            # A hack to make the notification item clickable
            # Onclick opens a modal dialog
            item$children[[1]] <- a(
                href = "#shiny-tab-dashboard", # Remove this then test
                "onclick" = paste0("clickFunction('", paste0(notif), "'); return false;"),
                list(item$children[[1]]$children)
            )
        } else {
            item <- notificationItem(
                text = notif,
                icon = shiny::icon(icon),
                status = status
            )
        }
        
        item <- list(item)
        
        notifications$notification_list <- append(item, notifications$notification_list)
        
        # Render the menu with appended notification list
        output$notifMenu <- renderMenu({
            
            dropdownMenu(type = "notifications", .list = notifications$notification_list)
        })
        
        
        # Animate the message next to the drop down menu with jquery
        setMessage <- function() {
            return(notif)
        }
        
        setIcon <- function() {
            return(status)
        }
        
        session$sendCustomMessage("background-color", setMessage())
        session$sendCustomMessage("set-icon", setIcon())
        
    }
    
    #### ANALYSIS ####
    # Sidebar menu: rules ----
    
    observeEvent(input$sidebarmenu, {
        
        if(input$sidebarmenu %in% c("structures")){
            updateNotifications("This feature is not yet available.","exclamation-triangle", "danger")
        }
        
        if(input$sidebarmenu == "blast"){
            updateNotifications("This feature is not yet available.","exclamation-triangle", "danger")
        }
        
        if(input$sidebarmenu == "news"){
            showModal(
                modalDialog(size = "l",
                            easyClose = T,
                            renderUI({
                                tags$iframe(
                                    src = paste0("doc/News.html"),
                                    width = "100%",
                                    height = "600px",
                                    frameborder = "0")})))
        }
        
        
        if(input$sidebarmenu == "interactions"){
            if(is.null(pathways$UPREGULATED_pathways)){
                updateNotifications("Run pathway analysis first.","exclamation-triangle", "danger")
            } else{
                updateTasks(text = "Run network analysis", value = 100, color = "green", i = 0007)
            }
        }
        
    })
    
    # Clicked notification item
    observeEvent(input$linkClicked, {
        
        linkCode <- as.character(input$linkClicked)
        
        i <- as.character(substr(linkCode, nchar(linkCode)-3, nchar(linkCode)))
        
        item <- tasks$tasks[id == i]
        
        showModal(
            modalDialog(title = paste0("Task: ", item$text),
                        h6(paste0("Progress: ", item$value, "%")),
                        renderUI({
                            tagList(tags$iframe(src = paste0("doc/", gsub(" ", "-", item$text), ".html"), width="100%", height="600px", frameborder="0"),
                                    if(i == "0001"){tags$iframe(src = paste0("https://qodb.shinyapps.io/generateDataset/"), width="100%", height="600px", frameborder="0")})
                            
                        }),
                        size = "l",
                        easyClose = T,
                        footer = list(modalButton("Dismiss"))
            )
        )
    })
    
    # Sidebar menu: items ----
    
    # Render "locked" menu items
    output$qualityrm <- renderMenu({ menuItem("Inspect data", icon = icon("lock"), tabName = "") })
    output$diffrm    <- renderMenu({ menuItem("Differential analysis", icon = icon("lock"), tabName = "") })
    output$enrichrm  <- renderMenu({ menuItem("Pathway analysis", icon = icon("lock"), tabName = "") })
    output$networkrm <- renderMenu({ menuItem("Network analysis", icon = icon("lock"), tabName = "") })
    
    unlock_menus <- function() {
        output$quality <- renderMenu({
            menuItem("Inspect data", icon = icon("check-circle"), href = NULL,
                     menuSubItem("PCA", tabName = "PCA", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F),
                     menuSubItem("Heatmap", tabName = "samplecorr", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F)
            )
        })
        output$differential <- renderMenu({
            menuItem("Differential analysis", class = 'btn-10', tabName = "differential", icon = icon("adjust"), href = NULL,
                     menuSubItem("Set contrasts", tabName = "contrasts", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F),
                     menuSubItem("Differential expression", tabName = "differentialexpression", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F)
            )
        })
        
        removeUI(selector = "#qualityrm")
        removeUI(selector = "#diffrm")
        
        updateSelectInput(session, "contrast1", choices = sampleinfo$group)
        
    }
    
    observeEvent(input$displayIdentifier, {
        
        i <- identifier(input$displayIdentifier)
        
        sampleinfo$sID <- i
        
        updateNotifications(paste0("Default ID set to ", i),"info-circle", "info")
        
    })
    
    
    
    
    
    # Settings: toggle tooltips
    
    observeEvent(input$toggleToolTip, {
        
        if(input$toggleToolTip == "Yes") shinyjs::runjs("$('*[title]').tooltip('enable');")
        else shinyjs::runjs("$('*[title]').tooltip('disable');")
        
        
        
    })
    
    # Build sample info ----
    
    sample_data <- function(d) {
        samples <- names(d)
        treatment <- as.factor(gsub('_.*', '', samples))
        replicate <- as.factor(gsub('.*_', '', samples))
        
        # We use data.frame instead of data.table here because Biobase requires it
        samples <- data.frame(samples, treatment, replicate)
        rownames(samples) <- samples$samples
        
        group <- sapply(as.character(unique(treatment)), function(x) paste("treatment", x, sep = ''))
        
        sampleinfo$samples <- samples
        # sampleinfo$treatment <- treatment
        # sampleinfo$replicate <- replicate
        
        sampleinfo$group <- group
        
        updateContrasts()
        
        return (FALSE)
    }
    
    subsample_data <- function(){
        samples <- sampleinfo$samples
        treatment <- samples$treatment
        
        group <- sapply(as.character(unique(treatment)), function(x) paste("treatment", x, sep = ''))
        sampleinfo$group <- group
    }
    
    
    # Filter NA ----
    
    # Remove rows where missing values per treatment is greater than set threshold
    subset_by_na <- function(dataset, treatments, threshold = 1) {
        
        dataset <- copy(dataset)
        
        # Get "treatment" names
        cols <- unique(treatments)
        
        # For each treatment, count missing values
        for(i in seq_along(cols)) {
            cn <- startsWith(names(dataset), as.character(cols[i]))
            dataset[, paste0("NAcount", i) := Reduce(`+`, lapply(.SD, function(x) is.na(x))), .SDcols = cn]
        }
        
        # Columns that keep NA sums
        NAsums <- dataset[, startsWith(names(dataset), "NAcount"), with = F]
        
        # Get treatment with greatest number of NAs
        dataset$NAcountMax <- apply(NAsums, 1, max)
        
        # Keep only rows where treatment with max NA is <= 1
        dataset <- dataset[NAcountMax <= threshold, -c(startsWith(names(dataset), "NAcount")), with = F]
        
        return(dataset)
        
    }
    
    subset_interactions <- function() {
        
        all_proteins <- maindata$udat@identifiers[, "UNIPROTID"][[1]]
        
        interactions <- data.table::fread(paste0("lib/", maindata$taxid, "/", maindata$taxid, ".string.interactions.txt.gz"))
        
        pathways$ints <- interactions[(Interactor1 %in% all_proteins) & (Interactor2 %in% all_proteins)]
        
        N <- maindata$udat@main[, .N]
       
        if(dplyr::between(N, 1500, 2999)) updateNumericInput(session = session, "interactionConfidence", label = "Interaction confidence", min = 0, max = 10, value = 7, step = 0.1)
        else if(dplyr::between(N, 3000, 4999)) updateNumericInput(session = session, "interactionConfidence", label = "Interaction confidence", min = 0, max = 10, value = 8, step = 0.1)
        else if(dplyr::between(N, 5000, 8999)) updateNumericInput(session = session, "interactionConfidence", label = "Interaction confidence", min = 0, max = 10, value = 9, step = 0.1)
        else if(dplyr::between(N, 9000, 14999)) updateNumericInput(session = session, "interactionConfidence", label = "Interaction confidence", min = 0, max = 10, value = 9.3, step = 0.1)
        else if(N >= 15000) updateNumericInput(session = session, "interactionConfidence", label = "Interaction confidence", min = 0, max = 10, value = 9.5, step = 0.1)
        
        
    }
    
    
    # File input: upload data function ----
    
    undup <- function(genes) {
        genes[!is.na(genes)]
    }
    
    infer_id <- function(d, k) {
        
        # Guess input ID based on successful conversions on small sample
        
        tr1 <- suppressWarnings(AnnotationDbi::mapIds(d, keys = k, column = "SYMBOL", keytype = "UNIPROTID", multiVals = "first"))
        tr2 <- suppressWarnings(AnnotationDbi::mapIds(d, keys = k, column = "UNIPROTID", keytype = "ENTREZID", multiVals = "first"))
        tr3 <- suppressWarnings(AnnotationDbi::mapIds(d, keys = k, column = "UNIPROTID", keytype = "SYMBOL", multiVals = "first"))
        tr4 <- suppressWarnings(AnnotationDbi::mapIds(d, keys = k, column = "UNIPROTID", keytype = "GENEID", multiVals = "first"))
        tr5 <- suppressWarnings(AnnotationDbi::mapIds(d, keys = k, column = "UNIPROTID", keytype = "PROTEINID", multiVals = "first"))
        
        trs <- list(tr1[!is.na(tr1)],
                    tr2[!is.na(tr2)],
                    tr3[!is.na(tr3)],
                    tr4[!is.na(tr4)],
                    tr5[!is.na(tr5)])
        
        # Our guess is whichever id type has the highest success rate
        i <- convertColumns[which.max(lapply(lapply(trs, lengths), sum))]
        
    }
    
    upload_data <- function(i = "auto") {
        
        data_wide <- maindata$data_wide_tmp
        
        data_wide <- data_wide[!duplicated(names(data_wide)[1])]
        
        # Convert all values to numeric
        for(j in seq_along(data_wide)){
            set(data_wide, i = which(data_wide[[j]] == 0 & is.numeric(data_wide[[j]])), j = j, value = NA)
        }
        
        # Remove completely empty rows
        empty_rows <- apply(data_wide[,2:ncol(data_wide)], 1, function(x) all(is.na(x)))
        data_wide <- data_wide[!empty_rows,]
        
        # We keep this variable as the originally uploaded dataset.
        # Used later to enable reverting to original from subsetted dataset.
        maindata$data_origin <- data_wide
        
        if(input$LogTransformData) data_wide[, 2:ncol(data_wide)] <- log2(data_wide[, 2:ncol(data_wide)])
        
        if(input$organism != "Other") {
            # 10 first rows for converting IDs
            keys <- data_wide[, as.character(.SD[[1L]])][1:10]
            
            i <- infer_id(d = maindata$organism, k = keys)
            
            if(i == "UNIPROTID") {
                tr_all <- data.table(AnnotationDbi::select(maindata$organism, keys = data_wide[, as.character(.SD[[1L]])], columns = convertColumns, keytype = i))
                tr_all <- tr_all[!duplicated(UNIPROTID)]
                
            } else {
                tr <- AnnotationDbi::mapIds(maindata$organism, keys = data_wide[, as.character(.SD[[1L]])], column = "UNIPROTID", keytype = i, multiVals = "first")
                tr <- tr[!is.na(tr)]
                tr <- tr[!duplicated(tr)]
                
                tr <- data.table(i = names(tr), "UNIPROTID" = unname(tr))
                
                tr_all <- data.table(AnnotationDbi::select(maindata$organism, keys = data_wide[, as.character(.SD[[1L]])], columns = convertColumns, keytype = i))
                tr_all <- tr_all[UNIPROTID %in% tr$UNIPROTID]
                tr_all <- tr_all[!duplicated(UNIPROTID)]
            }
            
            names(data_wide)[1] <- i
            
            setkeyv(tr_all, i)
            
            setkeyv(data_wide, i)
            
            tr_all <- tr_all[data_wide[, ..i], on = i]
            
            if(all(tr_all[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = 2:5] == tr_all[, .N])) all_missing <- T
            else all_missing <- F
            
            tr_all$ENTREZID <- as.character(tr_all$ENTREZID)
            
            tr_all[is.na(ENTREZID), ENTREZID := paste0("MISSING_", seq(1:length(is.na(ENTREZID))))]
            tr_all[is.na(SYMBOL), SYMBOL := paste0("MISSING_", seq(1:length(is.na(SYMBOL))))]
            tr_all[is.na(UNIPROTID), UNIPROTID := paste0("MISSING_", seq(1:length(is.na(UNIPROTID))))]
            tr_all[is.na(GENEID), GENEID := paste0("MISSING_", seq(1:length(is.na(GENEID))))]
            tr_all[is.na(PROTEINID), PROTEINID := paste0("MISSING_", seq(1:length(is.na(PROTEINID))))]
            
            pdesc <- data.table::fread(paste0("lib/", maindata$taxid, "/", maindata$taxid, ".protein.info.txt.gz"))
            
            setkey(tr_all, "UNIPROTID")
            setkey(pdesc, "UNIPROTID")
            
            pdesc <- pdesc[tr_all]
        } else {
            pdesc <- data.table()
            tr_all <- data.table("UNIPROTID" = data_wide[, 1][[1]],
                                 "ENTREZID" = "",
                                 "SYMBOL" = "",
                                 "GENEID" = "",
                                 "PROTEINID" = "")
            all_missing <- F
        }
            
       
        
        udat <- new("userdata",
                   raw            = data_wide[, 2:ncol(data_wide)],
                   main           = data_wide[, 2:ncol(data_wide)],
                   rawidentifiers = tr_all,
                   identifiers    = tr_all,
                   descriptions   = pdesc,
                   deoutput       = data.table())
        
        maindata$udat <- udat
        
        sample_data(udat@main)
        
        if(all_missing) return(T)
        else return(F)
        
    }

    # File input: Demo data ----
    observeEvent(input$useDemoData, {
        
        # Sample 1
        
        updateCheckboxInput(session = session, inputId = "LogTransformData", value = T)
        
        sample_1_exp <- "data/donors.uniprot.csv"
        maindata$data_wide_tmp <- fread(sample_1_exp, header = T)
        
        updateNotifications(paste0("Loading human annotation library..."), "info-circle", "info")
        library(EnsDb.Hsapiens.v86)
        
        maindata$organism <- EnsDb.Hsapiens.v86
        maindata$taxid <- 9606
        
        updateNotifications(paste0("Loading human interaction library..."), "info-circle", "info")
        upload_data()
        subset_interactions()
        
        app_meta$palette <- brewer.pal(uniqueN(sampleinfo$samples$treatment), "Accent")
        
        updateNotifications("Demo data uploaded.","check-circle", "success")
        updateTasks(text = "Upload a dataset", value = 100, color = "green", i = 0001)
        updateTasks(text = "Upload annotation data", value = 100, color = "green", i = 0002)
        
    })
    
    
    
    # Data summary plot ----
    
    nproteins <- reactive({
        if(!is.null(maindata$data_wide_tmp)) maindata$data_wide_tmp[, .N]
    })
    
    nsamples <- reactive({
        if(!is.null(maindata$data_wide_tmp)) (ncol(maindata$data_wide_tmp) - 1)
    })
    
    ntreatments <- reactive({
        if(!is.null(maindata$data_wide_tmp)) names(maindata$data_wide_tmp)[2:ncol(maindata$data_wide_tmp)] %>% sub('_.*', '', .) %>% uniqueN()
    })
    
    output$dataDetails <- renderUI({
        n_proteins <- paste(strong("Proteins:"), nproteins())
        n_samples <- paste(strong("Samples:"), nsamples())
        n_treatments <- paste(strong("Treatments:"), ntreatments())
        HTML(paste(n_proteins, n_samples, n_treatments, sep = "&nbsp;&nbsp;&nbsp;&nbsp;"))
    })
    
    
    output$datainfoBox <- renderInfoBox({
        if(!is.null(sampleinfo$samples)) infoBox("Proteins", maindata$udat@main[, .N])
        else infoBox("Proteins", 0)
    })
    
    output$sampleinfoBox <- renderInfoBox({
        if(!is.null(sampleinfo$samples)) infoBox("Samples", nrow(sampleinfo$samples))
        else infoBox("Samples", 0)
    })
    
    output$conditioninfoBox <- renderInfoBox({
        if(!is.null(sampleinfo$samples)) infoBox("Treatments", length(unique(sampleinfo$samples$treatment)))
        else infoBox("Treatments", 0)
    })
    
    observe({
        
        samples <- sampleinfo$samples$samples
        
        if(!is.null(samples)) {
            
            samples <- as.character(samples[order(samples)])
            
            shinyWidgets::updatePickerInput(
                session = session,
                inputId = "includesamples",
                choices = samples,
                selected = samples
            )
        }
        
    })
    
    observeEvent(input$confirmexclude, {
        
        # Check that we still have at least two conditions
        
        if(!is.null(input$includesamples)) {
            if(sum(sapply(levels(sampleinfo$samples$treatment), FUN = function(x) any(startsWith(input$includesamples, x)))) >= 2) {
                maindata$udat@main <- maindata$udat@main[, as.character(input$includesamples), with = F]
                maindata$udat@main <- maindata$udat@main[apply(maindata$udat@main, 1, FUN = function(x) !all(is.na(x))),]
                
                maindata$udat@raw <- maindata$udat@raw[, as.character(input$includesamples), with = F]
                maindata$udat@raw <- maindata$udat@raw[apply(maindata$udat@raw, 1, FUN = function(x) !all(is.na(x))),]
                
                sample_data(maindata$udat@main)
                
            } else {
                updateNotifications("At least two groups needed for comparison.","exclamation-triangle", "danger")
            }
        } else {
            updateNotifications("At least two groups needed for comparison.","exclamation-triangle", "danger")
        }
        
        
        
    })
    
    
    output$identifierinfo <- renderTable({
        if(!is.null(maindata$udat)) {
            
            missingids <- maindata$udat@identifiers
            missingidspc <- apply(missingids, 2, FUN = function(x) paste0(round(((sum(!startsWith(x, "MISSING")) / missingids[, .N]) * 100), 2), "% (", sum(!startsWith(x, "MISSING")), ")") )
            as.data.frame(missingidspc)
        }
    }, rownames = T, colnames = F)
    
    
    make_violin <- reactive({
        
        dat <- stack(as.data.frame(maindata$udat@main))
        dat <- dat[!is.na(dat$values), ]
        dat$sample <- sub("_.*", "", dat$ind)
        dat$ind <- factor(dat$ind, levels = levels(dat$ind)[order(levels(dat$ind))])
        dat$values <- dat$values
        
        return(dat)
        
    })
    
    # Remake reactive
    output$violinplot <- renderPlot({
        
        if(!is.null(maindata$udat) & !is.null(sampleinfo$samples)){
            
            dat <- make_violin()
            
            if(input$violintype == 1) {
                
                ggplot(dat, aes(x = sample, y = values, fill = sample)) + 
                    geom_violin(trim = FALSE, scale = "width") +
                    geom_boxplot(width=0.1, fill="white") +
                    ggthemes::theme_clean() +
                    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.position = "none") +
                    xlab("") + ylab("Log2 Expression") +
                    scale_fill_brewer(palette = "BuGn")
                
            } else {
                
                ggplot(dat, aes(x = ind, y = values, fill = sample)) + 
                    geom_violin(trim = FALSE, scale = "width") +
                    geom_boxplot(width=0.1, fill="white") +
                    ggthemes::theme_clean() +
                    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.position = "none") +
                    xlab("") + ylab("Log2 Expression") +
                    scale_fill_brewer(palette = "BuGn")
            }
            
        }

    })
    
    # Missing values: frequency plot ----
    
    nafrequencies <- reactive({
            
        data_wide <- maindata$udat@main
        
        # Which elements are NA?
        allNA <- is.na(data_wide)
        
        # Summary of how many TRUEs there are in each row
        NA_frequency <- table(rowSums(allNA))
        
        naf <- as.data.frame(NA_frequency)
        
        
    })
    
    output$nafreq <- renderPlot({
        
        if(!is.null(maindata$udat)){
            naf <- nafrequencies()
            
            # Draw plot
            ggplot(naf, aes(x = Var1, y = Freq, color = Var1)) +
                geom_bar(stat = "identity", width = .9, fill = "white") +
                labs(x = "Number of missing values in at least one sample", y = "Number of proteins") +
                ggthemes::theme_clean() +
                theme(axis.text.x = element_text(vjust = 0.6), legend.position = "none") +
                scale_color_grey()
        
        }
        
    })
    
    observeEvent(input$setcutoff, {
        
        if(!is.null(maindata$udat)) {
            data_wide <- maindata$udat@main
            
            subsample_data()
            
            d <- cbind(maindata$udat@rawidentifiers, maindata$udat@raw)
            data_wide <- subset_by_na(dataset = d, treatment = sampleinfo$samples$treatment, threshold = input$missingvalues)
            
            ids <- data_wide[, 1:5]
            data_wide <- data_wide[, 6:ncol(data_wide)]
            
            unlock_menus()
            
            
            maindata$udat@main <- data_wide
            maindata$udat@identifiers <- ids

            updateTasks(text = "Set a filter", value = 100, color = "green", i = 0003)
            updateNotifications(paste0("NA cutoff set to ", input$missingvalues, ".") ,"check-circle", "success")
        } else {
            updateNotifications("Upload a dataset first.","exclamation-triangle", "danger")
        }
        
    })
    
    
    
    # Differential expression: set contrasts ----
    observeEvent(input$contrast1, {
        subsample_data()
        updateContrasts()
    })
    
    updateContrasts <- function(){
        
        group <- sampleinfo$group
        
        if (is.null(group)){
            return(NULL)
        } else {
            
            if(input$contrast1 %in% group) {
                cont1 <- input$contrast1
            } else {
                cont1 <- group[1]
            }
            
            cont2 <- group[group != cont1]
            
            updateSelectInput(session, "contrast2", choices = cont2)
        }
        
    }
    
    observeEvent(input$setContrast, {
        if(!input$organism == "Other") {
            output$enrichment <- renderMenu({
                menuItem("Enrichment analysis", icon = icon("flask"), href = NULL,
                         menuSubItem("Pathway enrichment", tabName = "pathwayenrichment", href = NULL, newtab = TRUE,
                                     icon = shiny::icon("angle-double-right"), selected = F),
                         menuSubItem("Pathway visualization", tabName = "pathwayvisualization", href = NULL, newtab = TRUE,
                                     icon = shiny::icon("angle-double-right"), selected = F)
                )
            })
            output$network <- renderMenu({
                menuItem("Network analysis", icon = icon("connectdevelop"),
                         menuSubItem("Interactions", tabName = "interactions", href = NULL, newtab = TRUE,
                                     icon = shiny::icon("angle-double-right"), selected = F))
            })
            

            
            removeUI(selector = "#enrichrm")
            removeUI(selector = "#networkrm")
            removeUI(selector = "#interactionsrm")
        }
        
        pairing <- input$diffexppairing
        
        contrasts <- paste(input$contrast1,input$contrast2, sep = '-')
        
        contrast <- diff_exp(contrasts, pairing)
        
        rcont$contrasts <- contrasts
        rcont$contrast <- contrast
        
        updateNotifications("A DE contrast has been set.","check-circle", "success")
        updateTasks(text = "Set contrast", value = 100, color = "green", i = 0005)   
        
    })
    
    observeEvent(input$diffexppairing, {
        
        if(input$diffexppairing == 1) {
            # Paired analysis can only be performed when at least one sample is repeated.
            sinfo <- sampleinfo$samples
            sinfo <- as.data.table(sinfo)
            
            if(length(Reduce(setdiff, sinfo[, .(list(unique(replicate))), treatment]$V1)) > 0){
                
                updateRadioButtons(session = session, inputId = "diffexppairing", selected = 2)
                
                updateNotifications("Paired analysis can only be performed when at least one sample is repeated.","exclamation-triangle", "danger")
            }
        }
        

        
        
    }, ignoreInit = T)
    
    # Differential expression: DEA function ----
    
    diff_exp <- function(coeff, pairing) {
        
        dw <- maindata$udat@main
        sinfo <- sampleinfo$samples
        treatment <- sampleinfo$samples$treatment
        repl <- sampleinfo$samples$replicate

        
        if(input$setDEengine == 2) {
            
            isolate(updateNotifications(paste0("This process may take a while. Please wait."), "info-circle", "info"))
            
            # Replace NA -> 0
            for(j in seq_along(dw)){
                set(dw, i = which(is.na(dw[[j]]) & is.numeric(dw[[j]])), j = j, value = 0)
            }
            
            dw <- as.matrix(dw)
            
            sid <- sampleinfo$sID
            
            rn <- as.character(maindata$udat@identifiers[, ..sid][[1]])
            
            rownames(dw) <- rn
            
            if(pairing == 1) {
                # Paired
                
                design <- DESeqDataSetFromMatrix(countData  = round(2**dw),
                                                 colData    = sinfo,
                                                 design     = ~ 0 + treatment + replicate)
                
            } else {
                # Unpaired
                design <- DESeqDataSetFromMatrix(countData  = round(2**dw),
                                                 colData    = sinfo,
                                                 design     = ~ 0 + treatment)
            }
            
            dds <- DESeq(design)
            
            contrast <- results(dds, contrast=c("treatment", sub("treatment", "", input$contrast1), sub("treatment", "", input$contrast2)))
            
            colnames(contrast) <- c("baseMean", "logFC", "logFC.SE", "stat", "P.Value", "adj.P.Val")
            
            contrast <- suppressWarnings(data.table::as.data.table(contrast, keep.rownames = T))
            
            contrast <- contrast[, -"stat"]
            
            setcolorder(contrast, c("rn", "logFC", "logFC.SE", "baseMean", "P.Value", "adj.P.Val"))
            
            contrast <- contrast[order(contrast[, 1])]
            
            maindata$udat@deoutput <- contrast[, 2:ncol(contrast)]
            
            # Reduce network load
            if(contrast[, .N] > 2000) {
                # x <- unname(quantile(abs(contrast$logFC)))
                # x <- round(x[5] - x[4])
                updateNumericInput(session = session, "fccutoff", value = 2)
            }
            
        } else if(input$setDEengine == 1) {
            
            dw <- as.matrix(dw)
            
            sid <- sampleinfo$sID
            
            rn <- as.character(maindata$udat@identifiers[, ..sid][[1]])
            
            rownames(dw) <- rn
            
            # Create Annotation data and expression set (Biobase)
            phenoData <- new("AnnotatedDataFrame", data = sampleinfo$samples)
            exampleSet <- ExpressionSet(assayData = dw, phenoData = phenoData)
            
            unpaired <- model.matrix( ~ 0 + treatment )
            paired <- model.matrix( ~ 0 + treatment + repl )
            
            if(pairing == 1) {
                design <- paired
            } else {
                design <- unpaired
            }
            
            # Fit the linear model
            fit <- lmFit(exampleSet, design)

            # Decide possible contrasts
            c <- expand.grid(sampleinfo$group, sampleinfo$group)
            cc <- factor(ifelse(c$Var1 != c$Var2, paste(c$Var1, c$Var2, sep = '-'), NA ))
            cc <- cc[!is.na(cc)]
            names(cc) <- gsub('-','', gsub('treatment','',cc))
            
            cont.matrix <- makeContrasts(contrasts = cc, levels = design) # All possible contrasts
            
            # Contrast groups, run empirical bayes statistics
            fit.cont <- contrasts.fit(fit, cont.matrix)
            fit.cont <- eBayes(fit.cont, robust = T)
            
            # Generate data frame with results from linear model fit, with confidence intervals.
            contrast <- topTable(fit.cont, number = Inf, coef = coeff)
            
            contrast <- data.table::as.data.table(contrast, keep.rownames = T)
            #names(contrast) <- c(sampleinfo$sID, names(contrast[, 2:ncol(contrast)]))
            
            contrast <- contrast[order(contrast[, 1])]
            
            maindata$udat@deoutput <- contrast[, 2:ncol(contrast)]
            
            # Reduce network load
            if(contrast[, .N] > 2000) updateNumericInput(session = session, "fccutoff", value = 2)
            
        }
        
        return( contrast )
        
    }
    
    output$diffexptable_summary <- renderTable({
        contrast <- maindata$udat@deoutput
        
        if(maindata$udat@deoutput[, .N] > 0){
            contrast <- contrast[abs(logFC) >= input$diffexp_limit_fc]
            contrast <- contrast[adj.P.Val < input$diffexp_limit_pval]
            contrast <- contrast[!is.na(logFC)]
            
            df <- data.frame(as.character(contrast[, .N]),
                             as.character(contrast[logFC >= 0, .N]),
                             as.character(contrast[logFC <= 0, .N]),
                             mean(contrast[, logFC]),
                             mean(contrast[logFC >= 0, logFC]),
                             mean(contrast[logFC <= 0, logFC]),
                             min(contrast[logFC >= 0, logFC]),
                             min(contrast[logFC <= 0, logFC]),
                             max(contrast[logFC >= 0, logFC]),
                             max(contrast[logFC <= 0, logFC])
            )
            
            df <- t(df)
            
            rownames(df) <- c("DE proteins (all)",
                              "DE proteins (up-regulated)",
                              "DE proteins (down-regulated)",
                              "Mean logFC (all)",
                              "Mean logFC (up-regulated)",
                              "Mean logFC (down-regulated)",
                              "Min. logFC (up-regulated)",
                              "Min. logFC (down-regulated)",
                              "Max. logFC (up-regulated)",
                              "Max. logFC (down-regulated)")
            
            maindata$diffexptable_summary <- df
            
            maindata$diffexptable_summary
        }
        
        
    }, rownames = T, colnames = F)
    
    
    output$diffexptable_up <- DT::renderDataTable({
        
        contrast <- maindata$udat@deoutput[, -c("t", "B")]
        sID <- sampleinfo$sID
        contrast <- cbind(maindata$udat@identifiers[, ..sID], contrast)
        
        if(maindata$udat@deoutput[, .N] > 0){
            contrast <- contrast[abs(logFC) >= input$diffexp_limit_fc]
            contrast <- contrast[adj.P.Val < input$diffexp_limit_pval]
            contrast <- contrast[order(logFC, decreasing = T)]
            maindata$diffexptable_up <- DT::datatable(contrast,
                                                      options = list(order = list(2, 'desc'))) %>% 
                formatRound(columns=2:(ncol(contrast)), digits=4)
            
            maindata$diffexptable_up
            
            
        }
        
    })
    
    output$diffexptable_down <- DT::renderDataTable({
        
        contrast <- maindata$udat@deoutput[, -c("t", "B")]
        sID <- sampleinfo$sID
        contrast <- cbind(maindata$udat@identifiers[, ..sID], contrast)
        
        if(maindata$udat@deoutput[, .N] > 0){
            contrast <- contrast[abs(logFC) >= input$diffexp_limit_fc]
            contrast <- contrast[adj.P.Val < input$diffexp_limit_pval]
            contrast <- contrast[order(logFC, decreasing = F)]
            df <- DT::datatable(contrast,
                                options = list(order = list(2, 'asc'))) %>% 
                formatRound(columns=2:(ncol(contrast)), digits=4)
            
        }
        
    })
    
    
    # Pathways ----
    
    # User sets minimum fold change from diff exp output to decide which proteins to include in pathway analysis
    observeEvent(input$min_fc, {
        
        output$number_of_genes <- renderUI({
            
            contrast <- maindata$udat@deoutput
            
            up <- paste("Up-regulated:", contrast[logFC >= input$min_fc & adj.P.Val < input$min_pval, .N], "proteins", sep = " ")
            down <- paste("Down-regulated:", contrast[logFC < (input$min_fc * -1) & adj.P.Val < input$min_pval, .N], "proteins", sep = " ")
            
            HTML(paste(up, down, '<br/>', sep = '<br/>'))
            
        })
    })
    
    # List protein names when clicked pathway
    show_selected <- function(p, g){
        
        output$selected_pathway <- renderUI({
            
            HTML(paste(p, g, '<br/>', sep = '<br/>'))
            
        })
        
    }
    
    observeEvent(input$upregulated_pathways_table_rows_selected, {
        UPREGULATED_pathways <- pathways$UPREGULATED_pathways
        i = input$upregulated_pathways_table_rows_selected[1]
        up <- UPREGULATED_pathways[i,]
        p <- paste("<b>", up[,Pathway_name], "</b>")
        sID <- sampleinfo$sID
        g <- paste(udat@identifiers[UNIPROTID %in% unlist(up[, genes]), ..sID][[1]], collapse = ", ")
        
        show_selected(p, g)
        
    })
    
    observeEvent(input$downregulated_pathways_table_rows_selected, {
        DOWNREGULATED_pathways <- pathways$DOWNREGULATED_pathways
        i = input$downregulated_pathways_table_rows_selected[1]
        down <- DOWNREGULATED_pathways[i,]
        p <- paste("<b>", down[,Pathway_name], "</b>")
        #g <- paste(unlist(down[,genes]), collapse = ", ")
        sID <- sampleinfo$sID
        g <- paste(udat@identifiers[UNIPROTID %in% unlist(down[, genes]), ..sID][[1]], collapse = ", ")
        
        show_selected(p, g)
        
        
    })
    
    generate_pathways <- reactive({
        
        if(maindata$udat@deoutput[, .N] == 0){
            updateNotifications("Run differential expression first.","exclamation-triangle", "danger")
            return(F)
        } else {
            contrast <- cbind(maindata$udat@identifiers[, "UNIPROTID"], maindata$udat@deoutput)
            
            
            # In future updates, have all id types for reactome db
            if(contrast[logFC >= input$min_fc & adj.P.Val < input$min_pval & !startsWith(UNIPROTID, "MISSING"), .N] < 3 | contrast[logFC < (input$min_fc * -1) & adj.P.Val < input$min_pval & !startsWith(UNIPROTID, "MISSING"), .N] < 3) {
                isolate(updateNotifications("Too few differential proteins.","exclamation-triangle", "danger"))
                return(F)
            } else {
                
                UPREGULATED_genes <- contrast[logFC >= input$min_fc & adj.P.Val < input$min_pval, UNIPROTID]
                DOWNREGULATED_genes <- contrast[logFC < (input$min_fc * -1) & adj.P.Val < input$min_pval, UNIPROTID]
                
                UPREGULATED_pathways <- ora(UPREGULATED_genes)
                DOWNREGULATED_pathways<- ora(DOWNREGULATED_genes)
                
                if(!is.null(DOWNREGULATED_pathways) & !is.null(UPREGULATED_pathways)) {
                    UPREGULATED_pathways <- ora(UPREGULATED_genes)@output
                    pathways$UPREGULATED_pathways <- UPREGULATED_pathways
                    
                    DOWNREGULATED_pathways<- ora(DOWNREGULATED_genes)@output
                    pathways$DOWNREGULATED_pathways <- DOWNREGULATED_pathways
                    
                    isolate(updateTasks(text = "Run pathway enrichment", value = 100, color = "green", i = 0006))
                    isolate(updateNotifications("Pathway analysis complete.","check-circle", "success"))
                    
                    return(T)
                    
                } else {
                    return(F)
                }
                
            }
        }
    })
    
    output$upregulated_pathways_table <- DT::renderDT(
        
        if(maindata$udat@deoutput[, .N] > 0) {
            gp <- generate_pathways()
            if(gp){
                UPREGULATED_pathways <- pathways$UPREGULATED_pathways
                DT::datatable(UPREGULATED_pathways[, -c("genes", "background")],
                              selection = 'single',
                              options = list(autoWidth = TRUE,
                                             scrollX=TRUE,
                                             columnDefs = list(
                                                 list(width = '100px', targets = c(1, 3)),
                                                 list(width = '60px', targets = c(6, 7))
                                             )))
            }

        }
    )
    
    output$downregulated_pathways_table <- DT::renderDT(
        
        if(maindata$udat@deoutput[, .N] > 0) {
            gp <- generate_pathways()
            if(gp){
            
                DOWNREGULATED_pathways <- pathways$DOWNREGULATED_pathways
                DT::datatable(DOWNREGULATED_pathways[, -c("genes", "background")],
                              selection = 'single',
                              options = list(autoWidth = TRUE,
                                             scrollX=TRUE,
                                             columnDefs = list(
                                                 list(width = '100px', targets = c(1, 3)),
                                                 list(width = '60px', targets = c(6, 7))
                                             )))
            }
        }
    )
    

    
    # Volcano plot
    
    
    
    pathway_vis <- reactive({
        
        if(is.null(pathways$UPREGULATED_pathways)){
            isolate(updateNotifications("Run pathway analysis first.","exclamation-triangle", "danger"))
        } else {
            contrast <- maindata$udat@deoutput[, -c("t", "B")]
            #sID <- sampleinfo$sID
            contrast <- cbind(maindata$udat@identifiers, contrast)
            
            diff_df <- contrast[, c(1:6, 9)]
            colnames(diff_df)[6:7] <- c("Fold-change", "FDR")
            
            dt_all <- rbindlist(list(pathways$UPREGULATED_pathways, pathways$DOWNREGULATED_pathways))
            dt_all <- dt_all[ lengths(genes) > 0L]
            
            dt_all <- get_top_pathways(dt_all)
            
            setkey(dt_all, "UNIPROTID")
            setkey(diff_df, "UNIPROTID")
            
            diff_df <- merge.data.table(x = diff_df, y = dt_all, by = "UNIPROTID", all = T)
            
            diff_df[is.na(ReactomeID), 9:10] <- "[No significant over-representation]"
            
            diff_df[Pathway.P.Adj >= 0.05, "TopReactomeName"] <- "[No significant over-representation]"
            diff_df[Pathway.P.Adj >= 0.05, "Pathway_name"] <- "[No significant over-representation]"
            
            diff_df$Pathway_name <- stringr::str_wrap(diff_df$Pathway_name, 30)
            
            diff_df$group <- "Not Sign."
            
            diff_df[FDR < 0.05 & abs(`Fold-change`) < 1.5, "group"] <- "Sign."
            
            diff_df[FDR >= 0.05 & abs(`Fold-change`) > 1.5, "group"] <- "FC > 1.5"
            
            diff_df[FDR < 0.05 & abs(`Fold-change`) > 1.5, "group"] <- "Sign.  & FC > 1.5"
            
            diff_df$group <- as.factor(diff_df$group)
            
            diff_df$group <- factor(diff_df$group, levels = c("Sign.  & FC > 1.5", "Sign.", "FC > 1.5", "Not Sign."))
            
            return(diff_df)
            
        }
    })
    
    output$volcano_plot <- renderPlotly(
        {
            if(!is.null(pathways$UPREGULATED_pathways)){
                diff_df <- pathway_vis()
                
                nb.cols <- uniqueN(diff_df[, get(input$volcanoAnnotation)])
                mycolors <- c(colorRampPalette(brewer.pal(8, "Accent"))(nb.cols))
                
                plot_ly(data = diff_df,
                        type = "scatter",
                        x = ~`Fold-change`,
                        y = ~-log10(FDR),
                        text = ~get(sampleinfo$sID),
                        mode = "markers",
                        color = ~get(input$volcanoAnnotation),
                        colors = mycolors) %>%
                    layout(yaxis = list(type = "log1p", showgrid = T, ticks = "outside", autorange = T))
            }
            
        })
    
    
    output$volcano_plot2 <- renderPlotly(
        {
            
            if(!is.null(pathways$UPREGULATED_pathways)){
                diff_df <- pathway_vis()
                
                nb.cols <- uniqueN(diff_df[, get(input$volcanoAnnotation2)])
                mycolors <- c(colorRampPalette(brewer.pal(8, "Accent"))(nb.cols))
                
                sorted_paths <- sort(unique(diff_df[, get(input$volcanoAnnotation2)]))
                sorted_path_names <- mycolors[seq_along(sorted_paths)]
                names(sorted_path_names) <- sorted_paths
                diff_df$col <- sorted_path_names[diff_df[, get(input$volcanoAnnotation2)]]
                
                plot_ly(data = diff_df,
                        type = "scatter",
                        x = ~`Fold-change`,
                        y = ~-log10(FDR),
                        text = ~get(sampleinfo$sID),
                        key =  ~get(sampleinfo$sID),
                        mode = "markers",
                        color = ~get(input$volcanoAnnotation2),
                        colors = mycolors,
                        height = 500) %>%
                    layout(yaxis = list(type = "log1p",
                                        showgrid = T,
                                        ticks="outside",
                                        autorange = T),
                           dragmode = "lasso") %>% 
                    config(scrollZoom = T)
                    
                
            }
            
        })
    
    
    
    observeEvent(input$searchclick, {
        
        searchProtein <- function() {
            return(input$seachinput)
        }
        
        session$sendCustomMessage("searchProtein", searchProtein())
    })
    
    # Network ----
    
    # Volcano plot
    
    output$hovered_node <- renderUI({
        if (is.null(input$hovered_node)) {
            "No node has been hovered yet"
        } else {
            
            name <- input$hovered_node
            
            description <- maindata$udat@descriptions[get(sampleinfo$sID) == input$hovered_node, annotation]
            
            HTML(
                    paste(paste(tags$strong("Name:"), name, sep = " "),
                          paste(tags$strong("Description:"), description, sep = " "),
                          sep = "<br/>")
            )
            
        }
    })
    
    output$sankey <- renderSankeyNetwork({
        res <- pathways$v$res
        
        df <- pathways$v$res[Pathway_name != "[No significant over-representation]"]
        
        UPREGULATED_genes <- df[logFC >= 0 & P.Value <= 0.05]
        DOWNREGULATED_genes <- df[logFC < 0 & P.Value <= 0.05]
        
        
        df <- df[, c("Pathway_name", "TopReactomeName")]
        df <- rbindlist(list(data.table(Pathway_name = UPREGULATED_genes[, TopReactomeName], TopReactomeName = "Up-regulated"),
                             data.table(Pathway_name = DOWNREGULATED_genes[, TopReactomeName], TopReactomeName = "Down-regulated"),
                             df))
        
        links <- setDT(df)[, .N, by = c(names(df))]
        colnames(links) <- c("target", "source", "value")
        
        nodes <- data.frame(name=unique(c(links$source, links$target)))
        
        links$source <- match(links$source, nodes$name) - 1
        links$target <- match(links$target, nodes$name) - 1
        
        plots$sankey <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
                                      Target = "target", Value = "value", NodeID = "name",
                                      fontFamily = "sans-serif",
                                      fontSize = 10, nodeWidth = 60, sinksRight = F)
        
        plots$sankey
        
        
    })

    
    make_network <- reactive({
        sID <- sampleinfo$sID
        contrast <- cbind(maindata$udat@identifiers, maindata$udat@deoutput)
        
        dir <- input$network_regulation
        
        if(dir == 1){
            proteins_ <- contrast[(abs(logFC) >= input$fccutoff) &
                                     (logFC > 0)]$UNIPROTID
        } else if(dir == 2){
            proteins_ <- contrast[(abs(logFC) >= input$fccutoff) &
                                     (logFC <= 0)]$UNIPROTID
        } else {
            proteins_ <- contrast[(abs(logFC) >= input$fccutoff)]$UNIPROTID
        }
        
        if(length(proteins_) <= 1) return(NULL)
        else {
            ints2 <- pathways$ints[(Interactor1 %in% proteins_) & (Interactor2 %in% proteins_)]
            ints2 <- ints2[Score > input$interactionConfidence]
            
            if(ints2[, .N] <= 1) return(NULL)
            else {
                mynodes <- unlist(event_data("plotly_selected")$key)
                
                mynodes <- contrast[get(sID) %in% mynodes, "UNIPROTID"][[1]]
                
                interacts <- function(i){
                    return(ints2[(Interactor1 == i & Interactor2 %in% mynodes) | (Interactor2 == i & Interactor1 %in% mynodes), .N])
                }
                
                if(input$interaction_behaviour == 1){
                    
                    if(!is.null(mynodes)) {
                        if(sum(unlist(lapply(mynodes, interacts))) > 0){
                            ints2 <- ints2[Interactor1 %in% mynodes & Interactor2 %in% mynodes]
                        }
                    }
                    
                } else if(input$interaction_behaviour == 2){
                    
                    if(any(mynodes %in% ints2$Interactor1 & mynodes %in% ints2$Interactor2)){
                        ints2 <- ints2[Interactor1 %in% mynodes | Interactor2 %in% mynodes]
                    }
                    
                }
                
                nodes <- data.table(id = unique(c(ints2$Interactor1, ints2$Interactor2)))
                nodes$label <- nodes$id
                
                
                edges <-
                    data.table(
                        from = ints2$Interactor1,
                        to = ints2$Interactor2
                    )
                
                g <- igraph::graph_from_data_frame(edges, directed = F, vertices = nodes)
                
                diff_df <- pathway_vis()
                
                nb.cols <- uniqueN(diff_df[, get(input$volcanoAnnotation2)])
                mycolors <- c(colorRampPalette(brewer.pal(8, "Accent"))(nb.cols))
                
                sorted_paths <- sort(unique(diff_df[, get(input$volcanoAnnotation2)]))
                sorted_path_names <- mycolors[seq_along(sorted_paths)]
                names(sorted_path_names) <- sorted_paths
                diff_df$col <- sorted_path_names[diff_df[, get(input$volcanoAnnotation2)]]
                
                diff_df <- diff_df[V(g)$label, on = "UNIPROTID", ]
                
                g <- set.vertex.attribute(g, name = "name",  value = as.vector (diff_df[, ..sID][[1]]))
                g <- set.vertex.attribute(g, name = "color", value = diff_df$col)
                g <- set.vertex.attribute(g, name = "group", value = diff_df$TopReactomeName)
                
                if(exists("g")){
                    V(g)$name <- diff_df[V(g)$name, on = sID, ..sID][[1]]
                } else {
                    names(V(g)) <- diff_df[V(g)$name, on = "UNIPROTID", "UNIPROTID"][[1]]
                }
                
                return(g)
            }
            
        }
        
    })
    
    output$interaction_network <- renderVisNetwork({
        
        if(!is.null(pathways$UPREGULATED_pathways)){
            
            g <- make_network()
            
            if(is.null(g)) isolate(updateNotifications(paste0("No data to display. Review network options."), "info-circle", "info"))
            else {
                visNetwork::visIgraph(g) %>%
                    visInteraction(hover = TRUE) %>%
                    visIgraphLayout(layout = "layout_nicely", randomSeed = 1) %>%
                    visOptions(selectedBy = list(variable = "group"), height = "475px") %>%
                    visEvents(
                        click = "function(nodes) {
                    Shiny.onInputChange('clicked_node', {node : nodes.node});
                    ;}",
                        hoverNode = "function(nodes) {
                    console.info('hover')
                    console.info(nodes)
                    Shiny.onInputChange('hovered_node', {node : nodes.node});
                    ;}"
                    )
            }
    

        }
    })
    
    
    observeEvent(input$highlight_nodes, {
        
        if(nchar(input$network_proteins) > 0){
            
            if(get_delim(input$network_proteins) != "None") {
                
                # Multiple items
                
                selection <- unlist(strsplit(input$network_proteins, get_delim(input$network_proteins)))
                
                visNetworkProxy("interaction_network") %>%
                    visUnselectAll() %>%
                    visSelectNodes(id = selection)
                
            }
            
            else {
                # Single item
                selection <- as.character(input$network_proteins)
                
                visNetworkProxy("interaction_network") %>%
                    visUnselectAll() %>%
                    visFocus(id = selection)
            }
            
        } else {
            
            updateNotifications("Please enter a valid identifier.","exclamation-triangle", "danger")
            
        }
        
        
    })
    
    #### ADDITIONAL TOOLS ####
    # Goodness-of-fit ----
    observeEvent(input$generatedistributions, {
        
        data_origin <- maindata$udat@raw
        
        if(!is.null(data_origin)){
            
            data_wide_NAex <- na.exclude(data_origin[1:100])
            
            fit.lnorm <- tryCatch( apply(data_wide_NAex, 1, function(x) fitdistrplus::fitdist(as.numeric(x), "lnorm")),
                                   error = function(e) print(e))
            
            fit.norm <- tryCatch( apply(data_wide_NAex, 1,  function(x) fitdistrplus::fitdist(as.numeric(x), "norm")),
                                  error = function(e) print(e))
            
            # Render "data type" distribution plots
            output$distributions <- renderPlot({
                
                m <- nrow(data_wide_NAex)
                
                updateSliderInput(session, inputId = "setdist", max = m)
                
                d1 <- fit.norm
                d2 <- fit.lnorm
                
                val <- input$setdist
                
                par(mfrow = c(2, 2))
                
                fitdistrplus::denscomp(list(d1[[val]], d2[[val]]))
                fitdistrplus::qqcomp(list(d1[[val]], d2[[val]]))
                fitdistrplus::cdfcomp(list(d1[[val]], d2[[val]]))
                fitdistrplus::ppcomp(list(d1[[val]], d2[[val]]))
                
            })
            
        } else {
            
            updateNotifications("Upload a dataset first.","exclamation-triangle", "danger")
        }
        
        
        
    }, ignoreInit = TRUE)
    
    # 
    # Validate identifiers ----
    
    observeEvent(input$listCandidates, {
        
        identifiers <- maindata$udat@rawidentifiers
        
        if(!is.null(identifiers)){
            
            updateNotifications(paste0("Checking for outdated IDs. Please wait."), "info-circle", "info")
            
            candidates <- identifiers[apply(
                identifiers,
                1,
                FUN = function(x)
                    all(startsWith(x[2:5], "MISSING"))
            ), "UNIPROTID"]
            
            
            status <- history(candidates$UNIPROTID)
            status <- as.data.table(status)
            
            setkey(status, "Accession")
            setkey(candidates, "UNIPROTID")
            
            candidates <- candidates[status]
            
            candidates$flag <- ifelse(grepl("merged", candidates$Event, ignore.case = T), "Obsolete", "Up-to-date")
            
            candidates <- candidates[order(candidates$flag, decreasing = F),]
            
            n <- candidates[flag == "Obsolete", .N]
            
            updateNotifications(paste0(n, " IDs are outdated."), "info-circle", "info")
            
            output$obsolete <- DT::renderDT({
                formatStyle(
                    datatable(candidates),
                    'flag', target = 'row', 
                    backgroundColor = styleEqual(c("Obsolete", "Up-to-date"), c('snow', 'snow2'))
                ) 
                
            })
            
        } else{
            
            updateNotifications(paste0("Upload a dataset first."), "exclamation-triangle", "danger")
            
        }
        
        
        
    })
    
    
    observeEvent(input$convertids, {
        
        if(input$idstoconvert != "") {
            
            keys <- gsub("\n", ";", input$idstoconvert)
            
            keys <- strsplit(gsub("[^[:alnum:] ]", " ", keys), " +")[[1]]
            
            tr1 <- AnnotationDbi::mapIds(maindata$organism, keys = keys, column = "SYMBOL", keytype = "UNIPROTID", multiVals = "first")
            tr2 <- AnnotationDbi::mapIds(maindata$organism, keys = keys, column = "UNIPROTID", keytype = "ENTREZID", multiVals = "first")
            tr3 <- AnnotationDbi::mapIds(maindata$organism, keys = keys, column = "UNIPROTID", keytype = "SYMBOL", multiVals = "first")
            tr4 <- AnnotationDbi::mapIds(maindata$organism, keys = keys, column = "UNIPROTID", keytype = "GENEID", multiVals = "first")
            tr5 <- AnnotationDbi::mapIds(maindata$organism, keys = keys, column = "UNIPROTID", keytype = "PROTEINID", multiVals = "first")
            
            trs <- list(tr1[!is.na(tr1)],
                        tr2[!is.na(tr2)],
                        tr3[!is.na(tr3)],
                        tr4[!is.na(tr4)],
                        tr5[!is.na(tr5)])
            
            i <- convertColumns[which.max(lapply(lapply(trs, lengths), sum))]
            
            tr_all <- data.table(AnnotationDbi::select(maindata$organism, keys = keys, columns = convertColumns, keytype = i))
            
            output$convertedids <- DT::renderDT({
                tr_all
            })
            
        }
        
    })
    
    
    # Contact ----
    
    observeEvent(input$sendcomment, {
        
        envelope <- paste(sep = "\n",
                          paste("Date:", date()),
                          paste("Name:", input$cname),
                          paste("Email:", input$email),
                          paste("Comment:", input$suggestion))
        
        if(nchar(envelope) > 5000){
            updateNotifications(paste0("Too many characters."), "exclamation-triangle", "danger")
            
        } else {
            fwrite(list(envelope), paste0("contact/", format(as.POSIXct(Sys.time()), "%F.%Hh%Mm%Ss"), ".txt"))
            
            updateNotifications(paste0("Contact form submitted. Thanks!"), "thumbs-up", "success") 
        }
        
    })
    
    # Generate report ----

    output$downloadDemo <- downloadHandler(
        filename <- function() {
            paste0("Sample dataset ", input$selectDemoData, ".csv")
        },
        
        content <- function(file) {
            file.copy("data/donors.uniprot.csv", file)
        }
    )
    
    
    
    output$downloadReport <- downloadHandler(
        
        filename = function() {
            paste("proteomill-report.html")
        },
        content = function(file) {
            
            cat(file=stderr(), paste0(getwd()))
            
            src <- normalizePath('reports/report.rmd')
            
            file.copy('reports/report.rmd', 'reports/tmp/report.rmd', overwrite = TRUE) 
            
            out <- rmarkdown::render('reports/report.rmd')
            file.rename(out, file)
        }
    )
    
    # Data import wizard ----
    
    observeEvent(input$ShowHide1, {
        shinyjs::toggle("showHideBox1")
    })
    
    observeEvent(input$ShowHide2, {
        shinyjs::toggle("showHideBox2")
    })
    
    # File input: Main data ----
    
    observeEvent(input$file1, {
        
        inFile <- input$file1
        
        if (is.null(inFile)) return(NULL)
        
        dw <- tryCatch({fread(input=inFile$datapath);},
                       error = function(err){return(err)} ,
                       warning = function(war){return(war)} ,silent=F)
        
        if("error" %in% class(dw)) {
            
            createAlert(session, "fileHasError", "exampleAlert",
                        content = dw$message, append = FALSE,
                        style = "danger")
            
            rm(dw)
            
        } else if ("warning" %in% class(dw)) {
            
            createAlert(session, "fileHasWarning", "exampleAlert",
                        content = dw$message, append = FALSE,
                        style = "warning")
            
        } else {
            
            names(dw) <- make.names(names = names(dw), unique = T)
            
            maindata$data_wide_tmp <- dw
            
            shinyjs::show("previewDTInfo")
            shinyjs::show("dataDetailsWrapper")
            
            maindata$isValid <- validate_data_format(dw)
            
        }
        
    })
    
    
    output$previewDT <- renderDT({
        
        if(!is.null(maindata$data_wide_tmp)){
            
            dt <- maindata$data_wide_tmp
            
            dw <- tryCatch({
                
                if(ncol(dt) >= 20) {
                    datatable(dt[1:10, 1:20], options = list(autoWidth = T, scrollX = T))
                } else if (ncol(dt) < 20) {
                    datatable(dt[1:10, 1:ncol(dt)], options = list(autoWidth = T, scrollX = T))
                }
                
                ;},
                error = function(err){return(err)},
                warning = function(war){return(war)}, silent = F)
            
            if("error" %in% class(dw)) NULL
            else dw
            
        }
        
        
        
    })
    
    observeEvent(input$organism, {
        
        if(input$organism == "Other") {
            createAlert(session, "unsupportedOrganism", "exampleAlert",
                        content = "Please note that proceeding with organism 'other' means Pathway and Network analysis is unavailable. We are continuously adding support for more species.", append = FALSE,
                        style = "info")
        }

        
    })
    
    
    observeEvent(input$EndStep1, {
        
        shinyjs::show("Modal1Spinner")
        
        if(input$organism == "Homo sapiens | HSA | 9606") {
            
            library(EnsDb.Hsapiens.v86)
            
            maindata$organism <- EnsDb.Hsapiens.v86
            maindata$taxid <- 9606
            
        } else if(input$organism == "Mus musculus | MMU | 10090") {
            
            library(EnsDb.Mmusculus.v79)
            
            maindata$organism <- EnsDb.Mmusculus.v79
            maindata$taxid <- 10090
            
        } else if(input$organism == "Rattus norvegicus | RNO | 10116") {
            
            library(EnsDb.Rnorvegicus.v79)
            
            maindata$organism <- EnsDb.Rnorvegicus.v79
            maindata$taxid <- 10116
            
        } else if(input$organism == "Other") {
            
            maindata$organism <- NULL
            maindata$taxid <- NULL
        }
        
        
        toggleModal(session, "ImportModal1", toggle = "toggle")
        shinyjs::hide("Modal1Spinner")
        toggleModal(session, "ImportModal2", toggle = "toggle")
        
        
    })
    
    
    observeEvent(input$EndStep2, {

        if(is.null(maindata$data_wide_tmp)) {
            createAlert(session, "selectAFile", "exampleAlert",
                        content = "Please select a file.", append = FALSE,
                        style = "danger")
        }
        else {
            
            if(maindata$isValid) {
                
                shinyjs::show("Modal2Spinner")
                
                all_missing <- upload_data(i = "auto")
                
                if(all_missing) {
                    createAlert(session, "noneConverted", "exampleAlert",
                                content = "Couldn't convert any IDs. Did you choose the correct organism?", append = FALSE,
                                style = "warning")
                }

                if(uniqueN(sampleinfo$samples$treatment) > 8) {
                    app_meta$palette <- colorRampPalette(brewer.pal(8, "Accent"))(uniqueN(sampleinfo$samples$treatment))
                } else if (uniqueN(sampleinfo$samples$treatment) <= 8) {
                    app_meta$palette <- brewer.pal(uniqueN(sampleinfo$samples$treatment), "Accent")
                }
                
                toggleModal(session, "ImportModal2", toggle = "toggle")
                
                if(input$organism != "Other") subset_interactions()
                
                shinyjs::hide("Modal2Spinner")
                
                updateNotifications("Dataset uploaded.","check-circle", "success")
                updateTasks(text = "Upload a dataset", value = 100, color = "green", i = 0001)
                
            } else {
                
                createAlert(session, "checkRequirements", "exampleAlert1",
                            content = "File format is not valid. Please make sure dataset fulfills checklist requirements.", append = F,
                            style = "info")
            }
        }
    })

    
    # Heatmap Updated ----
    
    prepare_heatmap <- reactive({
        
        dw_cor <- maindata$udat@main
        
        if(!is.null(dw_cor)){
            dw_cor[is.na(dw_cor)] <- 0 # Impute
            dw_cor <- round(cor(dw_cor, method = input$corMethod), 2)
        }
    })
    
    output$heatmap <- renderPlotly({
        
        dw_cor <- prepare_heatmap()
        
        if(input$showGrid) gg = 0.75
        else gg = 0
        
        heatmaply_cor(dw_cor,
                      limits = c(min(dw_cor), 1),
                      col = viridis(n = 100),
                      fontsize_row = 8,
                      fontsize_col = 8,
                      grid_gap = gg,
                      plot_method = "plotly",
                      column_text_angle = 270,
                      show_dendrogram = c(T, F),
                      row_dend_left = T
        )
        
    })
    
    # PCA Updated ----
    
    observeEvent(input$multilevelPCA, {
        
        # Multilevel analysis can only be performed when at least one sample is repeated.
        sinfo <- sampleinfo$samples
        sinfo <- as.data.table(sinfo)
        
        if(length(Reduce(setdiff, sinfo[, .(list(unique(sinfo$replicate))), sinfo$treatment]$V1)) > 0){
            
            updateCheckboxInput(session = session, inputId = "multilevelPCA", value = F)
            
            updateNotifications("Multilevel analysis can only be performed when at least one sample is repeated.","exclamation-triangle", "danger")
        }
        
        
    }, ignoreInit = T)
    
    prepare_pca <- reactive({
        pca.data <- maindata$udat@main
        
        if(!is.null(pca.data)) {

            pca.data[is.na(pca.data)] <- 0 # "Impute" missing values as 0
            pca.data <- IMIFA::pareto_scale(pca.data)
            
            if(input$multilevelPCA == T) {
                
                ndim <- input$pcaDims[2] - input$pcaDims[1] + 1
                p.pca = mixOmics::pca(X = t(pca.data), ncomp = NULL, multilevel = sampleinfo$samples$replicate, logratio = 'none', scale = F, center = T)
                
                cat(file=stderr(), ndim)
                cat(file=stderr(), p.pca)
                
            } else {
                
                p.pca <- prcomp(t(pca.data), center = TRUE, scale. = F)
            }
            p.pca
        }
    })
    
    pca2d <- reactive({
        
        p.pca <- prepare_pca()
        
        si_treatment <- sampleinfo$samples$treatment
        si_replicate <- sampleinfo$samples$replicate

        mydf <- data.frame(pc1 = p.pca$x[, input$pcaDims[1]],
                           pc2 = p.pca$x[, input$pcaDims[2]],
                           si_treatment)
        
        poly.df <- mydf %>% 
            group_by(si_treatment) %>%
            do(.[chull(.$pc1, .$pc2),]) 
        
        if(input$showPolygons) {
            suppressWarnings(
                ggplot(mydf, aes(pc1, pc2, colour = as.factor(si_treatment))) +
                    geom_polygon(data = poly.df, fill = "grey", alpha = .15) +
                    geom_point(size = 5, aes(text = rownames(mydf))) +
                    xlab(paste0("PC", input$pcaDims[1])) +
                    ylab(paste0("PC", input$pcaDims[2])) +
                    scale_color_manual(values = app_meta$palette) +
                    theme_light() + theme(legend.title = element_blank())) -> p
        } else {
            suppressWarnings(
                ggplot(mydf, aes(pc1, pc2, colour = as.factor(si_treatment))) +
                    geom_point(size = 5, aes(text = rownames(mydf))) +
                    xlab(paste0("PC", input$pcaDims[1])) +
                    ylab(paste0("PC", input$pcaDims[2])) +
                    scale_color_manual(values = app_meta$palette) +
                    theme_light() + theme(legend.title = element_blank())) -> p
        }
    })
    
    output$PCAplots <- renderPlotly({
        
        if(input$pcaDims[2] - input$pcaDims[1] == 1) {
            pcaplot <- pca2d()
            ggplotly(pcaplot, tooltip = "text") %>% layout(legend = list(title=list(text='<b> Treatment </b>')))
        } else if (input$pcaDims[2] - input$pcaDims[1] == 2) {
            p.pca <- prepare_pca()
            plotly::plot_ly(x = p.pca$x[,input$pcaDims[1]],
                                 y = p.pca$x[,(input$pcaDims[2] - 1)],
                                 z = p.pca$x[,input$pcaDims[2]],
                                 text = rownames(p.pca$x),
                                 hoverinfo = "text",
                                 color = sampleinfo$samples$treatment,
                                 colors = brewer.pal(length(unique(sampleinfo$samples$treatment)), "Accent")
            ) %>%
                plotly::add_markers(marker=list(size = 15, line=list(width = 1, color = "black"))) %>%
                plotly::layout(scene = list(xaxis = list(title = paste0("PC", input$pcaDims[1])),
                                            yaxis = list(title = paste0("PC", (input$pcaDims[2] - 1))),
                                            zaxis = list(title = paste0("PC", input$pcaDims[2]))),
                               legend = list(title=list(text='<b> Treatment </b>')))
        }
        
        
    })
    

    output$scree <- renderPlot({
        
        p.pca <- prepare_pca()
        
        var_explained_df <- data.frame(PC = paste0("PC",1:length(sampleinfo$samples$samples)),
                                       var_explained=(p.pca$sdev)^2/sum((p.pca$sdev)^2))
        
        var_explained_df$PC <- factor(var_explained_df$PC, levels = var_explained_df$PC)
        
        var_explained_df %>%
            ggplot(aes(x = PC, y = var_explained))+
            geom_col() +
            labs(title="Scree plot: PCA on scaled data") +
            ggthemes::theme_clean() +
            theme(axis.text.x = element_text(vjust = 0.6), legend.position = "none") +
            scale_color_grey()
        
    })
    
    
    observeEvent(input$pcaDims, {
        
        if( (input$pcaDims[2] - input$pcaDims[1] > 2) | (input$pcaDims[2] - input$pcaDims[1] < 1)) {
            updateSliderInput(session = session, inputId = "pcaDims", value = c(input$pcaDims[1], (input$pcaDims[1] + 2)))
        }
        
    })
    
    # Differential expression volcano plot ----
    
    diffv <- reactive({
        
        if(c("t", "B") %in% colnames(maindata$udat@deoutput)) contrast <- maindata$udat@deoutput[, -c("t", "B")]
        else contrast <- maindata$udat@deoutput
        
        sID <- sampleinfo$sID
        contrast <- cbind(maindata$udat@identifiers[, ..sID], contrast)
        
        if(maindata$udat@deoutput[, .N] > 0) {
            
            diff_df <- contrast[, c(1, 2, 5)]
            colnames(diff_df)[2:3] <- c("Fold-change", "FDR")
            
            diff_df$group <- "Not Sign."
            
            diff_df[FDR < 0.05 & abs(`Fold-change`) < 1.5, "group"] <- "Sign."
            
            diff_df[FDR >= 0.05 & abs(`Fold-change`) > 1.5, "group"] <- "FC > 1.5"
            
            diff_df[FDR < 0.05 & abs(`Fold-change`) > 1.5, "group"] <- "Sign.  & FC > 1.5"
            
            diff_df$group <- as.factor(diff_df$group)
            
            diff_df$group <- factor(diff_df$group, levels = c("Sign.  & FC > 1.5", "Sign.", "FC > 1.5", "Not Sign."))
            
            as.data.frame(diff_df)
        }
    })
    
    output$dea_volcano <- renderPlotly({
        
        if(maindata$udat@deoutput[, .N] > 0) {
            diff_df <- diffv()
            
            nb.cols <- uniqueN(diff_df[, "group"])
            mycolors <- c(colorRampPalette(brewer.pal(8, "Accent"))(nb.cols))
            
            p <- plot_ly(data = diff_df,
                         type = "scatter",
                         x = ~`Fold-change`,
                         y = ~(-log10(FDR)),
                         text = ~get(sampleinfo$sID),
                         mode = "markers",
                         color = ~group,
                         colors = mycolors) %>%
                layout(yaxis = list(type = "log1p", showgrid = T,  ticks = "outside", autorange = T))
            p
        }

    })

    
    # Data validation ----
    
    validate_data_format <- function(dt) {
        
        cnames <- names(dt)[2:ncol(dt)]
        
        # Throw error if any false
        
        # Columns check
        test_col <- ncol(dt) >= 5
        test_row <- dt[, .N] >= 5
        test_uniq_cols <- all(duplicated(cnames) == F)
        test_col_sep <- all(grepl('_', cnames))
        test_starts_with <- all(grepl("^([[:alnum:]])", cnames))
        test_ends_with <- all(grepl("([[:alnum:]]$)", cnames))
        
        treatment <- as.factor(gsub('_.*', '', cnames))
        test_mult_treat <- length(unique(treatment)) >= 2
        test_mult_reps <- all(summary(treatment) >= 2 )
        
        test_dupl_acc <- uniqueN(dt[, 1]) == dt[, .N]
        
        test_res <- data.table(
            Test = c(
                "Number of columns",
                "Number of rows",
                "Duplicate columns",
                "Has separator",
                "Starts alpha-numeric",
                "Ends alpha-numeric",
                "Multiple treatments exist",
                "Multiple replicates exist",
                "Duplicate accessions"
            ),
            ErrMsg = c(
                "The dataset has too few (<5) or too many (>150) columns for ProteoMill's capacity.",
                "The dataset has too few (<5) or too many (>50.000) rows for ProteoMill's capacity.",
                "Duplicate column names is not allowed.",
                "Missing '_'-separator between treatment and replicate in column names.",
                "First character in column names must be alpha-numeric",
                "Last character in column names must be alpha-numeric",
                "Multiple treatments/groups required for differential analysis.",
                "Multiple replicates needed to estimate variance.",
                "The dataset has duplicate identifiers."
                
            ),
            Passed = c(
                test_col,
                test_row,
                test_uniq_cols,
                test_col_sep,
                test_starts_with,
                test_ends_with,
                test_mult_treat,
                test_mult_reps,
                test_dupl_acc
            )
        )
        
        if(any(test_res$Passed == F)) {
            createAlert(session, "fileFormattingError", "exampleAlert",
                        content = test_res[Passed==F, ErrMsg][1], append = FALSE,
                        style = "danger")
            
            return(F)
        } else {
            return(T)
        }
        
    }
    
    
    
    
    
    
}




