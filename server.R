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
#require(mixOmics)

# Parsing and reshaping
require(stringr)
require(dplyr)
require(data.table)
require(XML)
require(rmarkdown)
require(R.utils)
require(knitr)



# Generic functions ----

dframe <- function(dt, r){
    dt <- dt[!is.na(get(r))]                    # Remove NA
    rn <- dt[, as.character(.SD[[r]])]          # Rownames to vector
    df <- as.data.frame(dt[,-..convertColumns]) # 
    rownames(df) <- rn
    
    return(df)
}

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
                    renderUI({
                        tags$iframe(
                            src = paste0("doc/News.html"),
                            width = "100%",
                            height = "600px",
                            frameborder = "0")})))
        }
        
        
        if(input$sidebarmenu == "interactions"){
            v <- pathways$v
            if(is.null("v")){
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
                     # menuSubItem("UMAP", tabName = "UMAP", href = NULL, newtab = TRUE,
                     #             icon = shiny::icon("angle-double-right"), selected = F),
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
                     # menuSubItem("Confidence intervals", tabName = "diffexpoutput", href = NULL, newtab = TRUE,
                     #             icon = shiny::icon("angle-double-right"), selected = F)
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
    
    
    
    
    # Build sample info ----
    
    sample_data <- function(data) {
        samples <- names(data[, -..convertColumns])
        treatment <- as.factor(gsub('_.*', '', samples))
        replicate <- as.factor(gsub('.*_', '', samples))
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
        
        proteins <- maindata$data_wide[, "UNIPROTID"][[1]]
        
        interactions <- data.table::fread("lib/interactions5.txt.gz")
        
        pathways$ints <- interactions[(protein1 %in% proteins) & (protein2 %in% proteins)]
        
        # assign("proteins", proteins, envir = .GlobalEnv)
        # assign("ints", pathways$ints, envir = .GlobalEnv)
        
        
    }
    
    
    # File input: upload data function ----
    
    undup <- function(genes){
        genes[!is.na(genes)]
    }
    
    upload_data <- function(i = "auto"){
        
        data_wide <- maindata$data_wide
        
        # maindata$data_wide <- data.table::fread(
        #     path,
        #     sep = sep,
        #     dec = ".",
        #     header = T)
        
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
        
        # 10 first rows for converting IDs
        keys <- data_wide[, as.character(.SD[[1L]])][1:10]
        
        
        # Peptide data input
        
        # Guess input ID based on successful conversions
        if(i == "auto") {
            
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
            
            if(i == "UNIPROTID"){
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
            
            
            
            
        } else {
            
            if(i == "UNIPROTID"){
                
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
        }
        
        names(data_wide)[1] <- i
        
        setkeyv(tr_all, i)
        
        setkeyv(data_wide, i)
        
        tr_all <- tr_all[data_wide[, ..i], on = i]
        
        data_wide <- data_wide[tr_all, nomatch = 0]
        
        data_wide$ENTREZID <- as.character(data_wide$ENTREZID)
        
        for(j in seq_along(data_wide)){
            data.table::set(data_wide, i = which(duplicated(data_wide[[j]]) & is.character(data_wide[[j]])), j = j, value = NA)
        }
        
        data_wide[is.na(ENTREZID), ENTREZID := paste0("MISSING_", seq(1:length(is.na(ENTREZID))))]
        data_wide[is.na(SYMBOL), SYMBOL := paste0("MISSING_", seq(1:length(is.na(SYMBOL))))]
        data_wide[is.na(UNIPROTID), UNIPROTID := paste0("MISSING_", seq(1:length(is.na(UNIPROTID))))]
        data_wide[is.na(GENEID), GENEID := paste0("MISSING_", seq(1:length(is.na(GENEID))))]
        data_wide[is.na(PROTEINID), PROTEINID := paste0("MISSING_", seq(1:length(is.na(PROTEINID))))]
        
        setcolorder(data_wide, c(convertColumns, names(maindata$data_origin[,2:ncol(maindata$data_origin)])))
        
        maindata$data_wide <- data_wide
        maindata$data_origin <- data_wide
        
        pdesc <- data.table::fread("lib/protein_descriptions.txt.gz")

        setkey(maindata$data_wide, "UNIPROTID")
        setkey(pdesc, "UNIPROTID")
        
        pdesc <- pdesc[maindata$data_wide[,1:5]]
        maindata$pdesc <- pdesc
        
        sample_data(maindata$data_wide)
        
    }
    
    # File input: Main data ----
    # observeEvent(input$infile, {
    #     
    #     maindata$inFile <- input$infile
    #     
    #     if (is.null(maindata$inFile))
    #         return(NULL)
    # 
    #     upload_data(maindata$inFile$datapath, separator(input$dataSep), identifier(input$dataIdentiferType))
    # 
    #     updateNotifications("Dataset uploaded.","check-circle", "success")
    #     updateTasks(text = "Upload a dataset", value = 100, color = "green", i = 0001)
    #     
    # })
    
    # File input: Annotation ----
    observeEvent(input$anno_infile, {
        maindata$inFile <- input$anno_infile
        
        if (is.null(maindata$inFile))
            return(NULL)
        
        maindata$data_annotation <- data.table::fread(
            maindata$inFile$datapath,
            sep = separator(input$annoSep),
            dec = ".",
            header = T)
        

        updateNotifications("Annotations uploaded.","check-circle", "success")
        updateTasks(text = "Upload annotation data", value = 100, color = "green", i = 0002)
        
        
    })
    
    # File input: Demo data ----
    observeEvent(input$useDemoData, {

        # Sample 1
        
        sample_1_exp <- "data/donors.uniprot.csv"
        maindata$data_wide <- fread(sample_1_exp, header = T)
        #sample_1_anno <- "data/donors.uniprot.annotation.tsv"
        
        updateNotifications(paste0("Loading human annotation library..."), "info-circle", "info")
        library(EnsDb.Hsapiens.v86)
        
        updateNotifications(paste0("Loading human interaction library..."), "info-circle", "info")
        upload_data()
        subset_interactions()
        
        
        # maindata$data_annotation <- data.table::fread(
        #     sample_1_anno,
        #     sep = "auto",
        #     dec = ".",
        #     header = T)
        
        updateNotifications("Demo data uploaded.","check-circle", "success")
        updateTasks(text = "Upload a dataset", value = 100, color = "green", i = 0001)
        updateTasks(text = "Upload annotation data", value = 100, color = "green", i = 0002)
        
    })
    
    
    
    # Data summary plot ----
    
    output$datainfoBox <- renderInfoBox({
        if(!is.null(sampleinfo$samples)) infoBox("Proteins", maindata$data_wide[, .N])
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
                maindata$data_wide <- maindata$data_wide[, c(convertColumns, as.character(input$includesamples)), with = F]
                maindata$data_wide <- maindata$data_wide[apply(maindata$data_wide[, -convertColumns, with = F], 1, FUN = function(x) !all(is.na(x))),]
                
                maindata$data_origin <- maindata$data_origin[, c(convertColumns, as.character(input$includesamples)), with = F]
                maindata$data_origin <- maindata$data_origin[apply(maindata$data_origin[, -convertColumns, with = F], 1, FUN = function(x) !all(is.na(x))),]
                
                sample_data(maindata$data_wide)
                
            } else {
                updateNotifications("At least two groups needed for comparison.","exclamation-triangle", "danger")
            }
        } else {
            updateNotifications("At least two groups needed for comparison.","exclamation-triangle", "danger")
        }
        

        
    })
    
    
    output$identifierinfo <- renderTable({
        if(!is.null(maindata$data_wide)) {
            
            missingids <- maindata$data_wide[, ..convertColumns]
            missingidspc <- apply(missingids, 2, FUN = function(x) paste0(round(((sum(!startsWith(x, "MISSING")) / missingids[, .N]) * 100), 2), "% (", sum(!startsWith(x, "MISSING")), ")") )
            as.data.frame(missingidspc)
        }
    }, rownames = T, colnames = F)
    
    
    output$violinplot <- renderPlot({

        if(!is.null(maindata$data_wide) & !is.null(sampleinfo$samples)){
            dat <- stack(as.data.frame(maindata$data_wide[, 6:ncol(maindata$data_wide)]))
            dat <- dat[!is.na(dat$values), ]
            dat$sample <- sub("_.*", "", dat$ind)
            dat$ind <- factor(dat$ind, levels = levels(dat$ind)[order(levels(dat$ind))])
            dat$values <- log2(dat$values)
            
            if(input$violintype == 1) {
                
                plots$violin1 <- ggplot(dat, 
                       aes(x = sample, y = values, fill = sample)) + 
                    geom_violin(trim = FALSE, scale = "width") +
                    geom_boxplot(width=0.1, fill="white") +
                    ggthemes::theme_clean() +
                    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.position = "none") +
                    xlab("") + ylab("Log2 Expression") +
                    scale_fill_brewer(palette = "BuGn")
                
                plots$violin1
                
            } else {
                
                plots$violin2 <- ggplot(dat, 
                       aes(x = ind, y = values, fill = sample)) + 
                    geom_violin(trim = FALSE, scale = "width") +
                    geom_boxplot(width=0.1, fill="white") +
                    ggthemes::theme_clean() +
                    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.position = "none") +
                    xlab("") + ylab("Log2 Expression") +
                    scale_fill_brewer(palette = "BuGn")
                
                plots$violin2
                
            }
            


        }

            
        
        
    })
    
    # Missing values: frequency plot ----
    
    renderNAfreq <- function(data_wide){
        output$nafreq <- renderPlot({

            # Which elements are NA?
            allNA <- is.na(data_wide)
            
            # Summary of how many TRUEs there are in each row
            NA_frequency <- table(rowSums(allNA))
            
            naf <- as.data.frame(NA_frequency)
            
            
            # Draw plot
            plots$mv <- ggplot(naf, aes(x = Var1, y = Freq, color = Var1)) +
                geom_bar(stat = "identity", width = .9, fill = "white") +
                labs(x = "Number of missing values in at least one sample", y = "Number of proteins") +
                ggthemes::theme_clean() +
                theme(axis.text.x = element_text(vjust = 0.6), legend.position = "none") +
                scale_color_grey()
            
            plots$mv
            
        })
    }
    
    
    
    observeEvent(input$loadfilterplot, {
        
        data_wide <- maindata$data_wide
        
        if(!is.null(data_wide)){
            renderNAfreq(data_wide)
        } else {
            updateNotifications("Upload a dataset first.","exclamation-triangle", "danger")
        }
    })
    
    observeEvent(input$setcutoff, {
        
        
        
        data_wide <- maindata$data_wide
        
        subsample_data()
        
        if(!is.null(data_wide)){
            data_origin <- maindata$data_origin
            data_wide <- subset_by_na(dataset = data_origin, treatment = sampleinfo$samples$treatment, threshold = input$missingvalues)
            
            #subset_interactions()
            
            unlock_menus()
            renderNAfreq(data_wide)
            
            # assign("sic", sampleinfo$samples$treatment, envir = .GlobalEnv)
            # assign("sir", sampleinfo$samples$replicate, envir = .GlobalEnv)
            # assign("sis", sampleinfo$samples, envir = .GlobalEnv)
            # assign("siSID", sampleinfo$sID, envir = .GlobalEnv)
            # assign("sgroup", sampleinfo$group, envir = .GlobalEnv)
            
            
            maindata$data_wide <- data_wide
            updateTasks(text = "Set a filter", value = 100, color = "green", i = 0003)
            updateNotifications(paste0("NA cutoff set to ", input$missingvalues, ".") ,"check-circle", "success")
        } else {
            updateNotifications("Upload a dataset first.","exclamation-triangle", "danger")
        }
        
    })
    
    
    
    # Quality control ----
    
    # Plot PCA function ----
    
    plotPCA <- function(contribs, ellipse, type) {
        
        data_origin <- maindata$data_origin
        dt <- dframe(subset_by_na(data_origin, treatment = sampleinfo$samples$treatment, threshold = input$missingvalues), sampleinfo$sID)
        # Biplot extension displaying top contributing proteins currently only available for 2D plot.
        
        if(type == '2d') {
            
            pca.data <- log2(dt)  # Log2 transform data
            pca.data[is.na(pca.data)] <- 0 # "Impute" missing values as 0
            
            pca.data <- IMIFA::pareto_scale(pca.data)
            
            #pca.data <- t(pca.data)        # Transpose dataset
            
            p.pca <- prcomp(t(pca.data), center = TRUE, scale. = F)
            
            plots$pcaplot2d <- factoextra::fviz_pca_biplot(p.pca, title = '', label = "var", habillage = sampleinfo$samples$treatment,
                                                   addEllipses = TRUE, ellipse.level = ellipse,
                                                   select.var = list(contrib = contribs), repel = TRUE) + theme_light()
            
            return (plots$pcaplot2d)
            
        } else if (type == '2dpaired') {
            
            
            pca.data <- dt
            pca.data[is.na(pca.data)] <- .5 # "Impute" missing values as 0

            pca.data <- pca.data[rowSums(pca.data) != .5,]
            
            pca.data <- t(pca.data)        # Transpose dataset
            
            pca.res = pca(X = pca.data,
                          multilevel = sampleinfo$samples$replicate, logratio = 'CLR')
            
            plots$pcaplot2d <- plotIndiv(
                pca.res,
                ind.names = sampleinfo$samples$replicate,
                group = sampleinfo$samples$treatment,
                title = "Multilevel PCA",
                legend = T,
                style = "ggplot2",
                ellipse = FALSE,
                ellipse.level = .9
            )
            return(plots$pcaplot2d)
            
            
            
        } else if (type == '3d') {
            
            pca.data <- log2(dt)  # Log2 transform data
            pca.data[is.na(pca.data)] <- 0 # "Impute" missing values as 0
            
            pca.data <- IMIFA::pareto_scale(pca.data)
            
            #pca.data <- t(pca.data)        # Transpose dataset
            
            p.pca <- prcomp(t(pca.data), center = TRUE, scale. = F)
            
            
            plots$pcaplot3d <- plotly::plot_ly(x = p.pca$x[,1],
                                       y = p.pca$x[,2],
                                       z = p.pca$x[,3],
                                       text = rownames(p.pca$x),
                                       hoverinfo = "text",
                                       color = sampleinfo$samples$treatment,
                                       colors = c("red","green","blue"),
                                       sizes = c(100, 150)) %>%
                plotly::add_markers() %>%
                plotly::layout(scene = list(xaxis = list(title = 'PC1'),
                                            yaxis = list(title = 'PC2'),
                                            zaxis = list(title = 'PC3')))
            
            return (plots$pcaplot3d)
        } else if (type == '3dpaired') {
            
            pca.data <- dt
            pca.data[is.na(pca.data)] <- .5 # "Impute" missing values as 0.5
            
            pca.data <- pca.data[rowSums(pca.data) != .5,]
            
            pca.data <- t(pca.data)        # Transpose dataset
            
            pca.res = pca(X = pca.data,
                          ncomp = 3,
                          multilevel = sampleinfo$samples$replicate, logratio = 'CLR')
            
            plots$pcaplot3d <- plotly::plot_ly(x = pca.res$x[,1],
                                       y = pca.res$x[,2],
                                       z = pca.res$x[,3],
                                       text = rownames(pca.res$x),
                                       hoverinfo = "text",
                                       color = sampleinfo$samples$treatment,
                                       colors = c("red","green","blue"),
                                       sizes = c(100, 150)) %>%
                plotly::add_markers() %>%
                plotly::layout(scene = list(xaxis = list(title = 'PC1'),
                                            yaxis = list(title = 'PC2'),
                                            zaxis = list(title = 'PC3')))

            return(plots$pcaplot3d)
            
        }  else { return (FALSE) }
        
    }
    
    # Render PCA plots
    
    observeEvent(input$loadPCAplots, {
        
        # If unpaired
        if(input$diffexppairing == 2){
            output$pca2dplot <- renderPlot({
                plotPCA(input$contribs, input$ellipse, '2d')
            })
        } else if(input$diffexppairing == 1) {
            output$pca2dplot <- renderPlot({
                plotPCA(0, input$ellipse, '2dpaired')

            })
        }

        # If unpaired
        if(input$diffexppairing == 2){
            output$pca3dplot <- plotly::renderPlotly({
                plotPCA(0,0, '3d')
            })
        } else if(input$diffexppairing == 1) {
            output$pca3dplot <- plotly::renderPlotly({
                plotPCA(0, 0, '3dpaired')
                
            })
        }

        updateTasks(text = "Inspect data", value = (tasks$tasks[id == "0004", value] + 100/3), color = "green", i = 0004)
    })
    

    
    # Render UMAP
    observeEvent(input$loadUMAP, {
        output$UMAPplot <- renderPlot({
            plotPCA(0,0, 'UMAP')
        })
        updateTasks(text = "Inspect data", value = (tasks$tasks[id == "0004", value] + 100/3), color = "green", i = 0004)
    })
    

    # Render heatmap
    
    observeEvent(input$generateheatmap, {
        cor_mat_raw_logged <- log2(maindata$data_origin[,-..convertColumns])
        cor_mat_raw_logged[is.na(cor_mat_raw_logged)] <- 0
        cor_mat_raw_logged <- cor(cor_mat_raw_logged)
        
        #annotation_col = data.frame(tmp)
        
        data_annotation <- maindata$data_annotation
        
        if(!is.null(data_annotation)){
            if(all(data_annotation[, as.character(.SD[[1L]])] %in% names(maindata$data_origin[,-..convertColumns]))){
                plots$hmap <- pheatmap::pheatmap(
                    annotation_col = dframe(data_annotation, "V1")[2:4],
                    cor_mat_raw_logged,
                    legend_breaks = c(min(cor_mat_raw_logged), 1),
                    legend_labels = c(0, 1),
                    silent = T
                )
            }
        } else {
            plots$hmap <- pheatmap::pheatmap(
                cor_mat_raw_logged,
                legend_breaks = c(min(cor_mat_raw_logged), 1),
                legend_labels = c(0, 1),
                silent = T
            )
        }
        
        output$samplecorrheatmap = renderPlot({
            plots$hmap
            })
        
        updateTasks(text = "Inspect data", value = (tasks$tasks[id == "0004", value] + 100/3), color = "green", i = 0004)
        
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
        output$enrichment <- renderMenu({
            menuItem("Enrichment analysis", icon = icon("flask"), href = NULL,
                     menuSubItem("Pathway enrichment", tabName = "pathwayenrichment", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F),
                     menuSubItem("Pathway visualization", tabName = "pathwayvisualization", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F)
            )
        })
        output$network <- renderMenu({
            menuItem("Network analysis", icon = icon("vector-square"),
                     menuSubItem("Interactions", tabName = "interactions", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F))
        })
        
        pairing <- input$diffexppairing
        
        
        
        contrasts <- paste(input$contrast1,input$contrast2, sep = '-')

        contrast <- diff_exp(contrasts, pairing)

        rcont$contrasts <- contrasts
        rcont$contrast <- contrast
        
        
        removeUI(selector = "#enrichrm")
        removeUI(selector = "#networkrm")
        removeUI(selector = "#interactionsrm")
        updateNotifications("A DE contrast has been set.","check-circle", "success")
        updateTasks(text = "Set contrast", value = 100, color = "green", i = 0005)
    })
    
    # Differential expression: DEA function ----
    
    diff_exp <- function(coeff, pairing) {
        
        dw <- maindata$data_wide
        sinf <- sampleinfo$samples

        if(input$setDEengine == 2) {
            
            dw <- maindata$data_wide
            
            # Replace NA -> 0
            for(j in seq_along(dw)){
                set(dw, i = which(is.na(dw[[j]]) & is.numeric(dw[[j]])), j = j, value = 0)
            }
            
            if(pairing == 1) {
                # Paired
                
                design <- DESeqDataSetFromMatrix(countData  = round(dframe(dw, sampleinfo$sID)),
                                                 colData    = sinf,
                                                 design     = ~ 0 + sampleinfo$samples$treatment + sampleinfo$samples$replicate)
                
            } else {
                # Unpaired
                design <- DESeqDataSetFromMatrix(countData  = round(dframe(dw, sampleinfo$sID)),
                                                   colData    = sinf,
                                                   design     = ~ 0 + sampleinfo$samples$treatment)
            }

            dds <- DESeq(design)
            
            contrast <- results(dds, contrast=c("treatment", sub("treatment", "", input$contrast1), sub("treatment", "", input$contrast2)))
            
            colnames(contrast) <- c("baseMean", "logFC", "logFC.SE", "stat", "P.Value", "adj.P.Val")
            
            # Confidence intervals used for plot, global var
            cint <- contrast
            cint$protein <- rownames(cint)
            cint$protein <- factor(cint$protein, levels = cint$protein[order(cint$logFC)])
            
            rcont$cint <- cint
            
            contrast <- contrast[order(contrast$P.Value, decreasing = F),]
            contrast <- data.table::as.data.table(contrast, keep.rownames = T)
            
            setkeyv(contrast, "rn")
            setkeyv(maindata$data_wide, sampleinfo$sID)
            contrast <- contrast[maindata$data_wide[,..convertColumns], nomatch = 0]
            names(contrast) <- c(sampleinfo$sID, names(contrast[, 2:ncol(contrast)]))
            setcolorder(contrast, c(convertColumns, "logFC", "logFC.SE", "baseMean", "stat", "P.Value", "adj.P.Val"))
            
            # Reduce network load
            if(contrast[, .N] > 2000) {
                # x <- unname(quantile(abs(contrast$logFC)))
                # x <- round(x[5] - x[4])
                updateNumericInput(session = session, "fccutoff", value = 2)
            }

        } else if(input$setDEengine == 1) {
            
            
            # Create Annotation data and expression set (Biobase)
            phenoData <- new("AnnotatedDataFrame", data = sampleinfo$samples)
            exampleSet <- ExpressionSet(assayData = as.matrix(log2(dframe(maindata$data_wide, sampleinfo$sID))), phenoData = phenoData)
            
            treatment <- sampleinfo$samples$treatment
            repl <- sampleinfo$samples$replicate
            
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
            
            # Confidence intervals used for plot, global var
            cint <- contrast
            cint$protein <- rownames(cint)
            cint$protein <- factor(cint$protein, levels = cint$protein[order(cint$logFC)])
            
            rcont$cint <- cint
            
            contrast <- contrast[order(contrast$P.Value, decreasing = F),]
            contrast <- data.table::as.data.table(contrast, keep.rownames = T)
            
            setkeyv(contrast, "rn")
            setkeyv(maindata$data_wide, sampleinfo$sID)
            contrast <- contrast[maindata$data_wide[,..convertColumns], nomatch = 0]
            names(contrast) <- c(sampleinfo$sID, names(contrast[, 2:ncol(contrast)]))
            setcolorder(contrast, c(convertColumns, "logFC", "t", "P.Value", "adj.P.Val", "B"))
            
            # Reduce network load
            if(contrast[, .N] > 2000) {
                # x <- unname(quantile(abs(contrast$logFC)))
                # x <- round(x[5] - x[4])
                updateNumericInput(session = session, "fccutoff", value = 2)
            }
            
        }
        
        #assign("contrast", contrast, envir = .GlobalEnv)
        
        return( contrast )
        
    }
    
    output$diffexptable_summary <- renderTable({
        contrast <- rcont$contrast
        
        if(!is.null(contrast)){
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
        
        contrast <- rcont$contrast
        
        if(!is.null(contrast)){
            contrast <- contrast[abs(logFC) >= input$diffexp_limit_fc]
            contrast <- contrast[adj.P.Val < input$diffexp_limit_pval]
            contrast <- contrast[order(logFC, decreasing = T)]
            maindata$diffexptable_up <- DT::datatable(dframe(contrast, sampleinfo$sID),
                                options = list(order = list(1, 'desc'))) %>% 
                formatRound(columns=1:(ncol(contrast) - length(convertColumns)), digits=4)
            
            maindata$diffexptable_up
            
            
        }

    })
    
    output$diffexptable_down <- DT::renderDataTable({
        
        contrast <- rcont$contrast
        
        if(!is.null(contrast)){
            contrast <- contrast[abs(logFC) >= input$diffexp_limit_fc]
            contrast <- contrast[adj.P.Val < input$diffexp_limit_pval]
            contrast <- contrast[order(logFC, decreasing = F)]
            df <- DT::datatable(dframe(contrast, sampleinfo$sID),
                                options = list(order = list(1, 'asc'))) %>% 
                formatRound(columns=1:(ncol(contrast) - length(convertColumns)), digits=4)
            # df %>% DT::formatSignif('logFC', digits = 2)
            # df %>% DT::formatSignif('CI.L', digits = 2)
            # df %>% DT::formatSignif('CI.R', digits = 2)
            # df %>% DT::formatSignif('adj.P.Val', digits = 2)
            
        }
        
    })
    
    # Differential expression: confidence intervals
    
    output$contrasttable <- renderPlot({
        
        if(!is.null(rcont$cint)){
            ggplot2::ggplot(rcont$cint, aes(x = protein, y = logFC, colour = logFC)) +
                coord_flip() +
                geom_errorbar(aes(ymin = as.numeric(CI.L), ymax = as.numeric(CI.R)), width = 2) +
                scale_color_viridis_c() +
                geom_line() +
                geom_point() +
                theme_bw() +
                theme(
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    axis.text.y = element_text(size = 7, family = "Helvetica")
                )
        }

    })
    
    
    # Pathways ----
    
    observeEvent(input$min_fc, {
      
      output$number_of_genes <- renderUI({
        
        contrast <- rcont$contrast
          
        up <- paste("Up-regulated:", contrast[logFC >= input$min_fc & adj.P.Val < input$min_pval, .N], "proteins", sep = " ")
        down <- paste("Down-regulated:", contrast[logFC < (input$min_fc * -1) & adj.P.Val < input$min_pval, .N], "proteins", sep = " ")
        
        HTML(paste(up, down, '<br/>', sep = '<br/>'))
        
      })
    })
    
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
        g <- paste(unlist(up[,genes]), collapse = ", ")
        
        show_selected(p, g)
        
    })
    
    observeEvent(input$downregulated_pathways_table_rows_selected, {
        DOWNREGULATED_pathways <- pathways$DOWNREGULATED_pathways
        i = input$downregulated_pathways_table_rows_selected[1]
        down <- DOWNREGULATED_pathways[i,]
        p <- paste("<b>", down[,Pathway_name], "</b>")
        g <- paste(unlist(down[,genes]), collapse = ", ")
        
        show_selected(p, g)
        
        
    })
    
    
    
    
    observeEvent(input$generate_pathways, {
        
        contrast <- rcont$contrast
        
        if(contrast[logFC >= input$min_fc & adj.P.Val < input$min_pval, .N] < 3 | contrast[logFC >= input$min_fc & adj.P.Val < input$min_pval, .N] < 3){
            updateNotifications("Too few differential proteins.","exclamation-triangle", "danger")
        } else {
            
            #updateNumericInput(session = session, inputId = "pvaluecutoff", label = "Maximum adj. Pvalue", value = input$min_pval)
            #updateNumericInput(session = session, inputId = "fccutoff", label = "Minimum abs. log2FC", value = input$min_fc)
            
            UPREGULATED_genes <- contrast[logFC >= input$min_fc & adj.P.Val < input$min_pval, UNIPROTID]
            DOWNREGULATED_genes <- contrast[logFC < (input$min_fc * -1) & adj.P.Val < input$min_pval, UNIPROTID]
            
            UPREGULATED_pathways <- ora(UPREGULATED_genes)@output
            pathways$UPREGULATED_pathways <- UPREGULATED_pathways
            
            DOWNREGULATED_pathways<- ora(DOWNREGULATED_genes)@output
            pathways$DOWNREGULATED_pathways <- DOWNREGULATED_pathways
            
            output$upregulated_pathways_table <- DT::renderDT(
                
                pathways$upregulated_pathways_table <- DT::datatable(UPREGULATED_pathways[, -c("genes", "background")],
                              selection = 'single',
                              options = list(autoWidth = TRUE,
                                             scrollX=TRUE,
                                             columnDefs = list(
                                                 list(width = '100px', targets = c(1, 3)),
                                                 list(width = '60px', targets = c(6, 7))
                                             )))
            )
            
            pathways$upregulated_pathways_table
            
            pathways$downregulated_pathways_table <- output$downregulated_pathways_table <- DT::renderDT(
                DT::datatable(DOWNREGULATED_pathways[, -c("genes", "background")],
                              selection = 'single',
                              options = list(autoWidth = TRUE,
                                             scrollX=TRUE,
                                             columnDefs = list(
                                                 list(width = '100px', targets = c(1, 3)),
                                                 list(width = '60px', targets = c(6, 7))
                                             )))
            )
            
            pathways$downregulated_pathways_table
            
            updateTasks(text = "Run pathway enrichment", value = 100, color = "green", i = 0006)
            updateNotifications("Pathway analysis complete.","check-circle", "success")
        }

        
      
    })
    
    # Volcano plot
    
    observeEvent(input$loadPathwayPlots, {
        
        UPREGULATED_pathways <- pathways$UPREGULATED_pathways
        DOWNREGULATED_pathways <- pathways$DOWNREGULATED_pathways

        if(!exists("UPREGULATED_pathways")){
            updateNotifications("Run pathway analysis first.","exclamation-triangle", "danger")
        } else {
            
            output$volcano_plot <- renderPlotly(
                {
                    contrast <- rcont$contrast
                    sID <- sampleinfo$sID
                    
                    res <- enrichment_results(UPREGULATED_pathways, DOWNREGULATED_pathways, contrast)
                    
                    
                    # # # #
                    
                    tba <- contrast[!(UNIPROTID %in% res$Gene), c("UNIPROTID", "logFC", "P.Value")]
                    colnames(tba) <- c("Gene", "logFC", "P.Value")
                    
                    tba <- tba[!is.na(logFC)]
                    
                    tba[, c("ReactomeID", "P.Adj")] <- NA
                    
                    tba[, "TopReactomeName"] <- "[No significant over-representation]"
                    tba[, "Pathway_name"] <- "[No significant over-representation]"
                    
                    setcolorder(tba, colnames(res))
                    
                    res <- rbindlist(list(res, tba))
                    
                    
                    # # # #
                    
                    setkeyv(res, "Gene")
                    setkeyv(contrast, "UNIPROTID")
                    res <- res[contrast[,..convertColumns], nomatch = 0]
                    
                    names(res) <- c("UNIPROTID", c("ReactomeID", "Pathway_name", "TopReactomeName", "P.Adj", "logFC", "P.Value"), convertColumns[-1])
                    setcolorder(res, c(convertColumns, c("ReactomeID", "Pathway_name", "TopReactomeName", "P.Adj", "logFC", "P.Value")))

                    res$Pathway_name <- stringr::str_wrap(res$Pathway_name, 50)

                    pathways$v <- volcano(res, "Lowest", sID)
                    

                    plots$volcano <- plotly::ggplotly(pathways$v$volcano_plot) %>%
                        layout(dragmode = "select")
                    
                    plots$volcano

                })

            
            output$volcano_plot2 <- renderPlotly(
                {
                    gg <- plotly::ggplotly(pathways$v$volcano_plot, width = 700, height = 500) %>%
                        layout(dragmode = "select") %>% 
                        config(scrollZoom = T)
                    gg$x$data[[1]]$visible <- 'legendonly'  
                    gg
                    
                })
            
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
            
            description <- maindata$pdesc[get(sampleinfo$sID) == input$hovered_node, annotation]
            
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
    
    # observeEvent(input$loadNetworkplots, {
    #     
    #     ints2 <- interactions[(protein1 %in% proteins) & (protein2 %in% proteins)]
    #     
    # })
    
    output$interaction_network <- renderVisNetwork({
        
        contrast <- rcont$contrast
        res <- pathways$v$res
        sID <- sampleinfo$sID
        
        dir <- input$network_regulation

        if(dir == 1){
            proteins <- contrast[(abs(logFC) >= input$fccutoff) &
                                     (logFC > 0)]$UNIPROTID
        } else if(dir == 2){
            proteins <- contrast[(abs(logFC) >= input$fccutoff) &
                                     (logFC <= 0)]$UNIPROTID
        } else {
            proteins <- contrast[(abs(logFC) >= input$fccutoff)]$UNIPROTID
        }
        
        
        ints2 <- pathways$ints[(protein1 %in% proteins) & (protein2 %in% proteins)]
        ints2 <- ints2[score > input$interactioncutoff]
        

        mynodes <- unlist(event_data("plotly_selected")$key)
        
        mynodes <- contrast[get(sID) %in% mynodes, "UNIPROTID"][[1]]

        
        # TBA: all interactions option
        
        interacts <- function(i){
            return(ints2[(protein1 == i & protein2 %in% mynodes) | (protein2 == i & protein1 %in% mynodes), .N])
        }
        
        if(input$interaction_behaviour == 1){
            
            if(!is.null(mynodes)) {
                if(sum(unlist(lapply(mynodes, interacts))) > 0){
                    ints2 <- ints2[protein1 %in% mynodes & protein2 %in% mynodes]
                }
            }
            
        } else if(input$interaction_behaviour == 2){
            
            if(any(mynodes %in% ints2$protein1 & mynodes %in% ints2$protein2)){
                ints2 <- ints2[protein1 %in% mynodes | protein2 %in% mynodes]
            }
            
        }
        
        nodes <- data.table(id = unique(c(ints2$protein1, ints2$protein2)))
        nodes$label <- nodes$id

        
        edges <-
            data.table(
                from = ints2$protein1,
                to = ints2$protein2
            )

        g <- igraph::graph_from_data_frame(edges, directed = F, vertices = nodes)
        
        r <- res[V(g)$label, on = "UNIPROTID", ]
 
        
        g <- set.vertex.attribute(g, name = "name",  value = as.vector (r[, ..sID][[1]]))
        g <- set.vertex.attribute(g, name = "color", value = r$col)
        g <- set.vertex.attribute(g, name = "group", value = r$TopReactomeName)
        

        layoutf <- function(type) {
            switch(type,
                   "1" = "layout_nicely",
                   "2" = "layout_in_circle",
                   "3" = "layout_on_grid",
                   "4" = "layout_on_sphere",
                   "5" = "layout_randomly",
                   "6" = "layout_DH")
        }

        if(exists("g")){
            V(g)$name <- res[V(g)$name, on = sID, ..sID][[1]]
        } else {
            names(V(g)) <- res[V(g)$name, on = "UNIPROTID", "UNIPROTID"][[1]]
        }
        
            
        visNetwork::visIgraph(g) %>%
            visInteraction(hover = TRUE) %>%
            visIgraphLayout(layout = layoutf(1), randomSeed = 1) %>%
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
        
        data_wide <- maindata$data_wide
        
        if(!is.null(data_wide)){
            
            data_wide_NAex <- na.exclude(dframe(maindata$data_origin[1:100], sampleinfo$sID))
            
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
        
        data_wide <- maindata$data_wide
        
        if(!is.null(data_wide)){
            
            updateNotifications(paste0("Checking for outdated IDs. Please wait."), "info-circle", "info")
            
            candidates <- data_wide[apply(
                data_wide[, ..convertColumns],
                1,
                FUN = function(x)
                    all(startsWith(x[2:5], "MISSING"))
            ), "UNIPROTID"]
            
            # EventTime <- Sys.time() + candidates[, .N] * 5
            # output$eventTimeRemaining <- renderText({
            #     invalidateLater(1000, session)
            #     paste("Estimated time remaining:", 
            #           round(difftime(EventTime, Sys.time(), units='secs')), 'secs')
            # })
            
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
    
    # getOverview <- reactive({
    #     data.table(
    #         "Variable" = c("Total proteins", "DE proteins (adj. P. < 0.05)"),
    #         "Value" = c(maindata$data_origin[, .N], contrast[adj.P.Val < 0.05, .N])
    #         
    #     )
    #     
    # })
    # 
    # 
    # getDT <- reactive({
    #     maindata$data_wide
    #     
    # })
    # 
    # 
    # getpca2d <- reactive({
    #     pca2d
    # })
    # 
    # 
    # getUpPathways <- reactive({
    #     UPREGPATH <- UPREGULATED_pathways[, -c("ReactomeID", "background")]
    #     data.table::setcolorder(UPREGPATH, c("Pathway_name", "TopReactomeName", "q", "m", "p", "p.adj", "genes"))
    #     UPREGPATH
    # })
    # 
    # getDownPathways <- reactive({
    #     DOWNREGPATH <- DOWNREGULATED_pathways[, -c("ReactomeID", "background")]
    #     data.table::setcolorder(DOWNREGPATH, c("Pathway_name", "TopReactomeName", "q", "m", "p", "p.adj", "genes"))
    #     DOWNREGPATH
    # })
    # 
    # getVolcano <- reactive({
    #     res <- enrichment_results(UPREGULATED_pathways, DOWNREGULATED_pathways)
    #     
    #     tba <- contrast[!(UNIPROTID %in% res$Gene), c("UNIPROTID", "logFC", "CI.L", "CI.R", "P.Value")]
    #     colnames(tba) <- c("Gene", "logFC", "CI.L", "CI.R", "P.Value")
    #     
    #     tba <- tba[!is.na(logFC)]
    #     
    #     tba[, c("ReactomeID", "P.Adj")] <- NA
    #     
    #     tba[, "TopReactomeName"] <- "[No significant over-representation]"
    #     tba[, "Pathway_name"] <- "[No significant over-representation]"
    #     
    #     setcolorder(tba, colnames(res))
    #     
    #     res <- rbindlist(list(res, tba))
    #     
    #     res$Pathway_name <- stringr::str_wrap(res$Pathway_name, 50)
    #     
    #     setkeyv(res, "Gene")
    #     setkeyv(contrast, "UNIPROTID")
    #     res <- res[contrast[,..convertColumns], nomatch = 0]
    #     names(res) <- c("UNIPROTID", names(tba[,2:ncol(tba)]), convertColumns[-1])
    #     
    #     setcolorder(res, c(convertColumns, names(tba[,2:ncol(tba)])))
    # 
    #     v <- volcano(res, abstraction = input$abstractionlevel)
    #     
    # })
    # 
    # getReg <- reactive({
    #     data.table("Upregulated" = contrast[logFC <= 0 & adj.P.Val < 0.05, .N], "Downregulated" = contrast[logFC > 0 & adj.P.Val < 0.05, .N])
    #     
    # })


    

    
    
    output$download <- downloadHandler(
        filename = function(){
            if(!is.null(maindata$inFile$name)) paste(gsub("treatment", "", rcont$contrasts), gsub(".csv", "", maindata$inFile$name), "csv", sep = ".")
            else paste0(gsub("treatment", "", rcont$contrasts), ".csv")
        }, 
        content = function(fname){
            fwrite(rcont$contrast, fname, sep = ";")
        }
    )
    
    # output$downloadDemo <- downloadHandler(
    #     filename = function(){
    #         if(!is.null(maindata$inFile$name)) paste(gsub("condition", "", rcont$contrasts), gsub(".csv", "", maindata$inFile$name), "csv", sep = ".")
    #         else paste0(input$selectDemoData, ".csv")
    #     }, 
    #     content = function(fname){
    #         sample_1_exp <- "data/donors.uniprot.csv"
    #         fwrite(sample_1_exp, fname, sep = ";")
    #     }
    # )
    
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
            #paste(gsub("condition", "", rcont$contrasts), gsub(".csv", "", maindata$inFile$name), "html", sep = ".")
            paste("proteomill-report.html")
        },
        content = function(file) {
            
            cat(file=stderr(), paste0(getwd()))
            
            src <- normalizePath('reports/report.rmd')
            
            # temporarily switch to the temp dir, in case you do not have write
            # permission to the current working directory
            # owd <- setwd(tempdir())
            # on.exit(setwd(owd))
            
            # cat(file=stderr(), paste0(getwd()))
            
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
    
    # myData <- reactive({
    #     inFile <- input$file1
    #     if (is.null(inFile)) return(NULL)
    #     maindata$data_wide <- fread(inFile$datapath, header = T)
    #     
    #     maindata$data_wide
    # })

    

    observeEvent(input$file1, {
        
        inFile <- input$file1

        if (is.null(inFile)) return(NULL)
        #maindata$data_wide <- fread(inFile$datapath, header = T)
        
        dw <- tryCatch({fread(input=inFile$datapath);},
                        error = function(err){return(err)} ,
                        warning = function(war){return(war)} ,silent=F)
        
        # dw <- tryCatch({fread(input="files/E-PROT-45-query-results.tsv");},
        #                error = function(err){return(err)} ,
        #                warning = function(war){return(war)} ,silent=F)
        
        
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
            
            maindata$data_wide <- dw
            
            shinyjs::show("previewDTInfo")
            
            maindata$isValid <- validate_data_format(dw)
            
        }
        

    })

    
    output$previewDT <- renderDT({
        
        if(!is.null(maindata$data_wide)){
            
            dt <- maindata$data_wide
            
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
    
    
    observeEvent(input$EndStep1, {
        
        shinyjs::show("Modal1Spinner")
        
        if(input$organism == "Homo sapiens") {
            
            library(EnsDb.Hsapiens.v86)
            
            maindata$organism <- EnsDb.Hsapiens.v86
            
        } else if(input$organism == "Mus musculus") {
            
            library(EnsDb.Mmusculus.v79)
            
            maindata$organism <- EnsDb.Mmusculus.v79
            
        } else if(input$organism == "Rattus norvegicus") {
            
            library(EnsDb.Rnorvegicus.v79)
            
            maindata$organism <- EnsDb.Rnorvegicus.v79
            
        }
        
        
        toggleModal(session, "ImportModal1", toggle = "toggle")
        shinyjs::hide("Modal1Spinner")
        toggleModal(session, "ImportModal2", toggle = "toggle")
        
        
    })
    
    
    observeEvent(input$EndStep2, {
        
        
        
        if(is.null(maindata$data_wide)) {
            createAlert(session, "selectAFile", "exampleAlert",
                        content = "Please select a file.", append = FALSE,
                        style = "danger")
        }
        else {
            
            if(maindata$isValid) {
                
                shinyjs::show("Modal2Spinner")
                
                upload_data(i = "auto")
                
                toggleModal(session, "ImportModal2", toggle = "toggle")
                shinyjs::hide("Modal2Spinner")
                toggleModal(session, "ImportModal3", toggle = "toggle")
                
            } else {
                
                createAlert(session, "checkRequirements", "exampleAlert1",
                            content = "File format is not valid. Please make sure dataset fulfills checklist requirements.", append = F,
                            style = "info")
                
            }
            

            
        }
        
        
        # View(maindata$data_origin)
        # View(maindata$data_wide)
        # 
        
        
        # check identifiers
        
        
        
        # subset interactions
        
        
        

        
        
    })
    
    observeEvent(input$EndStep3, {
        
        # subset_interactions()
        
        toggleModal(session, "ImportModal3", toggle = "toggle")
        updateNotifications("Dataset uploaded.","check-circle", "success")
        updateTasks(text = "Upload a dataset", value = 100, color = "green", i = 0001)
        
        
    })
    
    # Heatmap Updated ----
    
    observeEvent(input$renderHeatmap, {
        
        dw <- maindata$data_wide
        
        dw_cor <- log2(dw[, -..convertColumns])
        dw_cor[is.na(dw_cor)] <- 0 # Impute
        dw_cor <- round(cor(dw_cor, method = input$corMethod), 2)
        
        output$heatmap <- renderPlotly({
            
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
        
    })
    
    # PCA Updated ----
    
    
    observeEvent(input$PCA, {
        
        dw <- maindata$data_wide
        si_condition <- sampleinfo$samples$treatment
        
        pca.data <- log2(dw[, -..convertColumns])
        pca.data[is.na(pca.data)] <- 0 # "Impute" missing values as 0
        pca.data <- IMIFA::pareto_scale(pca.data)
        
        p.pca <- prcomp(t(pca.data), center = TRUE, scale. = F)
        
        
        # if 2 dims
        
        if(input$pcaDims[2] - input$pcaDims[1] == 1) {
            
            mydf <- data.frame(pc1 = p.pca$x[, input$pcaDims[1]],
                               pc2 = p.pca$x[, input$pcaDims[2]],
                               si_condition)
            
            poly.df <- mydf %>% 
                group_by(si_condition) %>%
                do(.[chull(.$pc1, .$pc2),]) 
            
            if(input$showPolygons) {
                ggplot(mydf, aes(pc1, pc2, colour = as.factor(si_condition))) +
                    geom_polygon(data = poly.df, fill = "grey", alpha = .15) +
                    geom_point(size = 5) +
                    xlab(paste0("PC", input$pcaDims[1])) +
                    ylab(paste0("PC", input$pcaDims[2])) +
                    scale_color_brewer(palette = "Accent") +
                    theme_light()  -> p
            } else {
                ggplot(mydf, aes(pc1, pc2, colour = as.factor(si_condition))) +
                    geom_point(size = 5) +
                    scale_color_brewer(palette = "Accent") +
                    theme_light()  -> p
            }
            
            
            output$PCAplots <- renderPlotly({
                
                ggplotly(p)
                
            })
            
        } else if (input$pcaDims[2] - input$pcaDims[1] == 2) {
            
            p <- plotly::plot_ly(x = p.pca$x[,input$pcaDims[1]],
                                 y = p.pca$x[,(input$pcaDims[2] - 1)],
                                 z = p.pca$x[,input$pcaDims[2]],
                                 text = rownames(p.pca$x),
                                 hoverinfo = "text",
                                 color = sampleinfo$samples$treatment,
                                 #colors = c("red","green","blue"),
                                 colors = brewer.pal(length(unique(sampleinfo$samples$treatment)), "Accent")
            ) %>%
                plotly::add_markers(marker=list(size = 15, line=list(width = 1, color = "black"))) %>%
                plotly::layout(scene = list(xaxis = list(title = paste0("PC", input$pcaDims[1])),
                                            yaxis = list(title = paste0("PC", (input$pcaDims[2] - 1))),
                                            zaxis = list(title = paste0("PC", input$pcaDims[2]))))
            
            output$PCAplots <- renderPlotly({
                p
            })
            
        }
        
        
        output$scree <- renderPlot({
            
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
        
    })
    
    observeEvent(input$pcaDims, {
        
        print(input$pcaDims)
        
        if( (input$pcaDims[2] - input$pcaDims[1] > 2) | (input$pcaDims[2] - input$pcaDims[1] < 1)) {
            updateSliderInput(session = session, inputId = "pcaDims", value = c(input$pcaDims[1], (input$pcaDims[1] + 2)))
        }
        
    })
    
    observeEvent(input$checkMemory, {
        
        print(ls())
        
        
        
        #print(ls())
        
        #updateNotifications(as.numeric(mem_used() / 1000000),"info-circle", "info")
        
        # env <- environment()
        # 
        # print(ls(env))
        # 
        # data.frame(
        #     object = ls(env),
        #     size = unlist(lapply(ls(env), function(x) {
        #         object.size(get(x, envir = env, inherits = FALSE))
        #     }))
        # )
        
        
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
        
        test_res <- data.table(
            Test = c(
                "Number of columns",
                "Number of rows",
                "Duplicate columns",
                "Has separator",
                "Starts alpha-numeric",
                "Ends alpha-numeric",
                "Multiple treatments exist",
                "Multiple replicates exist"
            ),
            ErrMsg = c(
                "The dataset has too few (<5) or too many (>150) columns for ProteoMill's capacity.",
                "The dataset has too few (<5) or too many (>50.000) rows for ProteoMill's capacity.",
                "Duplicate column names is not allowed.",
                "Missing '_'-separator between treatment and replicate in column names.",
                "First character in column names must be alpha-numeric",
                "Last character in column names must be alpha-numeric",
                "Multiple treatments/groups required for differential analysis.",
                "Multiple replicates needed to estimate variance."
                
            ),
            Passed = c(
                test_col,
                test_row,
                test_uniq_cols,
                test_col_sep,
                test_starts_with,
                test_ends_with,
                test_mult_treat,
                test_mult_reps
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




