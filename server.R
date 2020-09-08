# Load packages ----
require(shiny)
require(shinydashboard)
require(knitr)
require(limma)
require(Biobase)
require(ggplot2)
require(ggrepel)
require(RColorBrewer)
require(dplyr)
require(plotly)
require(data.table)
require(DT)
require(AnnotationDbi)
require(EnsDb.Hsapiens.v86)
require(networkD3)
require(XML)
require(mixOmics)
require(rhandsontable)
require(stringr)
require(factoextra)
require(pheatmap)
require(rmarkdown)
require(fitdistrplus)
require(igraph)
require(visNetwork)
require(R.utils)
require(umap)



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
    
    # Define reactive variables ----
    notifications <- reactiveValues()
    tasks <- reactiveValues()
    sampleinfo <- reactiveValues()
    maindata <- reactiveValues()
    rcont <- reactiveValues()
    pathways <- reactiveValues()
    
    sampleinfo$sID <- identifier(2)
    
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
        updateTasks(0, 0, 0, 0)
    }, ignoreNULL = T, ignoreInit = F)
    
    
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
        
        if(input$sidebarmenu == "about"){
            updateNotifications("You must gather your party before venturing forth.","exclamation-triangle", "danger")
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
                        # if(i == "0001"){
                        #     renderUI({
                        #     tags$iframe(src = paste0("https://qodb.shinyapps.io/generateDataset/"), width="100%", height="600px", frameborder="0")
                        #     })},
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
                     menuSubItem("UMAP", tabName = "UMAP", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F),
                     menuSubItem("Heatmap", tabName = "samplecorr", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F)
            )
        })
        output$differential <- renderMenu({
            menuItem("Differential analysis", class = 'btn-10', tabName = "differential", icon = icon("adjust"), href = NULL,
                     menuSubItem("Contrasts", tabName = "contrasts", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F),
                     menuSubItem("Differential expression", tabName = "differentialexpression", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F),
                     menuSubItem("Confidence intervals", tabName = "diffexpoutput", href = NULL, newtab = TRUE,
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
    
    
    
    
    # Build sample info ----
    
    #group <- list()
    sample_data <- function(data) {
        samples <- names(data[, -..convertColumns])
        condition <- as.factor(gsub('_.*', '', samples))
        replicate <- as.factor(gsub('.*_', '', samples))
        samples <- data.frame(samples, condition, replicate)
        rownames(samples) <- samples$samples
        
        group <- sapply(as.character(unique(condition)), function(x) paste("condition", x, sep = ''))
        
        sampleinfo$samples <- samples
        sampleinfo$condition <- condition
        sampleinfo$replicate <- replicate
        
        sampleinfo$group <- group
        
        updateContrasts()
        
        return (FALSE)
    }
    
    subsample_data <- function(){
        samples <- sampleinfo$samples
        condition <- samples$condition
        
        group <- sapply(as.character(unique(condition)), function(x) paste("condition", x, sep = ''))
        sampleinfo$group <- group
    }
    
    
    # Filter NA ----
    filter_na <- function(data_origin, threshold) {

        # Which elements are NA?
        allNA <- is.na(data_origin[, -..convertColumns])

        # Summary of how many TRUEs there are in each row
        NA_frequency <- table(rowSums(allNA))

        # Subset to NA threshold ----
        
        subset_NA <- function(condition)
        {
            # Subset columns by condition
            condition_subset <- data_origin[, grep(condition,names(data_origin)), with = F]
            
            # Determine if rows pass NA threshold
            rows_to_keep <- rowSums(is.na(condition_subset)) <= threshold
            
            # Subset rownames
            keep <- data_origin[rows_to_keep, UNIPROTID]
            
            return(keep)
            
        }
        
        # assign("group", sampleinfo$group, envir = .GlobalEnv)
        # assign("subset_NA", subset_NA, envir = .GlobalEnv)
        
        # Apply function to all regions
        condition_sub <- lapply(names(sampleinfo$group), subset_NA)
        
        # Reduce to shared proteins
        condition_sub2 <- Reduce(intersect, condition_sub)
        
        d <- data_origin[UNIPROTID %in% condition_sub2,]
        
        return(d)
        
    }
    
    
    # File input: upload data function ----
    
    undup <- function(genes){
        genes[!is.na(genes)]
    }
    
    upload_data <- function(path, sep, i){
        
        maindata$data_wide <- data.table::fread(
            path,
            sep = sep,
            dec = ".",
            header = T)
        
        maindata$data_wide <- maindata$data_wide[!duplicated(names(maindata$data_wide)[1])]
        #maindata$data_wide <- maindata$data_wide[apply(maindata$data_wide[, 2:ncol(maindata$data_wide)], 1, function(x) sum(x, na.rm = T) > 500)]
        
        for(j in seq_along(maindata$data_wide)){
            set(maindata$data_wide, i = which(maindata$data_wide[[j]] == 0 & is.numeric(maindata$data_wide[[j]])), j = j, value = NA)
        }
        
        empty_rows <- apply(maindata$data_wide[,2:ncol(maindata$data_wide)], 1, function(x) all(is.na(x)))
        maindata$data_wide <- maindata$data_wide[!empty_rows,]
        
        maindata$data_origin <- maindata$data_wide
        
        convertColumns <- c("UNIPROTID", "ENTREZID", "SYMBOL", "GENEID", "PROTEINID")
        
        assign('convertColumns', convertColumns, envir = .GlobalEnv)
        
        keys <- maindata$data_wide[, as.character(.SD[[1L]])][1:10]
        
        if(i == "auto") {
            
            tr1 <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = keys, column = "SYMBOL", keytype = "UNIPROTID", multiVals = "first")
            tr2 <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = keys, column = "UNIPROTID", keytype = "ENTREZID", multiVals = "first")
            tr3 <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = keys, column = "UNIPROTID", keytype = "SYMBOL", multiVals = "first")
            tr4 <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = keys, column = "UNIPROTID", keytype = "GENEID", multiVals = "first")
            tr5 <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = keys, column = "UNIPROTID", keytype = "PROTEINID", multiVals = "first")
            
            trs <- list(tr1[!is.na(tr1)],
                        tr2[!is.na(tr2)],
                        tr3[!is.na(tr3)],
                        tr4[!is.na(tr4)],
                        tr5[!is.na(tr5)])
            
            i <- convertColumns[which.max(lapply(lapply(trs, lengths), sum))]
            
            if(i == "UNIPROTID"){
                
                tr_all <- data.table(AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = maindata$data_wide[, as.character(.SD[[1L]])], columns = convertColumns, keytype = i))
                tr_all <- tr_all[!duplicated(UNIPROTID)]
                
            } else {
                tr <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = maindata$data_wide[, as.character(.SD[[1L]])], column = "UNIPROTID", keytype = i, multiVals = "first")
                tr <- tr[!is.na(tr)]
                tr <- tr[!duplicated(tr)]
                
                tr <- data.table(i = names(tr), "UNIPROTID" = unname(tr))
                
                tr_all <- data.table(AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = maindata$data_wide[, as.character(.SD[[1L]])], columns = convertColumns, keytype = i))
                tr_all <- tr_all[UNIPROTID %in% tr$UNIPROTID]
                tr_all <- tr_all[!duplicated(UNIPROTID)]
            }
            
            
            
            
        } else {
            
            if(i == "UNIPROTID"){
                
                tr_all <- data.table(AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = maindata$data_wide[, as.character(.SD[[1L]])], columns = convertColumns, keytype = i))
                tr_all <- tr_all[!duplicated(UNIPROTID)]
                
            } else {
                tr <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = maindata$data_wide[, as.character(.SD[[1L]])], column = "UNIPROTID", keytype = i, multiVals = "first")
                tr <- tr[!is.na(tr)]
                tr <- tr[!duplicated(tr)]
                
                tr <- data.table(i = names(tr), "UNIPROTID" = unname(tr))
                
                tr_all <- data.table(AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = maindata$data_wide[, as.character(.SD[[1L]])], columns = convertColumns, keytype = i))
                tr_all <- tr_all[UNIPROTID %in% tr$UNIPROTID]
                tr_all <- tr_all[!duplicated(UNIPROTID)]
            }
        }
        
        names(maindata$data_wide)[1] <- i
        
        setkeyv(tr_all, i)
        
        setkeyv(maindata$data_wide, i)
        
        tr_all <- tr_all[maindata$data_wide[, ..i], on = i]
        
        maindata$data_wide <- maindata$data_wide[tr_all, nomatch = 0]
        
        maindata$data_wide$ENTREZID <- as.character(maindata$data_wide$ENTREZID)
        
        for(j in seq_along(maindata$data_wide)){
            set(maindata$data_wide, i = which(duplicated(maindata$data_wide[[j]]) & is.character(maindata$data_wide[[j]])), j = j, value = NA)
        }
        
        maindata$data_wide[is.na(ENTREZID), ENTREZID := paste0("MISSING_", seq(1:length(is.na(ENTREZID))))]
        maindata$data_wide[is.na(SYMBOL), SYMBOL := paste0("MISSING_", seq(1:length(is.na(SYMBOL))))]
        maindata$data_wide[is.na(UNIPROTID), UNIPROTID := paste0("MISSING_", seq(1:length(is.na(UNIPROTID))))]
        maindata$data_wide[is.na(GENEID), GENEID := paste0("MISSING_", seq(1:length(is.na(GENEID))))]
        maindata$data_wide[is.na(PROTEINID), PROTEINID := paste0("MISSING_", seq(1:length(is.na(PROTEINID))))]
        
        setcolorder(maindata$data_wide, c(convertColumns, names(maindata$data_origin[,2:ncol(maindata$data_origin)])))
        
        maindata$data_origin <- maindata$data_wide
        

        setkey(maindata$data_wide, "UNIPROTID")
        setkey(pdesc, "UNIPROTID")
        
        pdesc <- pdesc[maindata$data_wide[,1:5]]
        maindata$pdesc <- pdesc
        
        sample_data(maindata$data_wide)
        
        #assign("thiss", maindata$data_wide, envir = .GlobalEnv)
        
    }
    
    # File input: Main data ----
    observeEvent(input$infile, {
        
        maindata$inFile <- input$infile
        
        if (is.null(maindata$inFile))
            return(NULL)

        upload_data(maindata$inFile$datapath, separator(input$dataSep), identifier(input$dataIdentiferType))

        updateNotifications("Dataset uploaded.","check-circle", "success")
        updateTasks(text = "Upload a dataset", value = 100, color = "green", i = 0001)
        
        
        
    })
    
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
        sample_1_anno <- "data/donors.uniprot.annotation.tsv"
        
        upload_data(path = sample_1_exp, sep = ";", i = "UNIPROTID")
        
        maindata$data_annotation <- data.table::fread(
            sample_1_anno,
            sep = "auto",
            dec = ".",
            header = T)
        
        updateNotifications("Demo data uploaded.","check-circle", "success")
        updateTasks(text = "Upload a dataset", value = 100, color = "green", i = 0001)
        updateTasks(text = "Upload annotation data", value = 100, color = "green", i = 0002)
        
    })
    
    
    
    # Data summary plot ----
    
    output$datainfoBox <- renderInfoBox({
        if(!is.null(sampleinfo$samples)) infoBox("Genes", maindata$data_wide[, .N])
        else infoBox("Genes", 0)
    })
    
    output$sampleinfoBox <- renderInfoBox({
        if(!is.null(sampleinfo$samples)) infoBox("Samples", nrow(sampleinfo$samples))
        else infoBox("Samples", 0)
    })
    
    output$conditioninfoBox <- renderInfoBox({
        if(!is.null(sampleinfo$samples)) infoBox("Treatments", length(unique(sampleinfo$samples$condition)))
        else infoBox("Treatments", 0)
    })
    
    observe({
        if(!is.null(input$hot)){
            sampleinfo$samples <- as.data.frame(hot_to_r(input$hot))
            output$hot <- renderRHandsontable({
                rhandsontable(sampleinfo$samples, width = 600, rowHeaders = NULL)
            })
        }
    })    
    
    output$hot <- renderRHandsontable({
        if(!is.null(sampleinfo$samples)) rhandsontable(sampleinfo$samples, width = 600, rowHeaders = NULL)
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
            
            if(nrow(sampleinfo$samples) > 30) {
                
                ggplot(dat, 
                       aes(x = sample, y = values, fill = sample)) + 
                    geom_violin(trim = FALSE, scale = "width") +
                    geom_boxplot(width=0.1, fill="white") +
                    # labs(title ="Distribution of M2 Macrophages", 
                    #      x = "Tissue Samples", y = "Cibersort Count") +
                    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.position = "none") +
                    xlab("Samples") + ylab("Log2 Expression") +
                    scale_fill_brewer(palette = "BuGn")
                
            } else {
                
                ggplot(dat, 
                       aes(x = ind, y = values, fill = sample)) + 
                    geom_violin(trim = FALSE, scale = "width") +
                    geom_boxplot(width=0.1, fill="white") +
                    # labs(title ="Distribution of M2 Macrophages", 
                    #      x = "Tissue Samples", y = "Cibersort Count") +
                    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.position = "none") +
                    xlab("Samples") + ylab("Log2 Expression") +
                    scale_fill_brewer(palette = "BuGn")
                
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
            
            theme_set(theme_classic())
            
            
            
            # Draw plot
            ggplot(naf, aes(x = Var1, y = Freq)) +
                geom_bar(stat = "identity", width = .9, fill = "tomato3") +
                labs(x = "Number of missing values in at least one sample", y = "Number of rows") +
                theme(axis.text.x = element_text(angle = 65, vjust = 0.6))
            
            #graphics::barplot(NA_frequency, xlab = "Na frequency", ylab = "Number of genes")
            
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
            data_wide <- filter_na(data_origin, input$missingvalues)
            unlock_menus()
            renderNAfreq(data_wide)
            
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
        
        dt <- dframe(maindata$data_origin, sampleinfo$sID)
        
        # Biplot extension displaying top contributing proteins currently only available for 2D plot.
        
        if(type == '2d') {
            
            pca.data <- log2(dt)  # Log2 transform data
            pca.data[is.na(pca.data)] <- 0 # "Impute" missing values as 0
            
            pca.data <- t(pca.data)        # Transpose dataset
            
            p.pca <- prcomp(pca.data, center = TRUE, scale. = TRUE)
            
            pcaplot <- factoextra::fviz_pca_biplot(p.pca, title = '', label = "var", habillage = sampleinfo$samples$condition,
                                                   addEllipses = TRUE, ellipse.level = ellipse,
                                                   select.var = list(contrib = contribs), repel = TRUE)
            
            return (pcaplot)
            
        } else if (type == '2dpaired') {
            
            assign("dt", dt, envir = .GlobalEnv)
            
            pca.data <- dt
            pca.data[is.na(pca.data)] <- .5 # "Impute" missing values as 0

            pca.data <- pca.data[rowSums(pca.data) != .5,]
            
            pca.data <- t(pca.data)        # Transpose dataset
            
            pca.res = pca(X = pca.data,
                          multilevel = sampleinfo$samples$replicate, logratio = 'CLR')
            
            pcaplot <- plotIndiv(
                pca.res,
                ind.names = sampleinfo$samples$replicate,
                group = sampleinfo$samples$condition,
                title = "Multilevel PCA",
                legend = T,
                style = "ggplot2",
                ellipse = FALSE,
                ellipse.level = .9
            )
            return(pcaplot)
            
            
            
        } else if (type == '3d') {
            
            pca.data <- log2(dt)  # Log2 transform data
            pca.data[is.na(pca.data)] <- 0 # "Impute" missing values as 0
            
            pca.data <- t(pca.data)        # Transpose dataset
            
            p.pca <- prcomp(pca.data, center = TRUE, scale. = TRUE)
            
            pcaplot <- plotly::plot_ly(x = p.pca$x[,1],
                                       y = p.pca$x[,2],
                                       z = p.pca$x[,3],
                                       color = sampleinfo$samples$condition,
                                       colors = c("red","green","blue"),
                                       sizes = c(100, 150)) %>%
                plotly::add_markers() %>%
                plotly::layout(scene = list(xaxis = list(title = 'PC1'),
                                            yaxis = list(title = 'PC2'),
                                            zaxis = list(title = 'PC3')))
            
            return (pcaplot)
        } else if (type == '3dpaired') {
            
            pca.data <- dt
            pca.data[is.na(pca.data)] <- .5 # "Impute" missing values as 0
            
            pca.data <- pca.data[rowSums(pca.data) != .5,]
            
            pca.data <- t(pca.data)        # Transpose dataset
            
            pca.res = pca(X = pca.data,
                          ncomp = 3,
                          multilevel = sampleinfo$samples$replicate, logratio = 'CLR')
            
            pcaplot <- plotly::plot_ly(x = pca.res$x[,1],
                                       y = pca.res$x[,2],
                                       z = pca.res$x[,3],
                                       color = sampleinfo$samples$condition,
                                       colors = c("red","green","blue"),
                                       sizes = c(100, 150)) %>%
                plotly::add_markers() %>%
                plotly::layout(scene = list(xaxis = list(title = 'PC1'),
                                            yaxis = list(title = 'PC2'),
                                            zaxis = list(title = 'PC3')))

            return(pcaplot)
            
        } else if (type == 'UMAP') {
            
            pca.data <- log2(dt)  # Log2 transform data
            pca.data[is.na(pca.data)] <- 0 # "Impute" missing values as 0
            
            pca.data <- t(pca.data)        # Transpose dataset

            um <- umap::umap(pca.data, n_neighbors = ncol(dt))
            
            df <- data.frame(x = um$layout[,1],
                             y = um$layout[,2],
                             Sample <- sampleinfo$samples$condition)
            
            ggplot(df, aes(x, y, colour = sampleinfo$samples$condition, shape = sampleinfo$samples$condition)) +
                geom_point(size = 4)
            
            
            
        } else { return (FALSE) }
        
    }
    
    # Render PCA plots
    
    observeEvent(input$loadPCAplots, {
        
        # If unpaired
        if(!any(duplicated(sampleinfo$samples$replicate))){
            output$pca2dplot <- renderPlot({
                plotPCA(input$contribs, input$ellipse, '2d')
            })
        } else {
            output$pca2dplot <- renderPlot({
                plotPCA(0, input$ellipse, '2dpaired')

            })
        }

        # If unpaired
        if(!any(duplicated(sampleinfo$samples$replicate))){
            output$pca3dplot <- plotly::renderPlotly({
                plotPCA(0,0, '3d')
            })
        } else {
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
                hmap <- pheatmap::pheatmap(
                    annotation_col = dframe(data_annotation, "V1")[2:4],
                    cor_mat_raw_logged,
                    legend_breaks = c(min(cor_mat_raw_logged), 1),
                    legend_labels = c(0, 1),
                    silent = T
                )
            }
        } else {
            hmap <- pheatmap::pheatmap(
                cor_mat_raw_logged,
                legend_breaks = c(min(cor_mat_raw_logged), 1),
                legend_labels = c(0, 1),
                silent = T
            )
        }
        
        output$samplecorrheatmap = renderPlot({
            hmap
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
        
        #pairing <- input$pairing
        
        pairing <- ifelse(any(duplicated(sampleinfo$samples$replicate)), 1, 2)
        
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
        
        best_fit = 'normal'
        
        if(best_fit == 'nbinom') {
            dds <- DESeqDataSetFromMatrix(countData  = maindata$data_wide,
                                          colData    = sampleinfo$samples,
                                          design     = ~ sampleinfo$condition + sampleinfo$replicate)
            dds <- DESeq(dds)
        }
        
        # Create Annotation data and expression set (Biobase)
        phenoData <- new("AnnotatedDataFrame", data = sampleinfo$samples)
        exampleSet <- ExpressionSet(assayData = as.matrix(log2(dframe(maindata$data_wide, sampleinfo$sID))), phenoData = phenoData)
        
        condition <- sampleinfo$condition
        replicate <- sampleinfo$replicate
        
        unpaired <- model.matrix( ~ 0 + condition )
        paired <- model.matrix( ~ 0 + condition + replicate )
        
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
        names(cc) <- gsub('-','', gsub('condition','',cc))
        
        cont.matrix <- makeContrasts(contrasts = cc, levels = design) # All possible contrasts
        
        # Contrast groups, run empirical bayes statistics
        fit.cont <- contrasts.fit(fit, cont.matrix)
        fit.cont <- eBayes(fit.cont, robust = T)
        
        # Generate data frame with results from linear model fit, with confidence intervals.
        contrast <- toptable(fit.cont, number = Inf, coef = coeff, confint = TRUE)
        
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
        setcolorder(contrast, c(convertColumns, "logFC", "CI.L", "CI.R", "t", "P.Value", "adj.P.Val", "B"))
        
        assign("contrast", contrast, envir = .GlobalEnv)
        
        return( contrast )
        
    }
    

    output$diffexptable_up <- DT::renderDataTable({
        
        contrast <- rcont$contrast
        
        if(!is.null(contrast)){
            contrast <- contrast[order(logFC, decreasing = T)]
            df <- DT::datatable(dframe(rcont$contrast, sampleinfo$sID),
                                options = list(autoWidth = TRUE,
                                               scrollX=TRUE,
                                               order = list(1, 'desc'))) %>% 
                formatRound(columns=c(1, 2, 3, 4, 5, 6, 7), digits=4)
            # df %>% DT::formatSignif('logFC', digits = 2)
            # df %>% DT::formatSignif('CI.L', digits = 2)
            # df %>% DT::formatSignif('CI.R', digits = 2)
            # df %>% DT::formatSignif('adj.P.Val', digits = 2)
            
        }

    })
    
    output$diffexptable_down <- DT::renderDataTable({
        
        contrast <- rcont$contrast
        
        if(!is.null(contrast)){
            contrast <- contrast[order(logFC, decreasing = F)]
            df <- DT::datatable(dframe(rcont$contrast, sampleinfo$sID),
                                options = list(autoWidth = TRUE,
                                               scrollX=TRUE,
                                               order = list(1, 'asc'))) %>% 
                formatRound(columns=c(1, 2, 3, 4, 5, 6, 7), digits=4)
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
          
        up <- paste("Up-regulated:", contrast[logFC >= input$min_fc & adj.P.Val < input$min_pval, .N], "genes", sep = " ")
        down <- paste("Down-regulated:", contrast[logFC < (input$min_fc * -1) & adj.P.Val < input$min_pval, .N], "genes", sep = " ")
        
        HTML(paste(up, down, '<br/>', sep = '<br/>'))
        
      })
    })
    
    observeEvent(input$generate_pathways, {
        
        contrast <- rcont$contrast
        
        if(contrast[logFC >= input$min_fc & adj.P.Val < input$min_pval, .N] < 3 | contrast[logFC >= input$min_fc & adj.P.Val < input$min_pval, .N] < 3){
            updateNotifications("Too few differential genes.","exclamation-triangle", "danger")
        } else {
            
            updateNumericInput(session = session, inputId = "pvaluecutoff", label = "Maximum adj. Pvalue", value = input$min_pval)
            updateNumericInput(session = session, inputId = "fccutoff", label = "Minimum abs. log2FC", value = input$min_fc)
            
            UPREGULATED_genes <- contrast[logFC >= input$min_fc & adj.P.Val < input$min_pval, UNIPROTID]
            DOWNREGULATED_genes <- contrast[logFC < (input$min_fc * -1) & adj.P.Val < input$min_pval, UNIPROTID]
            
            UPREGULATED_pathways <- ora(UPREGULATED_genes)@output
            pathways$UPREGULATED_pathways <- UPREGULATED_pathways
            
            DOWNREGULATED_pathways<- ora(DOWNREGULATED_genes)@output
            pathways$DOWNREGULATED_pathways <- DOWNREGULATED_pathways
            
            output$upregulated_pathways_table <- DT::renderDT(
                
                DT::datatable(UPREGULATED_pathways[, -c("genes", "background")],
                              options = list(autoWidth = TRUE,
                                             scrollX=TRUE,
                                             columnDefs = list(
                                                 list(width = '100px', targets = c(1, 3)),
                                                 list(width = '60px', targets = c(6, 7))
                                             ))), server = F
                
            )
            
            output$downregulated_pathways_table <- DT::renderDT(
                DT::datatable(DOWNREGULATED_pathways[, -c("genes", "background")],
                              options = list(autoWidth = TRUE,
                                             scrollX=TRUE,
                                             columnDefs = list(
                                                 list(width = '100px', targets = c(1, 3)),
                                                 list(width = '60px', targets = c(6, 7))
                                             ))),
                server = F
            )
            
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

                    tba <- contrast[!(UNIPROTID %in% res$Gene), c("UNIPROTID", "logFC", "CI.L", "CI.R", "P.Value")]
                    colnames(tba) <- c("Gene", "logFC", "CI.L", "CI.R", "P.Value")
                    
                    tba <- tba[!is.na(logFC)]
                    
                    tba[, c("ReactomeID", "P.Adj")] <- NA
                    
                    tba[, "TopReactomeName"] <- "[No significant over-representation]"
                    tba[, "Pathway_name"] <- "[No significant over-representation]"
                    
                    setcolorder(tba, colnames(res))
                    
                    res <- rbindlist(list(res, tba))
                    
                    res$Pathway_name <- stringr::str_wrap(res$Pathway_name, 50)
                    
                    setkeyv(res, "Gene")
                    setkeyv(contrast, "UNIPROTID")
                    res <- res[contrast[,..convertColumns], nomatch = 0]
                    names(res) <- c("UNIPROTID", names(tba[,2:ncol(tba)]), convertColumns[-1])
                    
                    setcolorder(res, c(convertColumns, names(tba[,2:ncol(tba)])))

                    pathways$v <- volcano(res, abstraction = input$abstractionlevel, sID)

                    plotly::ggplotly(pathways$v$volcano_plot) %>%
                        layout(dragmode = "select")

                })

            
            output$volcano_plot2 <- renderPlotly(
                {
                    plotly::ggplotly(pathways$v$volcano_plot, width = 700, height = 500) %>%
                        layout(dragmode = "select") %>% 
                        config(scrollZoom = T)
                    
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
        
        assign("res", pathways$v$res, envir = .GlobalEnv)

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

        
        sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
                      Target = "target", Value = "value", NodeID = "name",
                      fontFamily = "sans-serif",
                      fontSize = 10, nodeWidth = 60, sinksRight = F)

        
    })
    
    output$xxxx <- renderVisNetwork({

        contrast <- rcont$contrast
        res <- pathways$v$res
        sID <- sampleinfo$sID
        
        dir <- input$network_regulation

        if(dir == 1){
            proteins <- contrast[(adj.P.Val <= input$pvaluecutoff) &
                                     (abs(logFC) >= input$fccutoff) &
                                     (logFC > 0)]$UNIPROTID
        } else if(dir == 2){
            proteins <- contrast[(adj.P.Val <= input$pvaluecutoff) &
                                     (abs(logFC) >= input$fccutoff) &
                                     (logFC <= 0)]$UNIPROTID
        } else {
            proteins <- contrast[(adj.P.Val <= input$pvaluecutoff) &
                                     (abs(logFC) >= input$fccutoff)]$UNIPROTID
        }

        # if(!exists("g")){
        #     ints2 <- interactions[(protein1 %in% proteins) & (protein2 %in% proteins)]
        #     ints2 <- ints2[score > input$interactioncutoff]
        # }
        
        ints2 <- interactions[(protein1 %in% proteins) & (protein2 %in% proteins)]
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
                #label = ints2$mode
                
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

                visNetworkProxy("xxxx") %>%
                    visUnselectAll() %>%
                    visSelectNodes(id = selection)
                
            }
            
            else {
                # Single item
                selection <- as.character(input$network_proteins)

                visNetworkProxy("xxxx") %>%
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
        
        if(exists("maindata$data_wide")){
            
            maindata$data_wide_NAex <- na.exclude(dframe(maindata$data_origin, sampleinfo$sID))
            
            fit.lnorm <- tryCatch( apply(maindata$data_wide_NAex, 1, function(x) fitdistrplus::fitdist(as.numeric(x), "lnorm")),
                                   error = function(e) print(e))
            
            fit.norm <- tryCatch( apply(maindata$data_wide_NAex, 1,  function(x) fitdistrplus::fitdist(as.numeric(x), "norm")),
                                  error = function(e) print(e))
            
            # Render "data type" distribution plots
            output$distributions <- renderPlot({
                
                updateSliderInput(session, inputId = "setdist", max = nrow(na.omit(maindata$data_origin)))
                
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
        
        if(!exists("data_wide")){
            
            updateNotifications(paste0("Upload a dataset first."), "exclamation-triangle", "danger")
            
        } else{
            
            updateNotifications(paste0("Checking for outdated IDs. Please wait."), "info-circle", "info")
            
            candidates <- data_wide[apply(
                data_wide[, ..convertColumns],
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
            
        }
        
        
        
    })
    
    # Generate report ----
    
    observeEvent(input$selectall, {

        updateCheckboxGroupInput (session, "export_filtering", "Filtering", selected = c("Missing values"))
        updateCheckboxGroupInput (session, "export_data_inspection", "Data inspection", selected = c("PCA 2D", "PCA 3D", "UMAP", "Heatmap"))
        updateCheckboxGroupInput (session, "export_de", "Differential analysis", selected = c("Fold-change", "P-values"))
        
        
        
    })
    
    observeEvent(input$deselectall, {
        
        updateCheckboxGroupInput(session, "export_filtering", "Filtering", inline = T, choices = c("Missing values"))
        updateCheckboxGroupInput(session, "export_data_inspection", "Data inspection", inline = T, choices = c("PCA 2D", "PCA 3D", "UMAP", "Heatmap"))
        updateCheckboxGroupInput(session, "export_de", "Differential analysis", inline = T, choices = c("Fold-change", "P-values"))
        
        
    })
    
    
    getOverview <- reactive({
        data.table(
            "Variable" = c("Total proteins", "DE proteins (adj. P. < 0.05)"),
            "Value" = c(maindata$data_origin[, .N], contrast[adj.P.Val < 0.05, .N])
            
        )
        
    })
    
    
    getDT <- reactive({
        maindata$data_wide
        
    })
    
    
    getpca2d <- reactive({
        pca2d
    })
    
    
    getUpPathways <- reactive({
        UPREGPATH <- UPREGULATED_pathways[, -c("ReactomeID", "background")]
        data.table::setcolorder(UPREGPATH, c("Pathway_name", "TopReactomeName", "q", "m", "p", "p.adj", "genes"))
        UPREGPATH
    })
    
    getDownPathways <- reactive({
        DOWNREGPATH <- DOWNREGULATED_pathways[, -c("ReactomeID", "background")]
        data.table::setcolorder(DOWNREGPATH, c("Pathway_name", "TopReactomeName", "q", "m", "p", "p.adj", "genes"))
        DOWNREGPATH
    })
    
    getVolcano <- reactive({
        res <- enrichment_results(UPREGULATED_pathways, DOWNREGULATED_pathways)
        
        tba <- contrast[!(UNIPROTID %in% res$Gene), c("UNIPROTID", "logFC", "CI.L", "CI.R", "P.Value")]
        colnames(tba) <- c("Gene", "logFC", "CI.L", "CI.R", "P.Value")
        
        tba <- tba[!is.na(logFC)]
        
        tba[, c("ReactomeID", "P.Adj")] <- NA
        
        tba[, "TopReactomeName"] <- "[No significant over-representation]"
        tba[, "Pathway_name"] <- "[No significant over-representation]"
        
        setcolorder(tba, colnames(res))
        
        res <- rbindlist(list(res, tba))
        
        res$Pathway_name <- stringr::str_wrap(res$Pathway_name, 50)
        
        setkeyv(res, "Gene")
        setkeyv(contrast, "UNIPROTID")
        res <- res[contrast[,..convertColumns], nomatch = 0]
        names(res) <- c("UNIPROTID", names(tba[,2:ncol(tba)]), convertColumns[-1])
        
        setcolorder(res, c(convertColumns, names(tba[,2:ncol(tba)])))

        v <- volcano(res, abstraction = input$abstractionlevel)
        
    })
    
    getReg <- reactive({
        data.table("Upregulated" = contrast[logFC <= 0 & adj.P.Val < 0.05, .N], "Downregulated" = contrast[logFC > 0 & adj.P.Val < 0.05, .N])
        
    })


    

    #reactiveFunction <- reactive({ ggplotly(v$volcano_plot) })
    
    # output$reactiveTable <- renderDataTable({ reactiveFunction() }, rownames = FALSE)
    # 
    # output$whatever <- renderUI({
    #     dataTableOutput("reactiveTable")
    # })
    
    
    output$download <- downloadHandler(
        filename = function(){
            if(!is.null(maindata$inFile$name)) paste(gsub("condition", "", rcont$contrasts), gsub(".csv", "", maindata$inFile$name), "csv", sep = ".")
            else paste0(gsub("condition", "", rcont$contrasts), ".csv")
        }, 
        content = function(fname){
            fwrite(rcont$contrast, fname, sep = ";")
        }
    )
    
    
    
    output$downloadReport <- downloadHandler(
        filename = function() {
            paste(gsub("condition", "", rcont$contrasts), gsub(".csv", "", maindata$inFile$name), "html", sep = ".")
        },
        content = function(file) {
            src <- normalizePath('report_file.Rmd')
            
            # temporarily switch to the temp dir, in case you do not have write
            # permission to the current working directory
            owd <- setwd(tempdir())
            on.exit(setwd(owd))
            file.copy(src, 'report_file.Rmd', overwrite = TRUE) 
            
            out <- rmarkdown::render('report_file.Rmd')
            file.rename(out, file)
        }
    )
    

    
    
}
