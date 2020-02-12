library(limma)
library(Biobase)
require(ggplot2)
require(ggrepel)
require(RColorBrewer)
require(dplyr)
library(plotly)
library(data.table)
require(AnnotationDbi)
library(EnsDb.Hsapiens.v86)

# Passing arguments from R Shiny to external jQuery scripts unfortunately requires som workarounds
# This function is actually necessary...
fade <- function(fadein) {
    return(fadein)
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


get_interactions <- function(){

  url <- "https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz"

  interactions <- data.table::fread(knee::collect(url),
                                    sep = ' ',
                                    header = T)
  return(interactions)
}



# Server ----
server <- function(session, input, output) {
    
    # Notifications ----
    
    # Initialize as empty list
    # Populate with updateNotifications() function
    notification_list <- list() 
    
    #task_list <- list(taskItem(text = "Upload a dataset", value = 0, color = "green"))
    
    tasks <- data.table(
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
    
    
    
    # Make notification_list global var
    assign("notification_list", notification_list, envir = .GlobalEnv)
    #assign("task_list", task_list, envir = .GlobalEnv)
    assign("tasks", tasks, envir = .GlobalEnv)
    
    # Initialize notification menu
    output$notifMenu <- renderMenu({
        dropdownMenu(type = "notifications", .list = notification_list)
    })
    
    # Initialize message menu
    output$helpMenu <- renderMenu({
        updateTasks("Run predictive analysis", value = 0, color = "green", i = 0008)
        dropdownMenu(type = "tasks", .list = task_list)
    })
    
    # A function to append task items to the menu
    updateTasks <- function(text, value, color, i) {
        
        i = sprintf("%04d", i)
        
        value <- ifelse(value > 99, 100, round(value, digits = 0))
        
        l = list(text = text, value = value, color = color, id = i)
        
        if(i %in% tasks$id){
            tasks[id == i, names(tasks) := l][]
        } else {
            # Create a new item
            
            tasks <- rbindlist(list(tasks, l))
            
        }

        task_list <- list()
        
        tmp <- tasks[value == 0]
        tmp <- tmp[order(id)][1:2]
        tmp <- tmp[!is.na(id)]
        tmp <- rbindlist(list(tasks[value > 0], tmp))
        
        
        for(x in 1:nrow(tmp)){
            item <- taskItem(text = tasks[x, text],
                             value = tasks[x, value],
                             color = tasks[x, color])
            
            # A hack to make the notification item clickable
            # Onclick opens a modal dialog
            item$children[[1]] <- a(
                "onclick" = paste0("clickFunction('", paste0(tasks[x, id]), "'); return false;"),
                list(item$children[[1]]$children))
            
            
            item <- list(item)
            
            task_list <- append(item, base::get("task_list"))
        }

        assign("task_list", task_list, envir = .GlobalEnv)
        assign("tasks", tasks, envir = .GlobalEnv)
        
        # Render the menu with appended notification list
        output$helpMenu <- renderMenu({
            
            dropdownMenu(type = "tasks", .list = task_list)
        })
        
    }
    
    
    
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
        
        notification_list <- append(item, base::get("notification_list", envir = .GlobalEnv))
        assign("notification_list", notification_list, envir = .GlobalEnv)
        
        # Render the menu with appended notification list
        output$notifMenu <- renderMenu({
            
            dropdownMenu(type = "notifications", .list = notification_list)
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
    
    observeEvent(input$updateAll, {
        data_wide <- data_wide[!duplicated(id_check[["updated"]]),]
        rownames(data_wide) <- id_check[["updated"]][!duplicated(id_check[["updated"]])]
        assign("data_wide", data_wide, envir = .GlobalEnv)
        
        background_data <- rownames(data_wide)
        input_names <- rownames(data_wide)
        assign("background_data", background_data, envir = .GlobalEnv)
        assign("input_names", input_names, envir = .GlobalEnv)
        
        message <- paste(nrow(id_check[["table"]]), " IDs updated successfully.")
        updateNotifications(message,"check-circle", "success")
        
    })
    
    observeEvent(input$sidebarmenu, {
        
        if(input$sidebarmenu %in% c("validateIDs", "structures")){
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
        
        if(input$sidebarmenu == "interactions" & !exists("v")){
            updateNotifications("Run pathway analysis first.","exclamation-triangle", "danger")
        } else if(input$sidebarmenu == "interactions" & exists("v")){
            updateTasks(text = "Run network analysis", value = 100, color = "green", i = 0007)
        }
    })
    
    # Clicked notification item
    observeEvent(input$linkClicked, {
        
        linkCode <- as.character(input$linkClicked)
        
        i <- as.character(substr(linkCode, nchar(linkCode)-3, nchar(linkCode)))
        
        item <- tasks[id == i]
        
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
    
    # Sidebar menu items ----
    
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
        
        updateSelectInput(session, "contrast1", choices = group)
        
    }
    
    observeEvent(input$displayIdentifier, {
        
        i <- identifier(input$displayIdentifier)
        
        assign("sID", i, envir = .GlobalEnv)
        
        updateNotifications(paste0("Default ID set to ", i),"info-circle", "info")
        
    }, ignoreInit = F)
    
    
    
    # File input ----
    
    # Main data file input
    observeEvent(input$infile, {
        
        inFile <- input$infile
        
        assign("inFile", inFile, envir = .GlobalEnv)
        
        if (is.null(inFile))
            return(NULL)

        upload_data(inFile$datapath, separator(input$dataSep), identifier(input$dataIdentiferType))

        updateNotifications("Dataset uploaded.","check-circle", "success")
        updateTasks(text = "Upload a dataset", value = 100, color = "green", i = 0001)
        
        
        
    })
    
    # Annotation file input
    observeEvent(input$anno_infile, {
        inFile <- input$anno_infile
        
        if (is.null(inFile))
            return(NULL)
        
        data_annotation <- data.table::fread(
            inFile$datapath,
            sep = separator(input$annoSep),
            dec = ".",
            header = T)
        
        assign('data_annotation', data_annotation, envir = .GlobalEnv)
        
        updateNotifications("Annotations uploaded.","check-circle", "success")
        updateTasks(text = "Upload annotation data", value = 100, color = "green", i = 0002)
        
        
    })
    
    # Demo data file input
    
    observeEvent(input$useDemoData, {
        
        # Sample 1
        
        sample_1_exp <- "data/donors.uniprot.csv"
        sample_1_anno <- "data/donors.uniprot.annotation.tsv"
        
        upload_data(sample_1_exp, "auto", "auto")
        
        data_annotation <- data.table::fread(
            sample_1_anno,
            sep = "auto",
            dec = ".",
            header = T)
        
        assign('data_annotation', data_annotation, envir = .GlobalEnv)

        updateNotifications("Demo data uploaded.","check-circle", "success")
        updateTasks(text = "Upload a dataset", value = 100, color = "green", i = 0001)
        updateTasks(text = "Upload annotation data", value = 100, color = "green", i = 0002)
        
    })
    
    # NA frequency plot ----
    
    renderNAfreq <- function(){
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
        
        if(exists("data_wide")){
            renderNAfreq()
        } else {
            updateNotifications("Upload a dataset first.","exclamation-triangle", "danger")
        }
    })
    
    observeEvent(input$setcutoff, {
        
        if(exists("data_wide")){
            data_wide <- filter_na(input$missingvalues)
            assign("data_wide", data_wide, envir = .GlobalEnv)
            
            unlock_menus()
            renderNAfreq()
            
            updateTasks(text = "Set a filter", value = 100, color = "green", i = 0003)
            updateNotifications(paste0("NA cutoff set to ", input$missingvalues, ".") ,"check-circle", "success")
        } else {
            updateNotifications("Upload a dataset first.","exclamation-triangle", "danger")
        }
        
    })
    
    
    
    # Distributions ----
    observeEvent(input$generatedistributions, {
        
        if(exists("data_wide")){
            
            data_wide_NAex <- na.exclude(dframe(data_origin, sID))
            
            fit.lnorm <- tryCatch( apply(data_wide_NAex, 1, function(x) fitdistrplus::fitdist(as.numeric(x), "lnorm")),
                                   error = function(e) print(e))
            
            fit.norm <- tryCatch( apply(data_wide_NAex, 1,  function(x) fitdistrplus::fitdist(as.numeric(x), "norm")),
                                  error = function(e) print(e))
            
            # Render "data type" distribution plots
            output$distributions <- renderPlot({
                
                updateSliderInput(session, inputId = "setdist", max = nrow(na.omit(data_origin)))
                
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
    
    # Quality control ----
    
    # Render PCA plots
    
    observeEvent(input$loadPCAplots, {
        output$pca2dplot <- renderPlot({
            c <- input$contribs
            e <- input$ellipse
            plotPCA(c,e, '2d')
        })
        
        output$pca3dplot <- plotly::renderPlotly({
            plotPCA(0,0, '3d')
        })
        
        updateTasks(text = "Inspect data", value = (tasks[id == "0004", value] + 100/3), color = "green", i = 0004)
    })
    

    
    # Render UMAP
    observeEvent(input$loadUMAP, {
        output$UMAPplot <- renderPlot({
            plotPCA(0,0, 'UMAP')
        })
        updateTasks(text = "Inspect data", value = (tasks[id == "0004", value] + 100/3), color = "green", i = 0004)
    })
    

    # Render heatmap
    
    observeEvent(input$generateheatmap, {
        cor_mat_raw_logged <- log2(data_origin[,-..convertColumns])
        cor_mat_raw_logged[is.na(cor_mat_raw_logged)] <- 0
        cor_mat_raw_logged <- cor(cor_mat_raw_logged)
        
        #annotation_col = data.frame(tmp)
        
        if(exists("data_annotation")){
            if(all(data_annotation[, as.character(.SD[[1L]])] %in% names(data_origin[,-..convertColumns]))){
                hmap <- pheatmap::pheatmap(
                    annotation_col = dframe(data_annotation, "V1")[2:4],
                    cor_mat_raw_logged,
                    legend_breaks = c(min(cor_mat_raw_logged), 1),
                    legend_labels = c(0, 1)
                )
            }
        } else {
            hmap <- pheatmap::pheatmap(
                cor_mat_raw_logged,
                legend_breaks = c(min(cor_mat_raw_logged), 1),
                legend_labels = c(0, 1)
            )
        }
        
        output$samplecorrheatmap = renderPlot({hmap})
        
        updateTasks(text = "Inspect data", value = (tasks[id == "0004", value] + 100/3), color = "green", i = 0004)
        
    })

    # Differential expression: set contrasts ----
    observeEvent(input$contrast1, {
        
        if (!exists("group", envir = .GlobalEnv)){
            return(NULL)
        } else if (exists("group", envir = .GlobalEnv)){
            cont1 <- input$contrast1
            cont2 <- group[group != cont1]
            
            updateSelectInput(session, "contrast2", choices = cont2)
        }
        
    })
    
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
                                 icon = shiny::icon("angle-double-right"), selected = F),
                     menuSubItem("Predictive analysis", tabName = "predictive", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F))
        })
        
        pairing <- input$pairing
        
        contrasts <- paste(input$contrast1,input$contrast2, sep = '-')
        
        contrast <- diff_exp(contrasts, pairing)
        assign('contrast', contrast, envir = .GlobalEnv)
        
        
        removeUI(selector = "#enrichrm")
        removeUI(selector = "#networkrm")
        removeUI(selector = "#interactionsrm")
        removeUI(selector = "#predictiverm")
        updateNotifications("A DE contrast has been set.","check-circle", "success")
        updateTasks(text = "Set contrast", value = 100, color = "green", i = 0005)
    })
    
    # Differential expression: output table
    
    observeEvent(input$loadDiffExpTable, {
        
        if(exists("contrast")){
            output$diffexptable <- DT::renderDataTable({
                
                df <- DT::datatable(dframe(contrast, sID),
                                    options = list(autoWidth = TRUE,
                                                   scrollX=TRUE))
                df %>% DT::formatSignif('adj.P.Val', digits = 2)
                
                
            })
        } else {
            updateNotifications("Set contrasts first.","exclamation-triangle", "danger")
        }
        
        
        
    })
    
    
    
    # Differential expression: confidence intervals
    
    output$contrasttable <- renderPlot({
        ggplot2::ggplot(cint, aes(x = protein, y = logFC, colour = logFC)) +
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
    })
    
    observeEvent(input$min_fc, {
      
      output$number_of_genes <- renderUI({
        
        up <- paste("Up-regulated:", contrast[logFC >= input$min_fc, .N], "genes", sep = " ")
        down <- paste("Down-regulated:", contrast[logFC < (input$min_fc * -1), .N], "genes", sep = " ")
        
        HTML(paste(up, down, '<br/>', sep = '<br/>'))
        
      })
    })
    
    observeEvent(input$generate_pathways, {
        
        UPREGULATED_genes <- contrast[logFC >= input$min_fc & adj.P.Val < 0.05, UNIPROTID]
        DOWNREGULATED_genes <- contrast[logFC < (input$min_fc * -1) & adj.P.Val < 0.05, UNIPROTID]
        
        UPREGULATED_pathways <- knee::ora(UPREGULATED_genes)@output
        assign("UPREGULATED_pathways", UPREGULATED_pathways, envir = .GlobalEnv)
        
        DOWNREGULATED_pathways<- knee::ora(DOWNREGULATED_genes)@output
        assign("DOWNREGULATED_pathways", DOWNREGULATED_pathways, envir = .GlobalEnv)
    
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
      
    })
    
    # Volcano plot
    
    observeEvent(input$loadPathwayPlots, {
        
        if(!exists("UPREGULATED_pathways")){
            updateNotifications("Run pathway analysis first.","exclamation-triangle", "danger")
        }else{
            
            output$volcano_plot <- renderPlotly(
                {
                    
                    
                    res <- knee::enrichment_results(UPREGULATED_pathways, DOWNREGULATED_pathways)
                    
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
                    
                    print(input$abstractionlevel)
                    v <- knee::volcano(res, abstraction = input$abstractionlevel)
                    
                    assign("v", v, envir = .GlobalEnv)
                    
                    plotly::ggplotly(v$volcano_plot) %>% layout(dragmode = "select")
                    
                })
            
            
            output$volcano_plot2 <- renderPlotly(
                {
                    plotly::ggplotly(v$volcano_plot, width = 700, height = 400) %>% layout(dragmode = "select")
                    
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
    
    
    
    
    output$xxxx <- renderVisNetwork({

        proteins <- contrast[adj.P.Val < input$pvaluecutoff & abs(logFC) >= input$fccutoff]$UNIPROTID

        # uniprot_to_string <- uniprot_to_string_src[up %in% proteins, c(3, 4, 6)]
        # 
        # 
        # if(!exists("g")){
        #     ints <- interactions[protein1 %in% uniprot_to_string$V3,]
        #     
        #     setkey(uniprot_to_string, V3)
        #     setkey(ints, protein2)
        #     ints <- ints[protein2 %in% uniprot_to_string$V3,]
        #     ints$combined_score <- ints$combined_score / 100
        #     
        #     setkey(uniprot_to_string, V3)
        #     setkey(ints, protein1)
        #     ints <- ints[uniprot_to_string, nomatch = 0]
        #     setkey(ints, protein2)
        #     ints <- ints[uniprot_to_string, nomatch = 0]
        #     ints <- ints[,c(5,7, 3)]
        #     colnames(ints) <- c("protein1", "protein2", "combined_score")
        # }
        

        
        mynodes <- unlist(event_data("plotly_selected")$key)
        
       # ints2 <- ints[combined_score > input$interactioncutoff]
        
        if(!exists("g")){
            ints2 <- interactions[(protein1 %in% proteins) & (protein2 %in% proteins)]
            ints2 <- ints2[score > input$interactioncutoff]
        }
        
        
        
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
        
        #edges$color <- palette(rainbow(7))[edges$value]
        
        g <- igraph::graph_from_data_frame(edges, directed = F, vertices = nodes)
        
        
        
        #g <- igraph::graph_from_data_frame(ints2, directed = T)
        r <- v$res[V(g)$name, on = "UNIPROTID", ]
        
        g <- set.vertex.attribute(g, name = "color",value = r$col)
        g <- set.vertex.attribute(g, name = "group",value = r$TopReactomeName)
        

        # g <- igraph::graph_from_data_frame(ints2, directed = F)
        # g <- igraph::simplify(g, remove.multiple = F, remove.loops = T)
        # 
        # r <- v$res[V(g)$name, on = "Gene", ]
        # 
        # g <- set.vertex.attribute(g, name = "color",value = r$col)
        # g <- set.vertex.attribute(g, name = "group",value = r$TopReactomeName)
        
        layout <- function(type) {
            switch(type,
                   "1" = "layout_nicely",
                   "2" = "layout_in_circle",
                   "3" = "layout_on_grid",
                   "4" = "layout_on_sphere",
                   "5" = "layout_randomly",
                   "6" = "layout_DH")
        }
        
        # visNetwork::visIgraph(g) %>%
        #     visIgraphLayout(layout = layout(input$network_layout_options), randomSeed = 1) %>%
        #     visEdges(arrows = list(to = list(enabled = TRUE))) %>%
        #     visOptions(selectedBy = list(variable = "group")) %>%
        #     visPhysics(stabilization = T)
        
        #  visNetwork(nodes, edges) %>% 
        #     # visEdges(arrows = list(to = list(enabled = TRUE))) %>%
        #      #visIgraphLayout(layout = layout("layout_nicely")) %>%
        #      visPhysics(enabled = FALSE)
        # # 
        visNetwork::visIgraph(g) %>%
            visIgraphLayout(layout = layout(input$network_layout_options), randomSeed = 1) %>%
            visOptions(selectedBy = list(variable = "group"))
        
        })
        

}
