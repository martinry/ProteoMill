library(limma)
library(Biobase)
require(ggplot2)
require(ggrepel)
require(RColorBrewer)
require(dplyr)


# Server ----
server <- function(session, input, output) {
    
    # Notifications ----
    
    # Initialize as empty list
    # Populate with updateNotifications() function
    notification_list <- list() 
    
    # Make notification_list global var
    assign("notification_list", notification_list, envir = .GlobalEnv)
    
    # Initialize notification menu
    output$notifMenu <- renderMenu({
        dropdownMenu(type = "notifications", .list = notification_list)
    })
    
    # A function to append notification items to the menu
    updateNotifications <- function(notif, icon, status) {
        
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
        
        message <- paste(nrow(id_check[["table"]]), " updated successfully.")
        updateNotifications(message,"check-circle", "success")
    })
    
    # Clicked notification item
    observeEvent(input$linkClicked, {
        
        showModal(
            modalDialog(title = "Obsolete IDs",
                        
                        actionButton("updateAll", "Update all"),
                        
                        
                        DT::renderDataTable({
                            id_check[["table"]]
                        }),
                        size = "l"
            )
        )
    })
    
    # Sidebar menu items ----
    
    # Render "locked" menu items
    output$qualityrm <- renderMenu({ menuItem("Quality control", icon = icon("lock"), tabName = "") })
    output$diffrm    <- renderMenu({ menuItem("Differential analysis", icon = icon("lock"), tabName = "") })
    output$enrichrm  <- renderMenu({ menuItem("Enrichment analysis", icon = icon("lock"), tabName = "") })
    output$networkrm <- renderMenu({ menuItem("Network analysis", icon = icon("lock"), tabName = "") })
    
    unlock_menus <- function() {
        output$quality <- renderMenu({
            menuItem("Quality control", icon = icon("check-circle"), href = NULL,
                     menuSubItem("PCA", tabName = "PCA", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F),
                     menuSubItem("Heatmap", tabName = "samplecorr", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F)
            )
        })
        output$differential <- renderMenu({
            menuItem("Differential analysis", class = 'btn-10', tabName = "differential", icon = icon("adjust"), href = NULL,
                     menuSubItem("ANOVA", tabName = "anova", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F),
                     menuSubItem("Contrasts", tabName = "contrasts", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F)
            )
        })
        
        removeUI(selector = "#qualityrm")
        removeUI(selector = "#diffrm")
        
        updateSelectInput(session, "contrast1", choices = groups)
        
    }



    # File input ----

    observeEvent(input$infile, {
        inFile <- input$infile
        
        if (is.null(inFile))
            return(NULL)
        
        if(input$sep == 1) {
            separator = ','
        } else if(input$sep == 2) {
            separator = ';'
        } else if (input$sep == 3) {
            separator = '\t'
        }
        
        read_file(inFile$datapath, separator)
        
        
        updateNotifications("Dataset successfully uploaded.","check-circle", "success")
        
        sample_data(data_wide)
        
    })
    
    
    
    
    # Identifiers ----
    
    observeEvent(input$verifyIDs, {
        
        id_check <- qob::update_obsolete(rownames(data_wide))
        
        obsolete <- id_check[["table"]][["Entry"]]
        
        assign("id_check", id_check, envir = .GlobalEnv)
        
        message <- paste(length(obsolete), " IDs are obsolete...")
        updateNotifications(message,"info-circle", "info")

    })
    
    
    # NA frequency plot ----
    
    renderNAfreq <- function(){
        output$nafreq <- renderPlot({
            
            # Which elements are NA?
            allNA <- is.na(data_wide)
            
            # Summary of how many TRUEs there are in each row
            NA_frequency <- table(rowSums(allNA))
            
            graphics::barplot(NA_frequency, xlab = "Na frequency", ylab = "Number of genes")
            
        })
    }
    
    renderNAfreq()
    
    
    observeEvent(input$setcutoff, {
        data_wide <- filter_na(input$missingvalues)
        assign("data_wide", data_wide, envir = .GlobalEnv)
        
        unlock_menus()
        renderNAfreq()
        
    })
    
    
    
    # Differential expression: set contrasts ----
    observeEvent(input$contrast1, {
        
        if (!exists("groups")) {return(NULL)}
        else {
            
            cont1 <- input$contrast1
            
            cont2 <- groups[groups != cont1]
            
            updateSelectInput(session, "contrast2", choices = cont2)
            
        }
        
    })
    
    observeEvent(input$setContrast, {
        output$enrichment <- renderMenu({
            menuItem("Enrichment analysis", icon = icon("flask"), href = NULL,
                     menuSubItem("Background data", tabName = "bgdata", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F),
                     menuSubItem("GO enrichment", tabName = "goenrichment", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F),
                     menuSubItem("Pathway enrichment", tabName = "pathwayenrichment", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F)
            )
        })
        output$network <- renderMenu({
            menuItem("Network analysis", icon = shiny::icon("connectdevelop"),
                     menuSubItem("Network", tabName = "network", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F))
        })
        
        pairing <- input$pairing
        
        contrasts <- paste(input$contrast1,input$contrast2, sep = '-')
        
        contrast <- diff_exp(contrasts, pairing)
        assign('contrast', contrast, envir = .GlobalEnv)
        
        
        removeUI(selector = "#enrichrm")
        removeUI(selector = "#networkrm")
        
    })
    
    # Enrichment ----
    
    observeEvent(input$addbg, {
        
        added_bg <- tissues[,input$Tissue]
        added_bg <- added_bg[!is.na(added_bg)] # Remove empty entries
        
        not_in_bg <- added_bg[!added_bg %in% background_data]
        
        exists_in_input <- added_bg[added_bg %in% input_names]
        
        
        background_data <- c(background_data, not_in_bg)
        
        assign("background_data", background_data, envir = .GlobalEnv)
        
        if(length(exists_in_input) > 0) {
            updateNotifications(paste(length(exists_in_input), " genes already in input dataset."), "info-circle", "info")
        }
        
        if(length(not_in_bg) > 0) {
            updateNotifications(paste("Added ", length(not_in_bg), " genes to background."), "check-circle", "success")
        }
        
    })
    
    output$pathtable <- DT::renderDataTable({
        
        bg <- background_data
        
        enrichment_output <- run_pathway_enrichment('REACTOME', bg)
        
        enriched_paths <- enrichment_output[[1]]
        
        interesting_paths <- enrichment_output[[2]]
        
        assign('interesting_paths', interesting_paths, envir = .GlobalEnv)
        
        df <- DT::datatable(enriched_paths[,1:5],
                            options = list(autoWidth = TRUE,
                                           scrollX=TRUE,
                                           columnDefs = list(
                                               list(width = '420px', targets = c(1))
                                           )))
        df %>% DT::formatSignif('Pvalue', digits = 2)
        
        
    })
    
    # Similarity plot
    output$similarity_plot <- renderPlot(
        {
            run_similarity_plot(interesting_paths)
        })
    
    # Volcano plot
    output$volcano_plot <- renderPlot(
        {
            run_volcano_plot(interesting_paths)
        })
    
    observeEvent(input$generatenetwork, {
        source("bin/networks.R")
        
        output$net <- networkD3::renderForceNetwork({
            
            top <- contrast[contrast$adj.P.Val < input$pvaluecutoff,]
            top <- rownames(top)
            
            merged3 <- merged2[merged2$Accession.x %in% top & merged2$Accession.y %in% top,]
            merged3$combined_score <- merged3$combined_score / 100
            
            merged3 <- merged3[merged3$combined_score > input$interactioncutoff, ]
            
            g <- graph_from_data_frame(merged3, directed = T)
            
            g <- simplify(g, remove.multiple = F, remove.loops = T)
            
            wt <- cluster_walktrap(g)
            members <- membership(wt)
            
            g2 <- igraph_to_networkD3(g, group = members)
            
            g2$nodes$group <- results[results$Gene %in% g2$nodes$name,"Rep.Path"]
            
            d3 <- forceNetwork(
                Links = g2$links,
                Nodes = g2$nodes,
                Source = "source",
                Target = "target",
                Value = "value",
                NodeID = "name",
                Group = "group",
                fontSize = 24,
                fontFamily = "sans-serif",
                opacity = 1,
                zoom = T,
                legend = T,
                charge = -20
                
            )
            
            d3
        })
        
        
    })
}