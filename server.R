library(limma)
library(Biobase)
library(qob)
require(ggplot2)
require(ggrepel)
require(RColorBrewer)
require(dplyr)
require(fitdistrplus)

# Passing arguments from R Shiny to external jQuery scripts unfortunately requires som workarounds
# This function is actually necessary...
fade <- function(fadein) {
    return(fadein)
}

# Server ----
server <- function(session, input, output) {
    
    # Notifications ----
    
    # Initialize as empty list
    # Populate with updateNotifications() function
    notification_list <- list() 
    
    task_list <- list(taskItem("Verify identifiers", value = 0, color = "aqua", href = NULL),
                      taskItem("Upload a dataset", value = 100, color = "aqua", href = NULL)
    )
    
    # Make notification_list global var
    assign("notification_list", notification_list, envir = .GlobalEnv)
    
    # Initialize notification menu
    output$notifMenu <- renderMenu({
        dropdownMenu(type = "notifications", .list = notification_list)
    })
    
    # Initialize message menu
    output$helpMenu <- renderMenu({
        dropdownMenu(type = "tasks", .list = task_list)
    })
    
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
    
    # Clicked notification item
    observeEvent(input$linkClicked, {
        showModal(
            modalDialog(title = "Obsolete IDs",
                        
                        DT::renderDataTable({
                            id_check[["table"]]
                        }),
                        size = "l",
                        footer = list(actionButton("updateAll", "Update all"), modalButton("Dismiss"))
            )
        )
    })
    
    # Sidebar menu items ----
    
    # Render "locked" menu items
    output$qualityrm <- renderMenu({ menuItem("Quality control", icon = icon("lock"), tabName = "") })
    output$diffrm    <- renderMenu({ menuItem("Differential analysis", icon = icon("lock"), tabName = "") })
    output$enrichrm  <- renderMenu({ menuItem("Pathway analysis", icon = icon("lock"), tabName = "") })
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
        
        updateSelectInput(session, "contrast1", choices = groups)
        
    }
    
    
    
    # File input ----
    
    # Main data file input
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
        
        read_file(inFile$datapath, separator, "main")
        
        updateNotifications("Dataset successfully uploaded.","check-circle", "success")
        
        sample_data(data_wide)
        
    })
    
    # Annotation file input
    observeEvent(input$anno_infile, {
        inFile <- input$anno_infile
        
        if (is.null(inFile))
            return(NULL)
        
        if(input$anno_sep == 1) {
            separator = ','
        } else if(input$anno_sep == 2) {
            separator = ';'
        } else if (input$anno_sep == 3) {
            separator = '\t'
        }
        
        read_file(inFile$datapath, separator, "anno")
        
        
        updateNotifications("Annotations successfully uploaded.","check-circle", "success")
        
    })
    
    # Identifiers ----
    
    # Verify IDs
    observeEvent(input$verifyIDs, {
        
        if(input$sourceIDtype == 1) {
        
            session$sendCustomMessage("fadeProcess", fade(12))
            
            id_check <- qob::update_obsolete(rownames(data_wide))
            
            obsolete <- id_check[["table"]][["Entry"]]
            
            session$sendCustomMessage("fadeProcess", fade(0))
            
            assign("id_check", id_check, envir = .GlobalEnv)
            
            input_names <- rownames(data_wide)
            assign("input_names", input_names, envir = .GlobalEnv)
            
            message <- paste(length(obsolete), " IDs are obsolete...")
            updateNotifications(message,"info-circle", "info", clickable = T)
            
        } else {
            message <- paste("Currently only UniProtKB IDs can be verified.")
            updateNotifications(message,"exclamation-triangle", "danger")
        }
        
    })
    
    # Convert IDs
    observeEvent(input$convertIDs, {
        if (input$sourceIDtype == 1) {
            si = 'UniprotAC'
        } else if (input$sourceIDtype == 2) {
            si = 'GeneID'
        } else if (input$sourceIDtype == 3) {
            si = 'Ensembl'
        } else if (input$sourceIDtype == 4) {
            si = 'Gene_Name'
        }
        
        if (input$identifiers == 1) {
            ti = 'GeneID'
        } else if (input$identifiers == 2) {
            ti = 'UniprotAC'
        } else if (input$identifiers == 3) {
            ti = 'Ensembl'
        } else if (input$identifiers == 4) {
            ti = 'Gene_Name'
        }
        
        dw <- qob::mapify(rownames(data_wide), source_id = si, target_id = ti)
        data_wide$mapped <- dw
        data_wide <- data_wide[!duplicated(data_wide$mapped),]
        data_wide <- data_wide[!rownames(data_wide) %in% data_wide$mapped,]
        rownames(data_wide) <- data_wide$mapped
        data_wide <- data_wide[,c(1:ncol(data_wide)-1)]
        
        assign("data_wide", data_wide, envir = .GlobalEnv)
        
        updateNotifications(paste(nrow(data_wide), "genes successfully mapped."),"check-circle", "success")

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
            ggplot(naf, aes(x=Var1, y=Freq)) +
                geom_bar(stat="identity", width=.9, fill="tomato3") + 
                theme(axis.text.x = element_text(angle=65, vjust=0.6))
            
            #graphics::barplot(NA_frequency, xlab = "Na frequency", ylab = "Number of genes")
            
        })
    }
    
    renderNAfreq()
    
    
    observeEvent(input$setcutoff, {
        data_wide <- filter_na(input$missingvalues)
        assign("data_wide", data_wide, envir = .GlobalEnv)
        
        unlock_menus()
        renderNAfreq()
        
    })
    
    
    
    # Distributions ----
    observeEvent(input$generatedistributions, {
        
        data_wide_NAex <- na.exclude(data_origin)
        
        fit.lnorm <- tryCatch( apply(data_wide_NAex, 1, function(x) fitdistrplus::fitdist(as.numeric(x), "lnorm")),
                               error = function(e) print("Can't do that."))
        
        fit.norm <- tryCatch( apply(data_wide_NAex, 1,  function(x) fitdistrplus::fitdist(as.numeric(x), "norm")),
                              error = function(e) print("Can't do that."))
        
        # Render "data type" distribution plots
        output$distributions <- renderPlot({
            
            updateSliderInput(session, inputId = "setdist", max = nrow(na.omit(data_wide)))
            
            d1 <- fit.norm
            d2 <- fit.lnorm
            
            val <- input$setdist
            
            par(mfrow = c(2, 2))
            
            denscomp(list(d1[[val]], d2[[val]]))
            qqcomp(list(d1[[val]], d2[[val]]))
            cdfcomp(list(d1[[val]], d2[[val]]))
            ppcomp(list(d1[[val]], d2[[val]]))
            
        })
        
    }, ignoreInit = TRUE)
    
    # Quality control ----
    
    # Render PCA plots
    
    output$pca2dplot <- renderPlot({
        c <- input$contribs
        e <- input$ellipse
        plotPCA(c,e, '2d')
    })
    
    output$pca3dplot <- plotly::renderPlotly({
        plotPCA(0,0, '3d')
    })
    
    # Render heatmap
    
    output$samplecorrheatmap <- renderPlot({
        cor_mat_raw_logged <- log2(data_origin)
        cor_mat_raw_logged[is.na(cor_mat_raw_logged)] <- 0
        cor_mat_raw_logged <- cor(cor_mat_raw_logged)
        
        #annotation_col = data.frame(tmp)
        pheatmap::pheatmap(
            #   annotation = annotation_col,
            cor_mat_raw_logged,
            legend_breaks = c(min(cor_mat_raw_logged), 1),
            legend_labels = c(0, 1)
        )
    })
    
    # Differential expression: set contrasts ----
    observeEvent(input$contrast1, {
        
        if (!exists("groups")){
            return(NULL)
        } else if (exists("groups") == T){
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
            menuItem("Network analysis", icon = icon("vector-square"),
                     menuSubItem("Interactions", tabName = "network", href = NULL, newtab = TRUE,
                                 icon = shiny::icon("angle-double-right"), selected = F))
        })
        
        pairing <- input$pairing
        
        contrasts <- paste(input$contrast1,input$contrast2, sep = '-')
        
        contrast <- diff_exp(contrasts, pairing)
        assign('contrast', contrast, envir = .GlobalEnv)
        
        
        removeUI(selector = "#enrichrm")
        removeUI(selector = "#networkrm")
        updateNotifications("Model fitted successfully.","check-circle", "success")
    })
    
    # Differential expression: output table
    
    output$diffexptable <- DT::renderDataTable({
        
        df <- DT::datatable(contrast,
                            options = list(autoWidth = TRUE,
                                           scrollX=TRUE))
        df %>% DT::formatSignif('adj.P.Val', digits = 2)
        
        
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
    
    # Enrichment ----
    observeEvent(input$resetbg, {
        background_data <- input_names
        assign("background_data", background_data, envir = .GlobalEnv)
        updateNotifications("Background data has been reset.", "info-circle", "info")
    })
    
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
    
    observeEvent(input$generatepathways, {
        
        session$sendCustomMessage("fadeProcess", fade(70))
        
        if(input$usebackground == 1) {
            bg <- input_names
        } else if (input$usebackground == 2) {
            bg <- background_data
        } else if (input$usebackground == 3) {
            bg <- unique(reactome$UniprotID)
        }
        
        if(input$abstractionlevel == 1) {
            abstraction = 'global'
        } else if(input$abstractionlevel == 2) {
            abstraction = 'lowest'
        }
        
        enrichment_output <- run_pathway_enrichment('REACTOME', bg, abstraction)
        
        enriched_paths <- enrichment_output[[1]]
        interesting_paths <- enrichment_output[[2]]
        
        contrast$Rep.Path.Top <- '* NOT SIGNFICANT'
        contrast$Rep.Path.Name <- '* NOT SIGNFICANT'
        contrast$Rep.Path.Score <- 0
        
        for (i in 1:nrow(interesting_paths)) {
            path_name <- interesting_paths[i,'Pathway_name']
            path_topname <- interesting_paths[i,'Pathway_topname']
            score <- interesting_paths[i,'iScore']
            
            prots <- unlist(interesting_paths[i,'Genes'])
            
            for (p in 1:length(prots)) {
                protein <- prots[p]
                
                if (score > contrast[rownames(contrast) == protein,'Rep.Path.Score']) {
                    contrast[rownames(contrast) == protein,'Rep.Path.Name'] <- path_name
                    contrast[rownames(contrast) == protein,'Rep.Path.Top'] <- path_topname
                    contrast[rownames(contrast) == protein,'Rep.Path.Score'] <- score
                }
            }
        }
        
        contrast$Gene <- rownames(contrast)
        results = dplyr::mutate(contrast, sig=ifelse(contrast$adj.P.Val<0.05, "FDR<0.05", "Not Sig"))
        
        assign('results', results, envir = .GlobalEnv)
        
        session$sendCustomMessage("fadeProcess", fade(0))
        
        output$pathtable <- DT::renderDataTable({
            
            df <- DT::datatable(enriched_paths[,1:6],
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
                mat <- run_similarity_plot(interesting_paths)
                pheatmap::pheatmap(mat, cluster_rows = F, cluster_cols = F, fontsize = 11)
            })
        
        # Volcano plot
        output$volcano_plot <- renderPlot(
            {
                run_volcano_plot(contrast, results)
            })
        
        
        
        output$sankey <- networkD3::renderSankeyNetwork({
            
            P <- run_sankey_diagram(results)
            
            # Plot
            s <- sankeyNetwork(Links = P$links, Nodes = P$nodes, Source = 'source',
                               Target = 'target', Value = 'value', NodeID = 'name',
                               fontSize = 12, fontFamily = 'sans-serif', nodeWidth = 60, sinksRight = F)
            
            s
        })
        
    })
    
    observeEvent(input$generatenetwork, {
        session$sendCustomMessage("fadeProcess", fade(60))
        
        source("bin/networks.R")
        
        output$net <- networkD3::renderForceNetwork({
            
            top <- contrast[contrast$adj.P.Val <= input$pvaluecutoff,]
            top <- top[abs(top$logFC) >= input$fccutoff,]
            
            if (input$direction == 1) {
                top <- top[top$logFC > 0,]
            } else if (input$direction == 2) {
                top <- top[top$logFC < 0,]
            }
            
            top <- rownames(top)

            interactions <- interactions3[(interactions3$protein1 %in% top) & (interactions3$protein2 %in% top), ]
            interactions$combined_score <- interactions$combined_score / 100
            interactions <- interactions[interactions$combined_score > input$interactioncutoff, ]
            
            g <- graph_from_data_frame(interactions, directed = F)
            
            g <- igraph::simplify(g, remove.multiple = F, remove.loops = T)
            
            wt <- igraph::cluster_walktrap(g)
            members <- igraph::membership(wt)
            
            g2 <- networkD3::igraph_to_networkD3(g, group = members)
            
            if(input$pathwaylevel == 1) {
                g2$nodes$group <- sapply(g2$nodes$name, FUN = function(x) results[results$Gene == x, "Rep.Path.Top"])
            } else if (input$pathwaylevel == 2) {
                g2$nodes$group <- sapply(g2$nodes$name, FUN = function(x) results[results$Gene == x, "Rep.Path.Name"])
            }
            
            session$sendCustomMessage("fadeProcess", fade(0))
            
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
            
            d3$x$nodes$hyperlink <- paste0(
                'https://www.uniprot.org/uniprot/',
                g2$nodes$name
            )
            
            d3$x$options$clickAction = 'window.open(d.hyperlink)'
            
            d3
        })
        
        
    })
    
    observeEvent(input$searchclick, {
        
        searchProtein <- function() {
            return(input$seachinput)
        }
        
        session$sendCustomMessage("searchProtein", searchProtein())
    })
    
    
    
}
