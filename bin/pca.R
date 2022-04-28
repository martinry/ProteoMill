
pcaUI <- function(id) {
	ns <- NS(id)
	
	fluidRow(
		column(width = 3,
			   shinydashboard::box(title = "PCA settings",
			   					status = "warning",
			   					width = NULL,
			   					shiny::sliderInput(ns("pcaDims"), "Dimensions", min = 1, max = 10, value = c(1, 2)),
			   					checkboxInput(ns("showPolygons"), "Show polygons", value = T),
			   					checkboxInput(ns("multilevelPCA"), "Multi-level", value = F),
			   					bsTooltip("pcaDims", "Should the plot be displayed in 2D or 3D? Which principal components should be visualized?",
			   							  "right", options = list(container = "body")),
			   					bsTooltip("showPolygons", "Should the area between samples be displayed?",
			   							  "right", options = list(container = "body")),
			   					bsTooltip("multilevelPCA", "Should PCA adjust for paired samples?",
			   							  "right", options = list(container = "body"))
			   )
		),
		column(width = 9,
			   shinydashboard::box(width = NULL,
			   					plotly::plotlyOutput(ns("PCAplots"), width = "95%", height = "500px"),
			   					shinycssloaders::withSpinner(plotOutput(ns("scree"), height = "250px"), type = 5, color = "#e80032dd", id = "ScreeSpinner", size = 1),
			   					bsTooltip("scree", "This plot tells us the number of components that are of relevance to our PCA. Each principal component in the screen plot is less informative than the previous one in terms of how much variation is explained. What we are typically looking for is a steep drop to be able to determine where which principal components are of value to our analysis.",
			   							  "left", options = list(container = "body"))
			   )
		)
	)
	

}

pcaServer <- function(id, pcadata, sampleinfo, palette) {
	
	moduleServer(
		id,
		function(input, output, session) {
			
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
				#pca.data <- maindata$udat@main
				
				if(!is.null(pcadata)) {
					
					pcadata[is.na(pcadata)] <- 0 # "Impute" missing values as 0
					pcadata <- IMIFA::pareto_scale(pcadata)
					
					if(input$multilevelPCA == T) {
						
						ndim <- input$pcaDims[2] - input$pcaDims[1] + 1
						p.pca = mixOmics::pca(X = t(pcadata), ncomp = NULL, multilevel = sampleinfo$samples$replicate, logratio = 'none', scale = F, center = T)
						
						names(p.pca)[2] <- "x"
						
					} else {
						
						p.pca <- prcomp(t(pcadata), center = TRUE, scale. = F)
					}
					p.pca
				}
			})
			
			pca2d <- reactive({
				
				p.pca <- prepare_pca()
				
				si_treatment <- sampleinfo$samples$treatment
				
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
							scale_color_manual(values = palette) +
							theme_light() + theme(legend.title = element_blank())) -> p
				} else {
					suppressWarnings(
						ggplot(mydf, aes(pc1, pc2, colour = as.factor(si_treatment))) +
							geom_point(size = 5, aes(text = rownames(mydf))) +
							xlab(paste0("PC", input$pcaDims[1])) +
							ylab(paste0("PC", input$pcaDims[2])) +
							scale_color_manual(values = palette) +
							theme_light() + theme(legend.title = element_blank())) -> p
				}
			})
			
			output$PCAplots <- renderPlotly({
				
				if(input$pcaDims[2] - input$pcaDims[1] == 1) {
					pcaplot <- pca2d()
					ggplotly(pcaplot, tooltip = "text") %>% layout(legend = list(title=list(text='<b> Treatment </b>'))) %>% 
						config(displaylogo = F,
							   showTips = F,
							   scrollZoom = F,
							   modeBarButtonsToRemove = list(
							   	'sendDataToCloud',
							   	'toImage',
							   	'autoScale2d',
							   	'resetScale2d',
							   	'hoverClosestCartesian',
							   	'hoverCompareCartesian',
							   	'pan2d',
							   	'zoomIn2d',
							   	'zoomOut2d'),
							   modeBarButtonsToAdd = list(
							   	'drawopenpath',
							   	'eraseshape'))
				} else if (input$pcaDims[2] - input$pcaDims[1] == 2) {
					p.pca <- prepare_pca()
					plotly::plot_ly(x = p.pca$x[,input$pcaDims[1]],
									y = p.pca$x[,(input$pcaDims[2] - 1)],
									z = p.pca$x[,input$pcaDims[2]],
									text = rownames(p.pca$x),
									hoverinfo = "text",
									color = sampleinfo$samples$treatment,
									colors = palette
					) %>%
						plotly::add_markers(marker=list(size = 15, line=list(width = 1, color = "black"))) %>%
						plotly::layout(scene = list(xaxis = list(title = paste0("PC", input$pcaDims[1])),
													yaxis = list(title = paste0("PC", (input$pcaDims[2] - 1))),
													zaxis = list(title = paste0("PC", input$pcaDims[2]))),
									   legend = list(title=list(text='<b> Treatment </b>'))) %>% 
						config(displaylogo = F,
							   showTips = F,
							   scrollZoom = T,
							   modeBarButtonsToRemove = list(
							   	'sendDataToCloud',
							   	'toImage',
							   	'autoScale2d',
							   	'resetScale2d',
							   	'hoverClosestCartesian',
							   	'hoverCompareCartesian',
							   	'pan2d',
							   	'zoomIn2d',
							   	'zoomOut2d'),
							   modeBarButtonsToAdd = list(
							   	'drawopenpath',
							   	'eraseshape'))
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

		}
	)
}
