
heatmapUI <- function(id) {
	ns <- NS(id)
	
	fluidRow(
		column(width = 3,
			   box(title = "Heatmap settings",
			   	width = NULL,
			   	status = "warning",
			   	selectInput(ns("corMethod"), "Corr. method", choices = list("Pearson" = "pearson",
			   																"Spearman" = "spearman"),
			   				selected = "pearson"),
			   	shiny::checkboxInput(ns("showGrid"), "Show grid", value = T),
			   	bsTooltip("corMethod", "If method is \\'spearman\\', Spearman\\'s rho statistic is used to estimate a rank-based measure of association. These are more robust and have been recommended if the data do not necessarily come from a bivariate normal distribution.",
			   			  "right", options = list(container = "body"))
			   )),
		
		column(width = 9,
			   box(title = "Sample-sample correlation heatmap", width = NULL,
			   	shinycssloaders::withSpinner(plotlyOutput(ns("heatmap"), height = "600px"), type = 5, color = "#e80032dd", id = "HeatmapSpinner", size = 1)
			   	)
			   )
		)
}

heatmapServer <- function(id, dw) {
	
	moduleServer(
		id,
		function(input, output, session) {

			prepare_heatmap <- reactive({
				
				if(!is.null(dw)){
					dw[is.na(dw)] <- 0 # Impute
					dw <- round(cor(dw, method = input$corMethod), 2)
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
		}
	)
}
