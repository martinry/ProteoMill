

contrastUI <- function(id) {
	ns <- NS(id)
	
	fluidRow(
		column(width = 4,
			   box(width = NULL,
			   	title = "Set contrasts", status = "primary", solidHeader = F,
			   	helpText("Which treatments should be compared?"),
			   	selectInput(ns("contrast1"), label = "Treatment 1",
			   				choices = list("Treatment1" = 1,
			   							   "Treatment2" = 2,
			   							   "Treatment3" = 3),
			   				selected = 1),
			   	selectInput(ns("contrast2"), label = "Treatment 2",
			   				choices = list("Treatment1" = 1,
			   							   "Treatment2" = 2,
			   							   "Treatment3" = 3),
			   				selected = 2),
			   	
			   	actionButton(ns("setContrast"), "Select")
			   )
		),
		column(width = 4,
			   box(width = NULL,
			   	selectInput(ns("setDEengine"),
			   				label = "Set engine",
			   				choices = list("Limma version 3.39.1" = 1,
			   							   "DESeq2 version 3.10" = 2),
			   				selected = 1),
			   	bsTooltip("setDEengine", "Limma and DESeq2 are common tools for differential expression. If you are unsure what to choose, a general guideline is that DESeq2 is better suited for RNA-seq data, while limma is more appropriate for MS-proteomics and microarray data.",
			   			  "left", options = list(container = "body")),
			   	radioButtons(ns("diffexppairing"), "Pairing", choices = list("Paired" = 1, "Unpaired" = 2), inline = T, selected = 2),
			   	bsTooltip("diffexppairing", "Should the statistical test adjust for paired samples?",
			   			  "left", options = list(container = "body"))
			   )
		)
	)
}


differentialServer <- function(id, maindata, sampleinfo) {
	
	moduleServer(
		id,
		function(input, output, session) {
			
			
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
				
				if(!is.null(contrast)){
					rcont$contrasts <- contrasts
					rcont$contrast <- contrast
					
					updateNotifications("A DE contrast has been set.","check-circle", "success")
					updateTasks(text = "Set contrast", value = 100, color = "green", i = 0005)
				}
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
					
					if(input$LogTransformData == F) {
						isolate(updateNotifications(paste0("DESeq2 needs raw counts as input."), "exclamation-triangle", "danger"))
						
						return(NULL)
					} else {
						
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
							
							design <- DESeqDataSetFromMatrix(countData  = round(dw),
															 colData    = sinfo,
															 design     = ~ 0 + treatment + replicate)
							
						} else {
							# Unpaired
							design <- DESeqDataSetFromMatrix(countData  = round(dw),
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
					
					# Generate data frame with results from linear model fit
					contrast <- topTable(fit.cont, number = Inf, coef = coeff)
					
					contrast <- data.table::as.data.table(contrast, keep.rownames = T)
					
					contrast <- contrast[order(contrast[, 1])]
					
					maindata$udat@deoutput <- contrast[, 2:ncol(contrast)]
					
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
			
		}
	)
}
