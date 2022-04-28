
tokenUI <- function(id) {
	ns <- NS(id)
	
		column(width = 6,
			   box(width = NULL,
			   	title = "Reproducibility token",
			   	textAreaInput(ns("repToken"), "Token", placeholder = "Your token will appear here."),
			   	actionButton(ns("createToken"), "Create token")
# 			   	bsModalNoClose(id = ns("tokenModal"),
# 			   				   title = "Please note", 
# 			   				   trigger = "noTrigger",
# 			   				   size = "medium",
# 			   				   p("This token is a snapshot of your current settings and databases currently used in ProteoMill.
#                                                   If you make any changes to the experiment after this token has been generated, the experiment may not be reproducible.
#                                                   To reproduce the experiment, an identical version of the input file must be used."),
# 			   				   p("Results may still deviate from the original experiment due to differences in package versions,
#                                                 or issues loading the original settings.")
# 			   	)
			   )
			   
		)
	
	
}

tokenServer <- function(id, maindata, appMeta) {
	
	moduleServer(
		id,
		function(input, output, session) {
			
			generateToken <- reactive({
				
				mylist <- list("date" = as.integer(format(Sys.Date(), "%Y%m%d")),
							   "taxid" = maindata$taxid,
							   "STRINGDB" = getFilemd5sum(getlibfPaths(taxid = maindata$taxid, lib = "STRINGDB")),
							   "REACTOMEDB" = getFilemd5sum(getlibfPaths(taxid = maindata$taxid, lib = "REACTOMEDB")),
							   "ORGDB" = getlibfPaths(taxid = maindata$taxid, lib = "ORGDB"),
							   "inFilemd5" = getFilemd5sum(maindata$fpath),
							   "LogTransformData" = input$LogTransformData,
							   "multilevelPCA" = input$multilevelPCA,
							   "corMethod" = input$corMethod,
							   "missingValues" = input$missingValues,
							   "contrast1" = input$contrast1,
							   "contrast2" = input$contrast2)
				
				mykey <- paste0(encrypt_object(mylist, key = "proteomill"), collapse = " ")
				
				
			})
			
			observeEvent(input$createToken, {
				
				if(!is.null(maindata$udat) & !is.null(maindata$taxid)) {
					
					toggleModal(session = session, "tokenModal")
					
					mykey <- generateToken()
					
					updateTextAreaInput(session = session, inputId = "repToken", value = mykey)
				} else {
					updateNotifications("Upload a dataset first.","exclamation-triangle", "danger")
				}
				
				
			})
			
			observeEvent(input$useInputToken, {
				
				if(is.null(input$inputToken) | nchar(input$inputToken) < 100) {
					createAlert(session, "invalidToken",
								content = "Invalid token. This token was not recognized.", append = T,
								style = "info")
				} else {
					
					token <- input$inputToken
					
					token <- tryCatch({as.raw(as.hexmode(unlist(strsplit(token, " "))));},
									  error = function(err){
									  	return(err)
									  },
									  warning = function(war){
									  	return(war)
									  }, silent=F)
					
					if("error" %in% class(token)) {
						
						createAlert(session, "invalidToken",
									content = "Invalid token. This token was not recognized.", append = T,
									style = "info")
					} else {
						token_decrypted <- decrypt_object(token, "proteomill")
						
						if(class(token_decrypted) == "list") {
							
							createAlert(session, "validToken",
										content = "Experiment settings have been loaded.", append = T,
										style = "success")
							
							orgs <-     c("Homo sapiens | HSA | 9606",
										  "Mus musculus | MMU | 10090",
										  "Rattus norvegicus | RNO | 10116",
										  "Other")
							
							taxid <- token_decrypted$taxid
							
							updateCheckboxInput(session = session, inputId = "LogTransformData", value = token_decrypted$LogTransformData)
							updateSelectInput(session = session, inputId = "organism", selected = orgs[endsWith(orgs, as.character(taxid))])
							
							appMeta$token_decrypted <- token_decrypted
							appMeta$hasToken <- T
							
						} else {
							createAlert(session, "invalidToken", title = "Invalid token",
										content = "This token was not recognized.", append = T,
										style = "info")
						}
					}
				}
				
				
				
				
			})
			
			

		}
	)
}
