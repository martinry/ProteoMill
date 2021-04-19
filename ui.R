# Load packages ----
require(shiny)
require(shinydashboard)
require(shinyWidgets)
require(shinyjs)
library(shinyBS)
require(visNetwork)
require(DT)
require(plotly)

# Credit to https://stackoverflow.com/questions/50368690/change-backdrop-for-a-bsmodal-in-shiny-app#answer-50570420
# for this solution
bsModalNoClose <-function(...) {
    b = bsModal(...)
    b[[2]]$`data-backdrop` = "static"
    b[[2]]$`data-keyboard` = "false"
    return(b)
}


# Notification menus ----
help <- shinydashboard::dropdownMenuOutput("helpMenu")
notifications <- shinydashboard::dropdownMenuOutput("notifMenu")

header <- dashboardHeader(
    help,
    notifications,
    title = list(span(class = "no-touch", tags$a(href = "https://proteomill.com/", tags$img(id = "mill", class = "normal", src = "mill.png"), "ProteoMill"))),
    tags$li(class = "dropdown",
            id = "notifications-wrapper",
            tags$div(id = 'load-process', style = 'display: none; position: absolute; margin-left: 6px',
                     tags$img(id = 'load-img', src = 'dna.svg', style = 'margin-top: -12px; margin-left: 5px; width: 45px; opacity: .9;'),
                     tags$span(id = "process-counter", 0, style = 'font-size: 14px; font-family: "Courier"; vertical-align: top; padding-left: 3px;')
            ),
            tags$span(class = "loading-menus", tags$text("Loading libraries, please wait...")),
            tags$i(id = "notif-icon"),
            tags$div(class = "ml11",
                     tags$span(class = "text-wrapper2",
                               tags$span(class = "line line1"),
                               tags$span(class = "letters2")
                     )
            )
    )
)

# Sidebar menu ----



sidebar <- dashboardSidebar(
    tags$div(id = 'test', "Analysis", style = "
             letter-spacing: 4px;
             line-height: 42px;
             text-transform: uppercase;
             text-align: center;
             background-color: #425664;
             color: #fff;"),
    sidebarMenu(id = "sidebarmenu",
                menuItem("Dataset options", icon = icon("table"),
                         menuSubItem("Upload data", tabName = "file-input", href = NULL, newtab = TRUE,
                                     icon = shiny::icon("caret-right"), selected = T),
                         menuSubItem("Data summary", tabName = "data-summary", href = NULL, newtab = TRUE,
                                     icon = shiny::icon("caret-right"), selected = F),
                         menuSubItem("Missing  values", tabName = "filters", href = NULL, newtab = TRUE,
                                     icon = shiny::icon("caret-right"), selected = F),
                         startExpanded = T
                         
                ),
                menuItemOutput("qualityrm"),
                menuItemOutput("quality"),
                menuItemOutput('diffrm'),
                menuItemOutput("differential"),
                menuItemOutput('enrichrm'),
                menuItemOutput('enrichment'),
                menuItemOutput('networkrm'),
                menuItemOutput('network'),
                menuItemOutput('interactionsrm'),
                menuItemOutput('interactions')
                
                
    ),
    tags$br(),
    tags$div(id = 'tools', "Additional tools", style = "
             letter-spacing: 3.3px;
             line-height: 42px;
             text-transform: uppercase;
             text-align: center;
             background-color: #425664;
             margin-top: 2px;
             color: #fff;"),
    sidebarMenu(id = "sidebarmenu",
                menuItem("Identifier tools", tabName = "validateIDs", icon = icon("font")),
                menuItem("Goodness-of-fit", tabName = "goodnessOfFit", icon = icon("chart-bar")),
                menuItem("Generate report", tabName = "file-export", icon = icon("file-alt"))
    ),
    tags$br(),
    tags$div(id = 'test', "Information", style = "
             letter-spacing: 3.3px;
             line-height: 42px;
             text-transform: uppercase;
             text-align: center;
             background-color: #425664;
             margin-top: 2px;
             color: #fff;"),
    sidebarMenu(id = "sidebarmenu",
                menuItem("News", tabName = "news", icon = icon("bullhorn"), badgeLabel = "Updates", badgeColor = "blue"),
                menuItem("Documentation", icon = icon("book"),
                         badgeLabel = "New", badgeColor = "red", href = "https://bookdown.org/martin_ryden/proteomill_documentation/"),
                menuItem("Settings", tabName = "settings", icon = icon("sliders-h"))
                
    ),
    useShinyjs()
)


# Main content ----

body <- dashboardBody(
    
    # Import css and js
    tags$head(
        tags$meta(name="google-site-verification", content="JC0Ph8rzlXWiAL6lWXnusIUEOhJSqf8u2yVzK5g2P04"),
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
        tags$link(rel = "stylesheet", type = "text/css", href = "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css"),
        tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css?family=Quicksand"),
        tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css?family=Patrick+Hand"),
        tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css?family=Spectral+SC"),
        tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css?family=Poiret+One"),
        tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css?family=Open+Sans"),
        tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/animejs/2.0.2/anime.min.js") # To do: keep local copy
        
    ),
    
    # Clickable notifications
    
    # Ugly but unavoidable hack. The onInputChange function only runs on unique links,
    # so we generate a random number to be able to click a notification item even after
    # closing the modal dialog.
    
    tags$script(HTML("function clickFunction(link){
                     var number = Math.random();
                     link = number + link;
                     Shiny.onInputChange('linkClicked',link);
                     }")),
    
    
    tags$script(src = "custom.js"),
    tags$script(src = "animate-notifs.js"),
    
    
    # Intro animation
    
    tags$div(class = "overlay",
             tags$h1(class = "ml9",
                     tags$span(class = "text-wrapper",
                               tags$span(class = "letters", "ProteoMill"))),
             tags$div(class = 'begindiv',
                      actionLink('removeBanner', label = 'LAUNCH')
             ),
             tags$div(style = "margin-top: 55px; opacity: .85;",
                      tags$video(playsinline = "playsinline", loop = "true", autoplay = "autoplay", muted = "muted", width="60%", height="auto",
                                 tags$source(src = "media1.mp4", type="video/mp4")))
    ),
    tags$script(src = "custom.js"),
    tags$script(src = "animate-notifs.js"),
    
    tabItems(
        # File input ----
        tabItem(tabName = "file-input",
                fluidRow(
                    column(width = 5,
                           tabBox(width = NULL,
                                  height = 250,
                                  tabPanel(title = "Data import wizard",
                                           p(
                                               helpText(
                                                   'Welcome! Get started by using the ',
                                                   tags$strong('Data import wizard'),
                                                   ' or download one our demo datasets to learn more about accepted file formats.'
                                               )
                                           ), 
                                           actionButton("ImportWizard", "Data import wizard...", icon = icon("magic")),
                                           
                                           bsModalNoClose(id ="ImportModal1",
                                                   title = "Data import wizard", 
                                                   trigger = "ImportWizard",
                                                   size = "large",
                                                   
                                                   bsAlert("unsupportedOrganism"),
                                                   
                                                   fluidRow(
                                                       column(width = 6,
                                                              box(width = NULL,
                                                                  title = "Sample information",
                                                                  
                                                                  
                                                                  # Currently only three organisms support UniprotID in ensembldb package
                                                                  # Will write to developer
                                                                  
                                                                  selectInput("organism",
                                                                              "Select organism",
                                                                              list("Homo sapiens | HSA | 9606",
                                                                                   "Mus musculus | MMU | 10090",
                                                                                   "Rattus norvegicus | RNO | 10116",
                                                                                   "Other")),
                                                                  
                                                                  bsTooltip("organism", "In order to use correct annotation databases, we need to know the species your data is derived from.",
                                                                            "right", options = list(container = "body"))
                                                                  
                                                                  # selectInput("DataType", "Type of data",
                                                                  #             list("Proteins",
                                                                  #                  "Peptides + proteins")),
                                                                  # bsTooltip("DataType", "ProteoMill supports peptide-level and protein-level data.",
                                                                  #           "right", options = list(container = "body"))
                                                              )
                                                       ),
                                                       
                                                       column(width = 6,
                                                              box(width = NULL,
                                                                  title = "Data configuration",
                                                                  
                                                                  checkboxInput("LogTransformData", "Log₂-transform data", value = T),
                                                                  bsTooltip("LogTransformData", "Throughout the analysis, log-transformed data will be used. If your data is already log₂-transformed, uncheck this box.",
                                                                            "right", options = list(container = "body"))
                                                                  
                                                              ))
                                                       
                                                   ),
                                                   
                                                   
                                                   tags$head(tags$style("#Modal1Spinner {display:none}")),
                                                   shinycssloaders::withSpinner(
                                                       textOutput("plot"), type = 6, color = "#e80032dd", id = "Modal1Spinner", size = 0.4, proxy.height = "50px"
                                                   ),
                                                   
                                                   tags$p(tags$hr()),
                                                   modalButton("Cancel"),
                                                   actionButton("EndStep1", "Next"),
                                                   tags$head(tags$style("#ImportModal1 .modal-footer{ display:none}"))
                                                   
                                           ),
                                           
                                           bsModalNoClose(id = "ImportModal2",
                                                   title = "Data import wizard",
                                                   trigger = "BeginStep2",
                                                   size = "large",
                                                   
                                                   tags$head(tags$style("#ImportModal2 .modal-footer{ display:none}")),
                                                   tags$head(tags$style("#showHideBox1,#showHideBox2 {display:none}")),
                                                   
                                                   bsAlert("selectAFile"),
                                                   bsAlert("noneConverted"),
                                                   bsAlert("fileHasError"),
                                                   bsAlert("fileHasWarning"),
                                                   bsAlert("fileFormattingError"),
                                                   bsAlert("checkRequirements"),
                                                   
                                                   
                                                   
                                                   
                                                   fluidRow(
                                                       column(width = 12,
                                                              tags$p(tags$strong("A minimal example."), actionLink(inputId = "ShowHide1", "Show/Hide")),
                                                              box(id = "showHideBox1",
                                                                  width = NULL,
                                                                  
                                                                  tags$p(img(src = "/img/yellow.png", width = "10px"), tags$strong("Protein IDs."), "In the first column, row two and onwards contains protein IDs (see accepted ID types below)."),
                                                                  
                                                                  tags$p(img(src = "/img/blue.png", width = "10px"), tags$strong("Sample names."), "In the first row, column two and onward contains the sample names, with treatments and replicates separated by an underscore character. Sample column order does not matter."),
                                                                  
                                                                  tags$p(img(src = "/img/green.png", width = "10px"), tags$strong("Quantitative values."), 'Only numeric values greater than or equal to 0, or NA. Period (.) should be used for decimal character.'),
                                                                  
                                                                  img(src = "/img/minimal_example.png", width = "75%")
                                                              ),
                                                              tags$p(tags$strong("Dataset checklist."), actionLink(inputId = "ShowHide2", "Show/Hide")),
                                                              box(id = "showHideBox2",
                                                                  width = NULL,
                                                                  
                                                                  tags$ul(
                                                                      tags$li("Experiment must contain at least two groups/treatments"),
                                                                      tags$li("Each treatment must contain at least two replicates"),
                                                                      tags$li("File must have a header row"),
                                                                      tags$li("First column should contain protein IDs"),
                                                                      tags$li("Column two and onwards should contain sample names"),
                                                                      tags$li("Sample names should be in the form treatment_replicate")
                                                                  )
                                                                  
                                                              )
                                                       )
                                                   ),
                                                   
                                                   fluidRow(
                                                       column(width = 4,
                                                              box(width = NULL,
                                                                  fileInput(inputId = "file1",
                                                                            label = "Upload a dataset",
                                                                            multiple = F,
                                                                            accept = c(
                                                                                "text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".csv")),
                                                                  
                                                                  bsTooltip("file1", "Upload a dataset of filetyp .csv, .tsv, or .txt.",
                                                                            "right", options = list(container = "bsModal"))
                                                                  
                                                              )
                                                       ),
                                                       column(width = 8,
                                                              tags$div(id = "dataDetailsWrapper", style = "display: none;",
                                                                       box(width = NULL,
                                                                           height = 119,
                                                                           title = "Dataset details",
                                                                           htmlOutput("dataDetails")
                                                                       )
                                                              )
                                                       )
                                                   ),
                                                   
                                                   tags$head(tags$style(
                                                       HTML("#DataTables_Table_0_length {visibility:hidden}"),
                                                       HTML("#DataTables_Table_0_filter {visibility:hidden}"),
                                                       HTML("#DataTables_Table_0_info {visibility:hidden}"),
                                                       HTML("#DataTables_Table_0_paginate {visibility:hidden}")
                                                       
                                                       
                                                   )),
                                                   # div(id = "previewDTInfo", tags$h4("Data preview"),
                                                   #   style = "display:none"),
                                                   div(id = "previewDTInfo", style = "display:none;",
                                                       box(
                                                           width = NULL,
                                                           title = "Data preview",
                                                           tags$div(DTOutput("previewDT"), style = "font-size: 58%; margin-top: -35px; max-width: 100%")
                                                       )
                                                   ),
                                                   
                                                   
                                                   tags$head(tags$style("#Modal2Spinner {display:none}")),
                                                   shinycssloaders::withSpinner(
                                                       textOutput("plot2"), type = 6, color = "#e80032dd", id = "Modal2Spinner", size = 0.4, proxy.height = "50px"
                                                   ),
                                                   
                                                   
                                                   tags$p(tags$hr()),
                                                   modalButton("Cancel"),
                                                   actionButton("EndStep2", "Next")
                                                   
                                           ),
                                           
                                           bsModalNoClose(id = "ImportModal3",
                                                   title = "Data import wizard",
                                                   trigger = "BeginStep3",
                                                   size = "large",
                                                   
                                                   DT::DTOutput("SamplePreview"),
                                                   
                                                   tags$head(tags$style("#ImportModal3 .modal-footer{ display:none}")),
                                                   
                                                   tags$p(tags$hr()),
                                                   modalButton("Cancel"),
                                                   actionButton("EndStep3", "Finish")
                                                   
                                           )
                                           
                                           
                                  ),
                                  tabPanel(
                                      title = 'Demo data',
                                      p("Explore ProteoMill with a demo datasets. You can also download a dataset to learn how to properly format your own datasets for use in ProteoMill."),
                                      selectInput("selectDemoData", label = "Select a dataset",
                                                  list(`Proteomics` = 
                                                           list("Sample dataset 1" = 1))),
                                      actionButton("useDemoData", label = "Use demo data"),
                                      downloadButton('downloadDemo',"Download")
                                  )
                           )
                    )
                )
        ),
        
        # Data summary ----
        tabItem(tabName = "data-summary",
                fluidRow(
                    column(width = 5,
                           box(title = "Mapped IDs",
                               width = NULL,
                               tableOutput("identifierinfo"),
                               bsTooltip("identifierinfo", "The number of protein IDs in your dataset that matches the database. I.e. - for how many proteins we can obtain annotation data. Proteins for which no match is found will still be included in differential expression analysis.",
                                         "right", options = list(container = "body"))
                           ),
                           box(title = "Include/exclude samples",
                               width = NULL,
                               pickerInput(
                                   inputId = "includesamples", 
                                   label = "Deselect a sample to exclude it", 
                                   choices = LETTERS, 
                                   options = list(
                                       `actions-box` = TRUE, 
                                       size = 10,
                                       `selected-text-format` = "count > 3"
                                   ), 
                                   multiple = TRUE
                               ),
                               actionButton("confirmexclude", label = "Confirm selection"))),
                    column(width = 7,
                           box(width = NULL,
                               infoBoxOutput("datainfoBox"),
                               infoBoxOutput("sampleinfoBox"),
                               infoBoxOutput("conditioninfoBox")
                           ),
                           box(title = "Expression levels",
                               width = NULL,
                               shinycssloaders::withSpinner(plotOutput("violinplot"), type = 5, color = "#e80032dd", id = "ViolinSpinner", size = 1),
                               radioButtons("violintype", "", choices = list("By condition" = 1, "By sample" = 2), inline = T)
                           ))
                )),
        
        # Filters: NA cutoff
        
        tabItem(tabName = "filters",
                fluidRow(
                    column(width = 4,
                           box(
                               width = NULL,
                               title = "Missing values cutoff", status = "primary", solidHeader = F,
                               helpText("Set maximum number allowed missing values for each treatment"),
                               br(),
                               helpText("To proceed without removing any proteins with missing values, set the cutoff to a value higher than the number of samples / treatments"),
                               br(),
                               numericInput("missingvalues", label = "",
                                            min = 0, max = 9999, value = 1), # max = number of samples / conditions
                               actionButton("setcutoff", "Set cutoff")
                           )
                    ),
                    column(width = 8,
                           box(
                               width = NULL,
                               status = "warning",
                               solidHeader = F,
                               plotOutput("nafreq")
                           )
                    )
                )
                
        ),
        
        
        # Identifiers
        
        tabItem(tabName = "validateIDs",
                fluidRow(
                    column(width = 6,
                           box(width = NULL,
                               title = "Validate IDs", status = "primary", solidHeader = F,
                               p(helpText(tags$strong('Please note: '), 'Depending on the number of IDs, this process may take a long time to run.')),
                               actionButton("listCandidates", label = "List outdated IDs"),
                               p(),
                               hr(),
                               p(),
                               DT::DTOutput("obsolete"))
                    ),
                    column(width = 6,
                           box(width = NULL,
                               title = "Convert IDs",
                               helpText("Insert a list of identifiers. Separators are automatically detected."),
                               textAreaInput("idstoconvert", "", height = "150px"),
                               actionButton("convertids", label = "Convert IDs"),
                               p(),
                               hr(),
                               p(),
                               DT::DTOutput("convertedids"))
                    )
                )
        ),
        
        
        # Data type: distributions
        
        tabItem(tabName = "goodnessOfFit",
                box(width = NULL,
                    title = "Assess goodness-of-fit", status = "primary", solidHeader = F,
                    selectInput(inputId = "fit", label = "Fit to distribution",
                                choices = list("Normal/Log-Normal" = 1),
                                selected = 1
                    ),
                    actionButton("generatedistributions", label = "Render distributions"),
                    tags$p(),
                    plotOutput("distributions"),
                    sliderInput(
                        "setdist",
                        label = "Select gene",
                        min = 1,
                        max = 1,
                        value = 1,
                        step = 1,
                        animate = animationOptions(interval = 600)
                    )
                )
        ),
        
        # BLAST
        
        tabItem(tabName = "blast",
                fluidRow(
                    column(width = 3,
                           box(width = NULL,
                               helpText("Find closest human homologue using BLAST.")
                           )
                    ),
                    column(width = 3,
                           box(width = NULL,
                               
                           )
                    )
                )
                
        ),
        
        # Generate report
        
        tabItem(tabName = "file-export",
                box(
                    title = "Generate report", status = "primary", solidHeader = F,
                    p(helpText("Exports generated content to an Rmarkdown document")),
                    pickerInput(
                        inputId = "p2",
                        label = "",
                        choices = list("Dataset options" = c("Data summary", "Missing values"),
                                       "Inspect data" = c("PCA 2D", "PCA 3D", "Heatmap"),
                                       "Differential analysis" = "Differential expression",
                                       "Enrichment analysis" = c("Pathway enrichment", "Volcano plot", "Sankey diagram")),
                        multiple = TRUE,
                        options = list(
                            `actions-box` = TRUE,
                            `deselect-all-text` = "Deselect all",
                            `select-all-text` = "Select all"
                        )),
                    
                    downloadButton("downloadReport", "Generate report")
                    
                )
        ),
        
        # Quality control: PCA: 2D, 3D
        
        tabItem(tabName = "PCA",
                fluidRow(
                    column(width = 3,
                           shinydashboard::box(title = "PCA settings",
                                               status = "warning",
                                               width = NULL,
                                               shiny::sliderInput("pcaDims", "Dimensions", min = 1, max = 10, value = c(1, 2)),
                                               checkboxInput("showPolygons", "Show polygons", value = T),
                                               checkboxInput("multilevelPCA", "Multi-level", value = F),
                                               bsTooltip("pcaDims", "Should the plot be displayed in 2D or 3D? Which principal components should be visualized?",
                                                         "right", options = list(container = "body")),
                                               bsTooltip("showPolygons", "Should the area between samples be displayed?",
                                                         "right", options = list(container = "body")),
                                               bsTooltip("multilevelPCA", "Should PCA adjust for paired samples?",
                                                         "right", options = list(container = "body"))
                           )
                    ),
                    column(width = 9,
                           shinydashboard::box(
                               #shinycssloaders::withSpinner(, type = 5, color = "#e80032dd", id = "PCASpinner", size = 1),
                               plotly::plotlyOutput("PCAplots", width = "95%", height = "500px"),
                               shinycssloaders::withSpinner(plotOutput("scree", height = "250px"), type = 5, color = "#e80032dd", id = "ScreeSpinner", size = 1),
                               width = NULL)
                    ))),
        
        tabItem(tabName = "samplecorr",
                fluidRow(
                    column(width = 3,
                           box(title = "Heatmap settings",
                               width = NULL,
                               status = "warning",
                               selectInput("corMethod", "Corr. method", choices = list("Pearson" = "pearson",
                                                                                       "Spearman" = "spearman",
                                                                                       "Kendall" = "kendall"),
                                           selected = "pearson"),
                               # bsTooltip("corMethod", "if method is 'kendall' or 'spearman', Kendall's tau or Spearman's rho statistic is used to estimate a rank-based measure of association. These are more robust and have been recommended if the data do not necessarily come from a bivariate normal distribution.",
                               #           placement = "bottom", trigger = "hover", options = list(container = "body")),
                               shiny::checkboxInput("showGrid", "Show grid", value = T),
                               bsTooltip("corMethod", "If method is \\'kendall\\' or \\'spearman\\', Kendall\\'s tau or Spearman\\'s rho statistic is used to estimate a rank-based measure of association. These are more robust and have been recommended if the data do not necessarily come from a bivariate normal distribution.",
                                         "right", options = list(container = "body"))
                           )),
                    
                    column(width = 9,
                           box(title = "Sample-sample correlation heatmap", width = NULL,
                               shinycssloaders::withSpinner(plotlyOutput("heatmap", height = "600px"), type = 5, color = "#e80032dd", id = "HeatmapSpinner", size = 1),
                               )
                    ))
                
        ),
        
        # Differential expression analysis
        
        tabItem(tabName = "contrasts",
                fluidRow(
                    column(width = 4,
                           box(width = NULL,
                               title = "Set contrasts", status = "primary", solidHeader = F,
                               helpText("Which treatments should be compared?"),
                               selectInput("contrast1", label = "Condition 1",
                                           choices = list("Condition1" = 1,
                                                          "Condition2" = 2,
                                                          "Condition3" = 3),
                                           selected = 1),
                               selectInput("contrast2", label = "Condition 2",
                                           choices = list("Condition1" = 1,
                                                          "Condition2" = 2,
                                                          "Condition3" = 3),
                                           selected = 2),
                               
                               actionButton("setContrast", "Select")
                           )
                    ),
                    column(width = 4,
                           box(width = NULL,
                               selectInput("setDEengine",
                                           label = "Set engine",
                                           choices = list("Limma version 3.39.1" = 1,
                                                          "DESeq2 version 3.10" = 2),
                                           selected = 1),
                               bsTooltip("setDEengine", "Limma and DESeq2 are common tools for differential expression. If you are unsure what to choose, a general guideline is that DESeq2 is better suited for RNA-seq data, while limma is more appropriate for MS-proteomics and microarray data.",
                                         "right", options = list(container = "body")),
                               radioButtons("diffexppairing", "Pairing", choices = list("Paired" = 1, "Unpaired" = 2), inline = T, selected = 2))
                           )
                )
                
        ),
        tabItem(tabName = "diffexpoutput",
                plotOutput("contrasttable", width = "800px", height = "1600px")),
        tabItem(tabName = "differentialexpression",
                fluidRow(
                    column(width = 6,
                           box(width = NULL,
                               helpText("Display only proteins that have:"),
                               br(),
                               numericInput("diffexp_limit_fc", "Abs. log2 fold-change greater than or equal to", min = 0, max = 50, value = 0, step = .5),
                               numericInput("diffexp_limit_pval", "Adj. P-value less than", min = 0, max = 1, value = 1, step = .1),
                               bsTooltip("diffexp_limit_fc", "Should PCA adjust for paired samples?",
                                         "right", options = list(container = "body")),
                               bsTooltip("diffexp_limit_pval", "Should PCA adjust for paired samples?",
                                         "right", options = list(container = "body"))
                           ),
                           tabBox(width = NULL,
                                  tabPanel("Summary statistics", tableOutput("diffexptable_summary")),
                                  tabPanel("Up-regulated proteins", DT::dataTableOutput("diffexptable_up")),
                                  tabPanel("Down-regulated proteins", DT::dataTableOutput("diffexptable_down")),
                                  bsTooltip("diffexptable_summary", "Should PCA adjust for paired samples?",
                                            "right", options = list(container = "body")),
                                  bsTooltip("diffexptable_up", "Should PCA adjust for paired samples?",
                                            "right", options = list(container = "body")),
                                  bsTooltip("diffexptable_down", "Should PCA adjust for paired samples?",
                                            "right", options = list(container = "body"))
                                  
                                  )),
                    column(width = 6,
                           box(width = NULL,
                               shinycssloaders::withSpinner(plotlyOutput("dea_volcano", height = "602px"),
                                                            type = 5, color = "#e80032dd", id = "VolcanoSpinner1", size = 1),
                               bsTooltip("dea_volcano", "Should PCA adjust for paired samples?",
                                         "bottom", options = list(container = "body"))
                           ))
                )
                
        ),
        
        # Pathway enrichment: Table, Similarity matrix, Volcano plot
        
        tabItem(tabName = "pathwayenrichment",
                fluidRow(
                    column(width = 3,
                           box(title = "Pathway settings", status = "warning", width = NULL,
                               selectInput(inputId = "pathdb", label = "Pathway database",
                                           choices = list("REACTOME" = 'REACTOME'),
                                           selected = 'REACTOME'),
                               numericInput("min_fc", "Min log2 fold change", value = 0, min = 0, max = 50, step = .5),
                               numericInput("min_pval", "Max adj. P-value", value = 0.05, min = 0, max = 1, step = .01),
                               htmlOutput("number_of_genes")
                           ),
                           
                           box(title = "Selected pathway info", width = NULL,
                               htmlOutput("selected_pathway")
                           )
                    ),
                    column(width = 9,
                           tabBox(width = NULL,
                                  tabPanel("Enrichment of up-regulated proteins", DT::dataTableOutput("upregulated_pathways_table")),
                                  tabPanel("Enrichment of down-regulated proteins", DT::dataTableOutput("downregulated_pathways_table")))
                    )
                )
                ,
                
        ),
        tabItem(tabName = "pathwayvisualization",
                fluidRow(
                    column(width = 3,
                           box(width = NULL,
                               radioButtons("volcanoAnnotation", "Annotation",
                                            choices = list("Differential expression groups" = "group",
                                                           "Enriched pathway (top-level)" = "TopReactomeName",
                                                           "Enriched pathway" = "Pathway_name"
                                                           ),
                                            selected = "Pathway_name")
                           ),
                    ),
                    column(width = 9,
                           tabBox(width = NULL,
                                  tabPanel("Volcano plot",
                                           shinycssloaders::withSpinner(plotlyOutput("volcano_plot", height = 750), type = 5, color = "#e80032dd", id = "volcanoSpinner2", size = 1)
                                           ),
                                  tabPanel("Sankey diagram",
                                           shinycssloaders::withSpinner(networkD3::sankeyNetworkOutput("sankey", height = 750), type = 5, color = "#e80032dd", id = "sankeySpinner", size = 1))
                           )
                    )
                    
                )),
        tabItem(tabName = "interactions",
                fluidRow(
                    box(
                        width = 6,
                        height = 520,#464,
                        shinycssloaders::withSpinner(visNetworkOutput("interaction_network"), type = 5, color = "#e80032dd", id = "interactionSpinner", size = 1)
                        
                    ),
                    box(
                        width = 6,
                        height = 520,
                        shinycssloaders::withSpinner(plotlyOutput("volcano_plot2", width = "auto"), type = 5, color = "#e80032dd", id = "volcanoSpinner3", size = 1)
                    )
                ),
                
                fluidRow(
                    column(width = 3,
                           tabBox(width = NULL, height = 255,
                                  tabPanel(
                                      title = "Network settings",
                                      
                                      radioButtons(
                                          "network_regulation",
                                          label = "Subset by up- or down-regulation",
                                          choices = list("Up" = 1, "Down" = 2, "Both" = 3),
                                          selected = 3, inline = T),
                                      
                                      numericInput(
                                          "fccutoff",
                                          label = "Minimum abs. log2FC",
                                          min = 0, max = 100, value = 0, step = 0.1
                                      ),
                                      
                                      numericInput(
                                          "interactionConfidence",
                                          label = "Interaction confidence",
                                          min = 0, max = 10, value = 7, step = 0.1
                                      )
                                  ),
                                  
                                  tabPanel(
                                      title = "Selection criteria",
                                      
                                      radioButtons("interaction_behaviour", label = "Selection subset",
                                                   choices = list("Strict" = 1,
                                                                  "Extended" = 2),
                                                   inline = T,
                                                   selected = 1),
                                      textAreaInput("network_proteins", label = "Select proteins"),
                                      actionButton("highlight_nodes", label = "Highlight")
                                  )
                           )
                           
                    ),
                    
                    column(width = 3,
                           box(height = 255,
                               title = "Protein info", width = NULL,
                               uiOutput("clicked_node"),
                               uiOutput("hovered_node", style = "height: 184px; overflow-y: scroll;")
                           )
                    ),
                    
                    column(width = 3,
                           box(width = NULL,
                               title = "Volcano annotation",
                               radioButtons("volcanoAnnotation2", "Annotation",
                                            choices = list("Differential expression groups" = "group",
                                                           "Enriched pathway (top-level)" = "TopReactomeName",
                                                           "Enriched pathway" = "Pathway_name"
                                            ), selected = "TopReactomeName")
                               # title = "Protein info", width = NULL,
                               # uiOutput("clicked_node"),
                               # uiOutput("hovered_node")
                           )
                    )
                )
        ),
        
        
        tabItem(tabName = "settings",
                fluidRow(
                    column(width = 5,
                           box(title = "Settings", width = NULL,
                               
                               h6("Display options"),
                               hr(),
                               radioButtons(inputId = "toggleToolTip", label = "Display tool-tip", choices = list("Yes", "No"), selected = "Yes", inline = T),
                               bsTooltip("toggleToolTip", "This is a tooltip and it\\'s currently visible.",
                                         "right", options = list(container = "body")),
                               br(),
                               
                               h6("Identifier type"),
                               hr(),
                               selectInput("displayIdentifier",
                                           label = "Display ID labels as",
                                           choices = list(
                                               "UniProtKB" = 2,
                                               "Entrez" = 3,
                                               "Gene Symbol" = 4,
                                               "Ensembl Gene ID" = 5,
                                               "Ensembl Protein ID" = 6
                                           ),
                                           selected = 1),
                               bsTooltip("displayIdentifier", "This option determines how identifiers are displayed (regardless of the identifier type used in your dataset). It can be changed at any time.",
                                         "right", options = list(container = "body")),
                           )
                    )
                )
        )
    ),
    div(class = "footer_wrapper",
        div(class = "sticky_footer",
            div(span(shiny::icon("github", class = "padded-icons"),
                     shiny::icon("twitter", class = "padded-icons"),
                     shiny::icon("linkedin", class = "padded-icons"),
                     shiny::icon("youtube", class = "padded-icons"),
                     shiny::icon("at", class = "padded-icons"),
                     style = "font-size: 14pt; line-height: 50px;")),
            div(class = "gradient no-touch", span("© 2021 · ProteoMill | Martin Rydén", style = "font-size: 12pt; line-height: 50px; padding-left: 25px;"))))
)

# Load dashboard page ----

dashboardPage(title = "ProteoMill | Differential expression pathway and network analysis tool", skin = 'black', header, sidebar, body)

