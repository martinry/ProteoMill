# Load packages ----
require(shiny)
require(shinydashboard)
require(shinyWidgets)
require(visNetwork)
require(DT)
require(knitr)
require(limma)
require(Biobase)
require(ggplot2)
require(ggrepel)
require(RColorBrewer)
require(dplyr)
require(plotly)
require(data.table)
require(AnnotationDbi)
require(EnsDb.Hsapiens.v86)
require(networkD3)
require(XML)
require(mixOmics)
require(stringr)
require(factoextra)
require(pheatmap)
require(rmarkdown)
require(fitdistrplus)
require(igraph)
require(R.utils)
require(umap)
    



# Notification menus ----

help <- shinydashboard::dropdownMenuOutput("helpMenu")
notifications <- shinydashboard::dropdownMenuOutput("notifMenu")

header <- dashboardHeader(help,
                          notifications,
                          title = list(tags$img(id = "mill", class = "normal", src = "mill.png"), "PROTEOMILL"),
                          tags$li(class = "dropdown",
                                  id = "notifications-wrapper",
                                  tags$div(id = 'load-process', style = 'display: none; position: absolute; margin-left: 6px',
                                    tags$img(id = 'load-img', src = 'dna.svg', style = 'margin-top: -12px; margin-left: 5px; width: 45px; opacity: .9;'),
                                    tags$span(id = "process-counter", 0, style = 'font-size: 14px; font-family: "Courier"; vertical-align: top; padding-left: 3px;')
                                  ),
                          tags$i(id = "notif-icon"),
                          tags$div(class = "ml11",
                                   tags$span(class = "text-wrapper2",
                                             tags$span(class = "line line1"),
                                             tags$span(class = "letters2")
                                             )
                                   )
                                 ),
                          tags$li(class = "dropdown", id = "test-button")
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
                             icon = shiny::icon("angle-double-right"), selected = T),
                 menuSubItem("Data summary", tabName = "data-summary", href = NULL, newtab = TRUE,
                             icon = shiny::icon("angle-double-right"), selected = F),
                 menuSubItem("Missing  values", tabName = "filters", href = NULL, newtab = TRUE,
                             icon = shiny::icon("angle-double-right"), selected = F),
                 startExpanded = T
                 
                 # menuSubItem("Data type", tabName = "datatype", href = NULL, newtab = TRUE,
                 #             icon = shiny::icon("angle-double-right"), selected = F)
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
        menuItem("Validate IDs", tabName = "validateIDs", icon = icon("check-square")),
        menuItem("Goodness-of-fit", tabName = "goodnessOfFit", icon = icon("chart-bar")),
        menuItem("Protein structures", tabName = "structures", icon = icon("fingerprint")),
        menuItem("Generate report", tabName = "file-export", icon = icon("file-download"))
    ),
    tags$br(),
    
    tags$div(id = 'test', "Documentation", style = "
             letter-spacing: 3.3px;
             line-height: 42px;
             text-transform: uppercase;
             text-align: center;
             background-color: #425664;
             margin-top: 2px;
             color: #fff;"),
    sidebarMenu(id = "sidebarmenu",
        menuItem("News", tabName = "news", icon = icon("book")),
        menuItem("Contact", tabName = "contact", icon = icon("at")),
        #menuItem("About", tabName = "about", icon = icon("book")),
        menuItem("Settings", tabName = "settings", icon = icon("sliders-h"))
        
    )
)


# Main content ----

body <- dashboardBody(
    
    # Import css and js
    tags$head(
        tags$meta(name="google-site-verification", content="JC0Ph8rzlXWiAL6lWXnusIUEOhJSqf8u2yVzK5g2P04"),
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
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
                        tabPanel(
                            title = "Dataset",
                            
                            p(helpText('Welcome! Click on',  tags$strong('Upload a dataset'), 'in the task menu ', shiny::icon("tasks"), ' or download one our demo datasets to learn about accepted file formats.')),
                            
                            selectInput("dataSep", label = 'Separator',
                                         choices = list("Auto detect" = 1, "Comma" = 2, "Semicolon" = 3, "Tab" = 4),
                                         selected = 1),
                            selectInput("dataIdentiferType", "Identifier type", choices = list("Auto detect" = 1, "UniProtKB" = 2, "Entrez" = 3, "Gene Symbol" = 4)),
                            fileInput("infile", "Select a file",
                                      accept = c(
                                          "text/csv",
                                          "text/comma-separated-values,text/plain",
                                          ".csv")),
                            helpText('Accepted filetypes are csv, tsv and txt.')
                        ),
                        tabPanel(
                            title = 'Annotation data',
                            selectInput("annoSep", label = 'Separator',
                                        choices = list("Auto detect" = 1, "Comma" = 2, "Semicolon" = 3, "Tab" = 4),
                                        selected = 1),
                            fileInput("anno_infile", "Select a file",
                                      accept = c(
                                          "text/csv",
                                          "text/comma-separated-values,text/plain",
                                          ".csv")),
                            helpText('Accepted filetypes are csv, tsv and txt.')
                        ),
                        tabPanel(
                            title = 'Demo data',
                            selectInput("selectDemoData", label = "Select a dataset",
                                        list(`Proteomics` = 
                                                 list("Sample dataset 1" = 1,
                                                      "Sample dataset 2" = 2,
                                                      "Sample dataset 3" = 3,
                                                      "Sample dataset 4" = 4),

                                             `Transcriptomics` = 
                                                 list("Sample dataset 5" = 5,
                                                      "Sample dataset 6" = 6,
                                                      "Sample dataset 7" = 7,
                                                      "Sample dataset 8" = 8))),
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
                               tableOutput("identifierinfo")
                           ),
                           box(title = "Include/exclude samples",
                               width = NULL,
                               #helpText("Deselect samples that you want to exclude."),
                               #h6(actionLink("samplesselectall","Select all"), strong(" | "), actionLink("samplesdeselectall","Deselect all")),
                               #checkboxGroupInput("includesamples", "", inline = T, choices = list()),
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
                               #rHandsontableOutput("hot"))),
                               #DT::dataTableOutput("sampletable"))),
                    column(width = 7,
                           box(width = NULL,
                               infoBoxOutput("datainfoBox"),
                               infoBoxOutput("sampleinfoBox"),
                               infoBoxOutput("conditioninfoBox")
                           ),
                           box(title = "Expression levels",
                               width = NULL,
                               plotOutput("violinplot"),
                               radioButtons("violintype", "", choices = list("By condition" = 1, "By sample" = 2), inline = T)
                           ))
                    )),



        # Data summary

        # tabItem(tabName = "data-summary",
        #         box(title = "NA frequencies", status = "warning", solidHeader = F,
        #             plotOutput("nafreq"))
        # ),


        
        # Filters: NA cutoff
        
        tabItem(tabName = "filters",
                fluidRow(
                    column(width = 4,
                           box(
                               width = NULL,
                               title = "Missing values cutoff", status = "primary", solidHeader = F,
                               helpText("Set maximum number allowed missing values for each condition."),
                               br(),
                               helpText("To proceed without removing any proteins with missing values, set the cutoff to a value higher than the number of samples / conditions"),
                               br(),
                               numericInput("missingvalues", label = "",
                                            min = 0, max = 9999, value = 1), # max = number of samples / conditions
                               actionButton("loadfilterplot", "Show plot"),
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


        # Data type: distributions
        
        tabItem(tabName = "validateIDs",
                fluidRow(
                    column(width = 6,
                        box(width = NULL,
                            title = "Validate IDs", status = "primary", solidHeader = F,
                            p(helpText("NB: Depending on the number of IDs, this process may take a long time to run.")),
                            actionButton("listCandidates", label = "List outdated IDs"),
                            p(),
                            hr(),
                            DT::DTOutput("obsolete")
                        )
                    ),
                    column(width = 6,
                           box(width = NULL,
                               title = "Convert IDs",
                               textAreaInput("idstoconvert", "Identifiers"),
                               actionButton("convertids", label = "Convert IDs")
                               )
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

        # Data type: distributions
        
        tabItem(tabName = "file-export",
                box(
                    title = "Generate report", status = "primary", solidHeader = F,
                    helpText("Exports generated content to an Rmarkdown document"),
                    h6(actionLink("selectall","Select all"), strong(" | "), actionLink("deselectall","Deselect all")),
                    checkboxGroupInput("export_filtering", "Filtering", inline = T, choices = c("Missing values")),
                    checkboxGroupInput("export_data_inspection", "Data inspection", inline = T, choices = c("PCA 2D", "PCA 3D", "UMAP", "Heatmap")),
                    checkboxGroupInput("export_de", "Differential analysis", inline = T, choices = c("Fold-change", "P-values")),
                    
                    textAreaInput("additionalNotes", "Additional notes"),
                    
                    downloadButton("downloadReport", "Generate report")
                    
                )
        ),
        
        # Quality control: PCA: 2D, 3D
        
        tabItem(tabName = "PCA",
                fluidRow(
                    column(width = 3,
                           box(title = "PCA settings",
                               status = "warning",
                               width = NULL,
                               sliderInput("contribs",
                                           "Number of contributors:",
                                           min = 1,  max = 50, value = 5),
                               sliderInput("ellipse",
                                           "Ellipse level:",
                                           min = 0,  max = 1, value = .75),
                               actionButton("loadPCAplots", "Load plots"),
                               p(),
                               hr(),
                               helpText("Please review pairing mode under Settings")
                           )
                          ),
                    column(width = 9,
                           tabBox(width = NULL,
                               tabPanel("PCA 2D", plotOutput("pca2dplot", width = "600px", height = "550px")),
                               tabPanel("PCA 3D", plotly::plotlyOutput("pca3dplot", width = "600px", height = "600px"))))
                           )
                    
                ),
                
                
        tabItem(tabName = "UMAP",
                fluidRow(
                    column(width = 3,
                        box(width = NULL,
                            status = "warning",
                            actionButton("loadUMAP","Load plot")
                        )),
                    column(width = 9,
                        box(title = "UMAP",
                            width = NULL,
                            plotOutput("UMAPplot"))
                    )
                )),
        tabItem(tabName = "samplecorr",
                fluidRow(
                    column(width = 3,
                        box(title = "Heatmap settings",
                            width = NULL,
                            status = "warning",
                            actionButton("generateheatmap","Load plot")
                        )),
                        
                    column(width = 9,
                        box(title = "Sample correlation", width = NULL,
                            plotOutput("samplecorrheatmap", width = "700px", height = "600px"))
                    ))

        ),
        
        # Differential expression analysis : ANOVA, Contrasts
                #DT::dataTableOutput("contrasttable") ),
        
        tabItem(tabName = "contrasts",
                box(
                    title = "Set contrasts", status = "primary", solidHeader = F,
                    helpText("Which groups should be differentially expressed?"),
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
                    # radioButtons("pairing", label = "Design",
                    #              choices = list("Paired" = 1,
                    #                             "Unpaired" = 2),
                    #              selected = 2,
                    #              inline = T),
                    actionButton("setContrast", "Select"),
                    p(),
                    hr(),
                    helpText("Please review pairing mode under Settings")
                )
        ),
        tabItem(tabName = "diffexpoutput",
                plotOutput("contrasttable", width = "800px", height = "1600px")),
        tabItem(tabName = "differentialexpression",
                fluidRow(
                    column(width = 12,
                    # box(width = NULL,
                    #     downloadButton('download',"Download"),
                    #     tags$p(),
                    #     tabPanel("Table", DT::dataTableOutput("diffexptable")))
                    
                    box(width = NULL,
                        downloadButton('download',"Download")),
                    
                    box(width = NULL,
                        helpText("Display only proteins that have:"),
                        br(),
                        numericInput("diffexp_limit_fc", "Abs. log2 fold-change greater than or equal to", min = 0, max = 50, value = 0, step = .5),
                        numericInput("diffexp_limit_pval", "Adj. P-value less than", min = 0, max = 1, value = 1, step = .1),
                        ),
                    
                    tabBox(width = NULL,
                           tabPanel("Summary statistics", tableOutput("diffexptable_summary")),
                           tabPanel("Up-regulated genes", DT::dataTableOutput("diffexptable_up")),
                           tabPanel("Down-regulated genes", DT::dataTableOutput("diffexptable_down")))
                    ))
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
                               htmlOutput("number_of_genes"),
                               actionButton(inputId = "generate_pathways", label = "Generate pathway data")
                               ),
                           
                           box(title = "Selected pathway info", width = NULL,
                               htmlOutput("selected_pathway")
                               )
                           ),
                    column(width = 9,
                           tabBox(width = NULL,
                                  tabPanel("Enrichment of up-regulated genes", DT::dataTableOutput("upregulated_pathways_table", width = 900)),
                                  tabPanel("Enrichment of down-regulated genes", DT::dataTableOutput("downregulated_pathways_table", width = 900)))
                           )
                )
                ,
                
        ),
        tabItem(tabName = "pathwayvisualization",
                fluidRow(
                    box(width = 3,
                        actionButton("loadPathwayPlots", "Load plots")
                    ),
                    tabBox(width = 9,
                    #tabPanel("Similarity matrix", plotOutput("similarity_plot", height = 850)),
                    tabPanel("Volcano plot", plotly::plotlyOutput("volcano_plot", height = 750)),
                    tabPanel("Sankey diagram", networkD3::sankeyNetworkOutput("sankey", height = 750)))
                )),
        tabItem(tabName = "interactions",
                fluidRow(
                    box(
                        width = 6,
                        height = 520,#464,
                        visNetworkOutput("xxxx")
                    ),
                    box(
                        width = 6,
                        height = 520,
                        plotly::plotlyOutput("volcano_plot2")
                    )
                ),
                
                fluidRow(
                    column(width = 3,
                           box(
                               title = "Network settings", width = NULL,
                               # radioButtons("network_layout_options", label = "Layout options",
                               #              choices = list(
                               #                  "Nicely" = 1,
                               #                  "Circle" = 2,
                               #                  "Grid" = 3,
                               #                  "Sphere" = 4,
                               #                  "Randomly" = 5
                               #              ), inline = T),
                               radioButtons(
                                   "network_regulation",
                                   label = "Subset by up- or down-regulation",
                                   choices = list("Up-regulated" = 1, "Down-regulated" = 2, "Both" = 3),
                                   selected = 3,
                                   inline = F),
                               # numericInput(
                               #     "pvaluecutoff",
                               #     label = "Maximum adj. Pvalue",
                               #     min = 0,
                               #     max = 1,
                               #     value = 0.05,
                               #     step = 0.0001
                               # ),
                               numericInput(
                                   "fccutoff",
                                   label = "Minimum abs. log2FC",
                                   min = 0,
                                   max = 100,
                                   value = 0,
                                   step = 0.1
                               ),
                               numericInput(
                                   "interactioncutoff",
                                   label = "Interaction confidence",
                                   min = 0,
                                   max = 10,
                                   value = 7,
                                   step = .1
                               )
                           )
                    ),
                    
                    column(width = 3,
                           box(
                               title = "Selection criteria", width = NULL,
                               
                               radioButtons("interaction_behaviour", label = "Selection subset",
                                            choices = list("Strict" = 1,
                                                           "Extended" = 2),
                                            inline = T,
                                            selected = 1),
                               textAreaInput("network_proteins", label = "Select proteins"),
                               actionButton("highlight_nodes", label = "Highlight")
                           )
                    ),
                    
                    column(width = 3,
                           box(
                               title = "Protein info", width = NULL,
                               uiOutput("clicked_node"),
                               uiOutput("hovered_node")
                           )
                           )
                )
        ),
        
        tabItem(tabName = "contact",
                fluidRow(
                    column(width = 5,
                           box(title = "Contact", width = NULL,
                               helpText("Found a bug? Use this form to contact me about anything: errors, feature requests, ideas or other comments."),
                               textInput("cname", "Name", placeholder = "Anonymous"),
                               textInput("email", "Email", placeholder = "Anonymous"),
                               textAreaInput("suggestion", "Comment"),
                               actionButton("sendcomment", "Send", icon = icon("paper-plane"))
                           )
                    )
                )
        ),
        
        
        tabItem(tabName = "settings",
                fluidRow(
                    column(width = 5,
                        box(title = "Settings", width = NULL,
                            h4("Display options"),
                            hr(),
                            radioButtons(
                                inputId = 'colorScheme',
                                label = 'Preferred color scheme',
                                choices = list(
                                    "Normal" = 1,
                                    "Colorblind friendly" = 2
                                ),
                                inline = T
                            ), 
                            radioButtons(inputId = 'textSize', label = 'Text size',
                                         choices = list("Small", "Medium", "Large"), selected = "Medium", inline = T),
                            h4("Target organism"),
                            hr(),
                            selectInput("species", "Select species",
                                        list("Species" = list("Human (Homo sapiens)"))),
                            h4("Identifier type"),
                            hr(),
                            selectInput("displayIdentifier",
                                        label = "Display ID labels as",
                                        choices = list(
                                            "UniProtKB" = 2,
                                            "Entrez" = 3,
                                            "Gene Symbol" = 4,
                                            "Ens. Gene ID" = 5,
                                            "Ens. Protein ID" = 6
                                        ),
                                        selected = 1),
                            h4("Differential expression"),
                            hr(),
                            selectInput("setDEengine",
                                        label = "Set engine",
                                        choices = list("Limma version 3.39.1" = 1,
                                                       "DESeq2 version 3.10" = 2),
                                        selected = 1),
                            radioButtons("diffexppairing", "Pairing", choices = list("Paired" = 1, "Unpaired" = 2), inline = T, selected = 2)
                            )
                    )
                    )
                )
    ),
    div(class = "footer_wrapper",
        div(class = "sticky_footer", span("© 2020 · ProteoMill | Martin Rydén", style = "line-height: 35px; float: right; margin-right: 20px;")))
)

# Load dashboard page ----

dashboardPage(title = "ProteoMill | Differential expression pathway and network analysis tool", skin = 'black', header, sidebar, body)

