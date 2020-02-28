library(shiny)
library(shinydashboard)

#setwd("~/qodb-shiny/")

#tissues <- list()

# Notification menus ----

help <- dropdownMenuOutput("helpMenu")
notifications <- dropdownMenuOutput("notifMenu")

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
                 menuSubItem("File input", tabName = "file-input", href = NULL, newtab = TRUE,
                             icon = shiny::icon("angle-double-right"), selected = T),
                 menuSubItem("Filters", tabName = "filters", href = NULL, newtab = TRUE,
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
        menuItemOutput('interactions'),
        menuItemOutput('predictiverm'),
        menuItemOutput('predictive')
        
        
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
        menuItem("About", tabName = "about", icon = icon("book")),
        menuItem("Settings", tabName = "settings", icon = icon("sliders-h"))
    )
)


# Main content ----

body <- dashboardBody(
    
    # Import css and js
    tags$head(
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
    
    #tags$div(class = "overlay",
    #         tags$h1(class = "ml9",
    #                 tags$span(class = "text-wrapper",
    #                           tags$span(class = "letters", "Quantitative Omics Discovery Base"))),
    #         tags$div(class = 'begindiv',
    #                  actionLink('removeBanner', label = 'LAUNCH')
    #         )
    #),
    tags$script(src = "custom.js"),
    tags$script(src = "animate-notifs.js"),

    tabItems(
        # File input
        tabItem(tabName = "file-input",
                fluidRow(
                    tabBox(width = 4,
                        tabPanel(
                            title = "Dataset",
                            footer = helpText('Accepted filetypes are csv, tsv and txt.'),
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
                            actionButton("useDemoData", label = "Use demo data")
                            #DT::dataTableOutput("contrasttable")
                        )
                        
                    )
                )
),
        
        # Filters: NA cutoff
        
        tabItem(tabName = "filters",
                box(
                    title = "Missing values cutoff", status = "primary", solidHeader = F, width = 3,
                    helpText("Set maximum number allowed missing values per condition."),
                    numericInput("missingvalues", label = "NA cutoff",
                                 min = 0, max = 50, value = 1), # max = number of samples / conditions
                    actionButton("loadfilterplot", "Load plot"),
                    actionButton("setcutoff", "Set cutoff")
                ),
                box(title = "NA frequencies", status = "warning", solidHeader = F,
                    plotOutput("nafreq"))
        ),


        # Data type: distributions
        
        tabItem(tabName = "validateIDs",
                box(
                    title = "Validate IDs", status = "primary", solidHeader = F,
                    helpText("NB: Depending on the number of IDs, this process may take a long time to run."),
                    actionButton("listCandidates", label = "List outdated IDs"),
                    br(),
                    hr(),
                    DT::DTOutput("obsolete")
                )
        ),
                
        
        # Data type: distributions
        
        tabItem(tabName = "goodnessOfFit",
                box(
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
                    actionButton("compile_data", "Generate")
                )
        ),
        
        # Quality control: PCA: 2D, 3D
        
        tabItem(tabName = "PCA",
                box(title = "PCA settings",
                    status = "warning",
                    width = 3,
                    sliderInput("contribs",
                                "Number of contributors:",
                                min = 1,  max = 300, value = 10),
                    sliderInput("ellipse",
                                "Ellipse level:",
                                min = 0,  max = 1, value = .75),
                    actionButton("loadPCAplots", "Load plots")
                ),
                tabBox(
                    tabPanel("PCA 2D", plotOutput("pca2dplot", width = "600px", height = "520px")),
                    tabPanel("PCA 3D", plotly::plotlyOutput("pca3dplot")))),
        tabItem(tabName = "UMAP",
                fluidRow(
                    box(width = 3,
                        status = "warning",
                        actionButton("loadUMAP","Load plot")
                    ),
                    box(title = "UMAP",
                        plotOutput("UMAPplot"))
                )),
        tabItem(tabName = "samplecorr",
                fluidRow(
                        box(title = "Heatmap settings",
                            width = 3,
                            status = "warning",
                            actionButton("generateheatmap","Generate Heatmap")
                        ),
                        
                        box(title = "Sample correlation", width = 6,
                            plotOutput("samplecorrheatmap", width = "700px", height = "600px"))
                    )

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
                    radioButtons("pairing", label = "Design",
                                 choices = list("Paired" = 1,
                                                "Unpaired" = 2),
                                 selected = 2,
                                 inline = T),
                    actionButton("setContrast", "Select")
                )
        ),
        tabItem(tabName = "diffexpoutput",
                plotOutput("contrasttable", width = "800px", height = "1600px")),
        tabItem(tabName = "differentialexpression",
                fluidRow(
                    box(actionButton("loadDiffExpTable","Load table"), width = 3),
                    box(tabPanel("Table", DT::dataTableOutput("diffexptable")))
                    )
                ),

        # Pathway enrichment: Table, Similarity matrix, Volcano plot
        
        tabItem(tabName = "pathwayenrichment",
                box(title = "Pathway settings", status = "warning", width = 3,
                    selectInput(inputId = "pathdb", label = "Pathway database",
                                choices = list("REACTOME" = 'REACTOME'),
                                selected = 'REACTOME'),
                    radioButtons(inputId = "abstractionlevel", label = "Level of abstraction",
                                 choices = list("Global", "Lowest"),
                                 selected = "Lowest", inline = T),
                    numericInput("min_fc", "Min. log2 fold change", value = 0, min = 0, max = 50, step = .5),
                    numericInput("min_pval", "Min. adj. P-value", value = 0.05, min = 0, max = 1, step = .01),
                    htmlOutput("number_of_genes"),
                    radioButtons(inputId = "usebackground", label = "Background genes",
                                 choices = list("My dataset" = 1, "Extended background" = 2, "No background (entire genome)" = 3),
                                 selected = 2),
                    actionButton(inputId = "generate_pathways", label = "Generate pathway data")
                    ),
                tabBox(width = 9,
                       tabPanel("Enrichment of up-regulated genes", DT::dataTableOutput("upregulated_pathways_table", width = 900)),
                       tabPanel("Enrichment of down-regulated genes", DT::dataTableOutput("downregulated_pathways_table", width = 900)))
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
                        height = 550,#464,
                        visNetworkOutput("xxxx")
                    ),
                    tabBox(width = 6,
                           height = 550,
                           tabPanel("volcano_network_tab", plotly::plotlyOutput("volcano_plot2")),
                           tabPanel("pca_network_tab", "PCA goes here")
                    )
                ),
                
                fluidRow(
                    column(width = 3,
                           box(
                               title = "Network settings", width = NULL,
                               radioButtons("network_layout_options", label = "Layout options",
                                            choices = list(
                                                "Nicely" = 1,
                                                "Circle" = 2,
                                                "Grid" = 3,
                                                "Sphere" = 4,
                                                "Randomly" = 5
                                            ), inline = T),
                               radioButtons(
                                   "network_regulation",
                                   label = "Subset by up- or down-regulation",
                                   choices = list("Up-regulated" = 1, "Down-regulated" = 2, "Both" = 3),
                                   selected = 3,
                                   inline = T),
                               numericInput(
                                   "pvaluecutoff",
                                   label = "Maximum adj. Pvalue",
                                   min = 0,
                                   max = 1,
                                   value = 0.05,
                                   step = 0.0001
                               ),
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
                                                           "Extended" = 2,
                                                           "Full range" = 3),
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
        tabItem(tabName = "predictive",
                fluidRow(
                    box(title = "Predictive network analysis",
                        width = 6
                    )
                )
            ),
        tabItem(tabName = "settings",
                fluidRow(
                    box(title = "Settings", width = 3,
                        # h5("Display options"),
                        # hr(),
                        # radioButtons(inputId = 'colorScheme', label = 'Color palette',
                        #              choices = list("Normal", "Colorblind friendly"), inline = T),
                        # radioButtons(inputId = 'textSize', label = 'Text size',
                        #              choices = list("Small", "Medium", "Large"), selected = "Medium", inline = T),
                        h5("Target organism"),
                        hr(),
                        selectInput("species", "Select species",
                                    list("Species" = list("Human (Homo sapiens)"))),
                        h5("Identifier type"),
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
                        h5("Differential expression"),
                        hr(),
                        selectInput("setDEengine",
                                    label = "Set engine",
                                    choices = list("Limma version 3.39.1" = 1,
                                                   "DESeq2 version 3.10" = 2),
                                    selected = 1)
                        )
                    )
                )
    )
    )

# Load dashboard page ----

dashboardPage(title = "ProteoMill | Differential expression pathway and network analysis tool", skin = 'black', header, sidebar, body)