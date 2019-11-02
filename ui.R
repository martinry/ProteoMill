library(shiny)
library(shinydashboard)

setwd("~/qodb-shiny/")

tissues <- list()

# Notification menus ----

help <- dropdownMenuOutput("helpMenu")
notifications <- dropdownMenuOutput("notifMenu")

header <- dashboardHeader(help,
                          notifications,
                          title = "qodb",
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
    sidebarMenu(
        menuItem("Overview", tabName = "overview", icon = icon("align-left"), selected = T),
        menuItem("Dataset options", icon = icon("table"),
                 menuSubItem("File input", tabName = "file-input", href = NULL, newtab = TRUE,
                             icon = shiny::icon("angle-double-right"), selected = F),
                 menuSubItem("ID Mapping", tabName = "id-mapping", href = NULL, newtab = TRUE,
                             icon = shiny::icon("angle-double-right"), selected = F),
                 menuSubItem("Filters", tabName = "filters", href = NULL, newtab = TRUE,
                             icon = shiny::icon("angle-double-right"), selected = F),
                 menuSubItem("Data type", tabName = "datatype", href = NULL, newtab = TRUE,
                             icon = shiny::icon("angle-double-right"), selected = F)
        ),
        menuItemOutput("qualityrm"),
        menuItemOutput("quality"),
        menuItemOutput('diffrm'),
        menuItemOutput("differential"),
        menuItemOutput('enrichrm'),
        menuItemOutput('enrichment'),
        menuItemOutput('networkrm'),
        menuItemOutput('network')
        
        
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
    sidebarMenu(
        menuItem("Structures", tabName = "structures", icon = icon("fingerprint"))
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
    sidebarMenu(
        menuItem("News", tabName = "", icon = icon("book")),
        menuItem("Functions", tabName = "", icon = icon("book")),
        menuItem("About", tabName = "", icon = icon("book"))
    )
)


# Main content ----

body <- dashboardBody(
    
    # Import custom stylesheet and javascript
    
    #<link rel="stylesheet" href="~/martiniry/LiteMol/css/LiteMol-plugin.css" type="text/css" />
    #<script src="~/martiniry/LiteMol/js/LiteMol-plugin.js"></script>
    
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
        tags$link(rel="stylesheet", href="https://fonts.googleapis.com/css?family=Quicksand"),
        tags$link(rel="stylesheet", href="https://fonts.googleapis.com/css?family=Poiret+One"),
        tags$link(rel="stylesheet", href="https://fonts.googleapis.com/css?family=Open+Sans"),
        tags$link(rel="stylesheet", href="css/LiteMol-plugin-light.css"),
        tags$script(src = "js/LiteMol-plugin.js?lmversion=10"),
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

    # Overview: settings: colors, font
    
    tabItems(
        tabItem(tabName = "overview",
                box(title = 'General settings', status = 'primary', solidHeader = F,
                    hr(),
                    radioButtons(inputId = 'colorScheme', label = 'Color palette',
                                 choices = list("Normal", "Colorblind friendly"), inline = T),
                    radioButtons(inputId = 'textSize', label = 'Text size',
                                 choices = list("Small", "Medium", "Large"), selected = "Medium", inline = T),
                    hr(),
                    selectInput("species", "Select species",
                                list("Species" = list("Human (Homo sapiens)")))
                    
                )
        ),
        
        # File input
        tabItem(tabName = "file-input",
                tabBox(
                    tabPanel(
                        title = "Dataset",
                        footer = helpText('Accepted filetypes are csv, tsv and txt.'),
                        radioButtons("sep", label = 'Separator',
                                     choices = list("Comma" = 1, "Semicolon" = 2, "Tab" =3),
                                     selected = 1,
                                     inline = T),
                        fileInput("infile", "Select a file",
                                  accept = c(
                                      "text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                        helpText('Accepted filetypes are csv, tsv and txt.')
                        ),
                    tabPanel(
                        title = 'Annotation data',
                        radioButtons("anno_sep", label = 'Separator',
                                     choices = list("Comma" = 1, "Semicolon" = 2, "Tab" =3),
                                     selected = 1,
                                     inline = T),
                        fileInput("anno_infile", "Select a file",
                                  accept = c(
                                      "text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                        helpText('Accepted filetypes are csv, tsv and txt.')
                        )
                    )),
        
        # Filters: NA cutoff
        
        tabItem(tabName = "filters",
                box(
                    title = "Missing values cutoff", status = "primary", solidHeader = F, width = 3,
                    helpText("Set maximum number allowed NA per condition."),
                    numericInput("missingvalues", label = "NA cutoff",
                                 min = 0, max = 50, value = 1), # max = number of samples / conditions
                    actionButton("setcutoff", "Set cutoff")
                ),
                box(title = "NA frequencies", status = "warning", solidHeader = F,
                    plotOutput("nafreq"))
        ),
        
        
        # Convert IDs
        
        tabItem(tabName = "id-mapping",
                box(title = "Identifiers", width = 3,
                selectInput("sourceIDtype",
                            label = "Source ID",
                            choices = list("UniProtKB" = 1,
                                           "Entrez" = 2,
                                           "Ensembl gene ID" = 3,
                                           "Gene symbol" = 4),
                            selected = 1),
                actionButton("verifyIDs", label = "Check IDs"), textOutput("runningprocesstext")
                ),
                box(
                    title = "Convert to", status = "primary", solidHeader = F, width = 3,
                    helpText("The gene identifier to be used."),
                    selectInput("identifiers", label = "Identifier type",
                                choices = list("Entrez" = 1,
                                               "UniProtKB" = 2,
                                               "Ensembl gene ID" = 3,
                                               "Gene symbol" = 4), 
                                selected = 1),
                    actionButton("convertIDs", "Convert")
                )
                
        ),
                
        
        # Data type: distributions
        
        tabItem(tabName = "datatype",
                box(
                    title = "Distributions", status = "primary", solidHeader = F,
                    selectInput(inputId = "fit", label = "Fit to distribution",
                                choices = list("Normal" = 1,
                                               "Poisson" = 2,
                                               "Negative binomial" = 3),
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
                    
                ),
                box(
                    title = "Data type and distribution assumptions", status = "primary", solidHeader = F,
                    helpText("What type of experiment best describes your data?"),
                    selectInput("groups", label = "Conditions",
                                choices = list("RNA-seq" = 1,
                                               "Single-cell RNA" = 2,
                                               "Microarray" = 3,
                                               "Mass spectrometry" = 4),
                                selected = 1),
                    actionButton("useIDs", "Select")
                )
        ),
        
        # Quality control: PCA: 2D, 3D
        
        ## To do: Make tab box. Fix contribs.
        
        tabItem(tabName = "PCA",
                box(title = "PCA settings",
                    status = "warning",
                    width = 3,
                    sliderInput("contribs",
                                "Number of contributors:",
                                min = 1,  max = 300, value = 10),
                    sliderInput("ellipse",
                                "Ellipse level:",
                                min = 0,  max = 1, value = .75)
                ),
                tabBox(
                    tabPanel("PCA 2D", plotOutput("pca2dplot", width = "600px", height = "520px")),
                    tabPanel("PCA 3D", plotly::plotlyOutput("pca3dplot")))),
        tabItem(tabName = "samplecorr",
                box(title = "Sample correlation",
                    plotOutput("samplecorrheatmap", width = "700px", height = "600px"))
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
                tabPanel("Table", DT::dataTableOutput("diffexptable", width = 1000))),
        
        # Select background data for enrichment

        tabItem(tabName = "bgdata",
                box(title = "Add custom background data",
                    selectInput("Tissue", "Select a tissue type:",
                                list(`Dataset: Santos, A et al. (2015)` = tissue_names)),
                    actionButton("addbg", "Add to background"),
                    actionButton("resetbg", "Reset background"))
        ),
        
        # Pathway enrichment: Table, Similarity matrix, Volcano plot
        
        tabItem(tabName = "pathwayenrichment",
                box(title = "Pathway settings", status = "warning", width = 3,
                    selectInput(inputId = "pathdb", label = "Pathway database",
                                choices = list("REACTOME" = 'REACTOME', "KEGG" = 'KEGG', "BIOCARTA" = 'BIOCARTA'),
                                selected = 'REACTOME'),
                    checkboxInput(inputId = "trimprefix", label = "Hide prefixes", value = T),
                    radioButtons(inputId = "abstractionlevel", label = "Level of abstraction",
                                 choices = list("Global" = 1, "Lowest" = 2),
                                 selected = 1, inline = T),
                    radioButtons(inputId = "usebackground", label = "Background genes",
                                 choices = list("My dataset" = 1, "Extended background" = 2, "No background (entire genome)" = 3),
                                 selected = 2),
                    actionButton(inputId = "generatepathways", label = "Generate pathway data")
                    ),
                tabBox(width = 9,
                       tabPanel("Table", DT::dataTableOutput("pathtable", width = 900)),
                       tabPanel("Similarity matrix", plotOutput("similarity_plot", height = 850)),
                       tabPanel("Volcano plot", plotOutput("volcano_plot", height = 750)),
                       tabPanel("Sankey diagram", networkD3::sankeyNetworkOutput("sankey", height = 750)))
        ),
        
        tabItem(tabName = "network",
                box(title = "Network settings", width = 3,
                    radioButtons("pathwaylevel", label = "Pathway annotation level",
                                 choices = list("Highest" = 1, "Lowest" = 2), inline = T),
                    helpText("Differential expression"),
                    hr(),
                    radioButtons("direction", label = "Regulation", choices = list("Up-regulation" = 1, "Down-regulation" = 2), inline = T),
                    numericInput(
                        "pvaluecutoff",
                        label = "Maximum adj. Pvalue",
                        min = 0,
                        max = 1,
                        value = 0.0001,
                        step = 0.0001
                    ), 
                    numericInput(
                        "fccutoff",
                        label = "Minimum log2FC",
                        min = 0,
                        max = 100,
                        value = 1,
                        step = 0.1
                    ), 
                    helpText("Interactions"),
                    hr(),
                    sliderInput("interactioncutoff", label = "Minimum interaction score", min = 0, max = 9.9, value = 3, step = .1),
                    actionButton("generatenetwork", label = "Generate network")),
                # box(title = "Structure",
                #     tags$div(id = "litemol", style = 'width: 640px; height: 480px; margin-top: 200px; position: relative')),
                box(title = "Interactions", width = 9,
                    networkD3::forceNetworkOutput("net", width = "100%", height = "750px"))),
        tabItem(tabName = "structures",
                # box(title = "Search protein", width = 2,
                #     textInput("proteinsearch", label = "Search"),
                #     actionButton("searchclick", "Find structure")
                #     ),
                box(title = "Structure", width = 10,
                    textInput("seachinput", "Search"),
                    actionButton("searchclick", "Find structure"),
                    tags$div(id = "litemol", style = 'width: 640px; margin-top: 10px; height: 480px; position: relative'),
                    tags$script(HTML("var plugin = LiteMol.Plugin.create({ target: '#litemol', layoutState: { hideControls: true } });")))
                )
    )
    )

# Load dashboard page ----

dashboardPage(skin = 'black', header, sidebar, body)