
library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
    h3("Design a dataset"),
    helpText("To learn how the input data file should be formatted, play around with the settings below and generate simulated datasets. When you have set the conditions and replicates according to your experiment, you can download the file and use the header in the original data file."),
    numericInput("conditions", "Number of conditions", min = 2, max = 6, value = 2),
    numericInput("replicates", "Number of replicates", min = 1, max = 20, value = 2),
    selectInput("identifiers", "Identifier type", choices = list(
        "UniProtKB" = 2,
        "Entrez" = 3,
        "Gene Symbol" = 4,
        "Ens. Gene ID" = 5,
        "Ens. Protein ID" = 6
    )),
    
    actionButton("generate", "Generate dataset"),
    p(),
    
    tableOutput("generatedDataset"),
    
    downloadLink('downloadData', 'Save to file'),
    p()
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$downloadData <- downloadHandler(
        filename = function() {
            paste("data-", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
            write.csv(n, file)
        }
    )
    
    UNIPROTID <- c("A1L4H1", "A5A3E0", "A6NMZ7", "A8TX70", "O00142", "O00299", "O00300", "O00339", "O00391", "O00622", "O14791", "O14818", "O14960", "O14974", "O15144", "O15230", "O15240", "O15335", "O15511", "O43364", "O43488", "O43491", "O43586", "O43707", "O43827", "O43852", "O43854", "O60565", "O60664", "O60687", "O60701", "O75339", "O75365", "O75368", "O75369", "O75462", "O75593", "O75874", "O76076", "O95631", "O95810", "O95843", "O95865", "P00325", "P00338", "P00352", "P00367", "P00387", "P00441", "P00450")
    
    ENTREZID <- c("284297", "728378", "131873", "256076", "7084", "1192", "4982", "4147", "5768", "3491", "8542", "5688","3950", "4659", "10109", "100616340", "7425", "1101", "10092", "3199", "8574", "2037", "9051", "81", "10218", "813", "10085", "26585", "10226", "27286", "7358", "8483", "11156", "6451", "2317", "9244", "8928", "3417", "8839", "9423", "8436", "9626", "23564", "125", "3939", "216", "2746", "1727","6647", "1356")
    
    SYMBOL <- c("SSC5D", "POTEF", "COL6A6", "COL6A5", "TK2", "CLIC1", "TNFRSF11B", "MATN2", "QSOX1", "CYR61", "APOL1", "PSMA7", "LECT2", "PPP1R12A", "ARPC2", "LAMA5", "VGF", "CHAD", "ARPC5", "HOXA2", "AKR7A2", "EPB41L2", "PSTPIP1", "ACTN4", "ANGPTL7", "CALU", "EDIL3", "GREM1", "PLIN3", "SRPX2", "UGDH", "CILP", "PTP4A3", "SH3BGRL", "FLNB", "CRLF1", "FOXH1", "IDH1", "WISP2", "NTN1", "SDPR", "GUCA1C", "DDAH2", "ADH1B", "LDHA", "ALDH1A1", "GLUD1", "CYB5R3", "SOD1", "CP")
    
    GENEID <- c("ENSG00000179954", "ENSG00000196604", "ENSG00000206384", "ENSG00000172752", "ENSG00000166548", "ENSG00000206394", "ENSG00000164761", "ENSG00000132561", "ENSG00000116260", "ENSG00000142871", "ENSG00000100342", "ENSG00000101182", "ENSG00000145826", "ENSG00000058272", "ENSG00000163466", "ENSG00000130702", "ENSG00000128564", "ENSG00000136457", "ENSG00000162704", "ENSG00000105996", "ENSG00000053371", "ENSG00000079819", "ENSG00000140368", "ENSG00000130402", "ENSG00000171819", "ENSG00000128595", "ENSG00000164176", "ENSG00000166923", "ENSG00000105355", "ENSG00000102359", "ENSG00000109814", "ENSG00000138615", "ENSG00000184489", "ENSG00000131171", "ENSG00000136068", "ENSG00000006016", "ENSG00000160973", "ENSG00000138413", "ENSG00000064205", "ENSG00000065320", "ENSG00000168497", "ENSG00000138472", "ENSG00000206395", "ENSG00000196616", "ENSG00000134333", "ENSG00000165092", "ENSG00000148672", "ENSG00000100243", "ENSG00000142168", "ENSG00000047457")
    
    PROTEINID <- c("ENSP00000467252", "ENSP00000386786", "ENSP00000351310", "ENSP00000309762", "ENSP00000440898", "ENSP00000382926", "ENSP00000297350", "ENSP00000429977", "ENSP00000356574", "ENSP00000398736", "ENSP00000380448", "ENSP00000359910", "ENSP00000274507", "ENSP00000389168", "ENSP00000327137", "ENSP00000252999", "ENSP00000249330", "ENSP00000423812", "ENSP00000352918", "ENSP00000222718", "ENSP00000235835", "ENSP00000434308", "ENSP00000452746", "ENSP00000252699", "ENSP00000366015", "ENSP00000442110", "ENSP00000296591", "ENSP00000478319", "ENSP00000221957", "ENSP00000362095", "ENSP00000319501", "ENSP00000261883", "ENSP00000428976", "ENSP00000362308", "ENSP00000420213", "ENSP00000376188", "ENSP00000366534", "ENSP00000260985", "ENSP00000361959", "ENSP00000173229", "ENSP00000305675", "ENSP00000261047", "ENSP00000382935", "ENSP00000306606", "ENSP00000395337", "ENSP00000297785", "ENSP00000277865", "ENSP00000338461", "ENSP00000270142", "ENSP00000264613")
    
    
    observeEvent(input$generate, {
        
        identifier <- function(i) {
            switch(i,
                   "2" = "UNIPROTID",
                   "3" = "ENTREZID",
                   "4" = "SYMBOL",
                   "5" = "GENEID",
                   "6" = "PROTEINID")
        }
        
        conditions <- input$conditions
        replicates <- input$replicates
        identifiers <- identifier(input$identifiers)
        
        n <- matrix(rnorm(50, 200, 50), 50, conditions * replicates)
        rownames(n) <- get(identifiers)
        colnames(n) <- paste0(rep(LETTERS[1:conditions], each=replicates), "_", rep(1:replicates))
        
        assign("n", n, envir = .GlobalEnv)
        
        output$generatedDataset <- renderTable({n}, rownames = T)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

