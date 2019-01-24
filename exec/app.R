library(shiny)
library(MIGTranscriptomeExplorer)
library(gridExtra)
library(RSQLite)

# Define UI for MIGTranscriptomeDB app
ui <- pageWithSidebar(
  
    # App title
    headerPanel("MIGTranscriptomeDB: Database of transcriptomic data sets"),
  
    # Sidebar panel for inputs

    sidebarPanel(
        actionButton("show.datasets", "show datasets"),
        textInput("gene", label="Gene:", value = ""),
        actionButton("gene.search", "get expression")),

    # Main panel for displaying outputs
    mainPanel(dataTableOutput('dataset.table'),
              plotOutput("gene.expression"))
)


db <- system.file("data/csvdb", package="MIGTranscriptomeExplorer")
conn <- connect(db=db)

# Define server logic
server <- function(input, output) {

    # displaying datasets
    data <- eventReactive(input$show.datasets,{

    df <- showDatasets(conn)
    df
    })

    output$dataset.table <- renderDataTable({
    data()
    })

    # gene expression
    expression <- eventReactive(input$gene.search, {

    statement <- 'SELECT dataset FROM reference'
    datasets <- dbGetQuery(conn, statement)$dataset

    grobs.list <- list()
    for (i in 1:length(datasets)){
        dataset <- datasets[i]
        expression <- getExpression(conn, dataset, input$gene) 
        metadata <- getMetadata(conn, dataset)
        p <- plotGeneOfInterest(dataset, expression, metadata, variable="treatment")
        grobs.list[[i]] <- p
    }
    grid.arrange(grobs=grobs.list, ncol=length(datasets))
    })

    output$gene.expression <- renderPlot({
    expression()
    })

}

shinyApp(ui, server)