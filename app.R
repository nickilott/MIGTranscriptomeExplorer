library(shiny)
library(MIGTranscriptomeExplorer)

# Define UI for MIGTranscriptomeDB app
ui <- pageWithSidebar(
  
    # App title
    headerPanel("MIGTranscriptomeDB: Database of transcriptomic data sets"),
  
    # Sidebar panel for inputs

    sidebarPanel(
        actionButton("show.datasets", "show datasets"),
        tableOutput('dataset.table'),
        textInput("gene", label="Gene:", value = ""),
        actionButton("gene.search", "get expression")),

    # Main panel for displaying outputs
    mainPanel()
)

# Define server logic
server <- function(input, output) {

    output$dataset.table <- renderDataTable({
        input$show.datasets
        showDatasets(connect())})
}

shinyApp(ui, server)