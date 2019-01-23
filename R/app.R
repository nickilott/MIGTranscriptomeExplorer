library(shiny)
library(MIGTranscriptomeExplorer)

# Define UI for MIGTranscriptomeDB app
ui <- pageWithSidebar(
  
  # App title
  headerPanel("MIGTranscriptomeDB: Database of transcriptomic data sets"),
  
  # Sidebar panel for inputs

  sidebarPanel(
       actionButton("show_datasets", "show datasets")
       textInput("word", label="Gene:", value = ""),
       actionButton("gene_search", "get expression")),
  
  # Main panel for displaying outputs
  mainPanel()
)

# Define server logic
server <- function(input, output) {
  
  
  

}

shinyApp(ui, server)