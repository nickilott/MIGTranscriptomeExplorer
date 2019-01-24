library(shiny)
library(MIGTranscriptomeExplorer)
library(gridExtra)
library(RSQLite)
library(dplyr)

# Define UI for MIGTranscriptomeDB app

db <- system.file("data/csvdb", package="MIGTranscriptomeExplorer")
conn <- connect(db=db)

ui <- pageWithSidebar(
  
    # App title
    headerPanel("MIGTranscriptomeExplorer: Exploring database of transcriptomic data sets"),
  
    # Sidebar panel for inputs

    sidebarPanel(
	h4("Available datasets"),
        actionButton("show.datasets", "show datasets"),

	h4("Search for gene in database"),
	textInput("gene", label="Gene:", value = ""),
        actionButton("gene.search", "get expression"),

	h4("Significant contrasts"),
	h5("Specify thresholds"),
	numericInput("lfc", label="lfc:", value = 0),
	numericInput("padj", label="padj", value = 0.05),
	actionButton("significant", "get significant"),

	h4("Choose dataset"),
	selectInput("choose.dataset", "Dataset:", choices=list(showDatasets(conn)$dataset)),

	h4("Principle Components Analysis"),
	actionButton("PCA", "PCA"),

	h4("Colour by"),
	uiOutput("variable"),

	h5("Plot expression vs. fold change"),
	actionButton("MA", "MA plot"),

	uiOutput("ma.contrasts"),
	h5("Thresholds"),
	numericInput("ma.lfc", label="lfc", value = 1)
	),

    # Main panel for displaying outputs
    mainPanel(dataTableOutput("dataset.table"),
              plotOutput("gene.expression"),
	      dataTableOutput("significant.results"),
	      plotOutput("PCA"),
	      plotOutput("MA"))
)

db <- system.file("data/csvdb", package="MIGTranscriptomeExplorer")
conn <- connect(db=db)

# Define server logic
server <- function(input, output) {


    ######################
    # displaying datasets
    ######################
    data <- eventReactive(input$show.datasets,{

    df <- showDatasets(conn)
    df
    })

    output$dataset.table <- renderDataTable({
    data()
    })

    #####################
    # gene expression
    #####################
    expression <- eventReactive(input$gene.search, {

    statement <- 'SELECT dataset FROM reference'
    datasets <- dbGetQuery(conn, statement)$dataset

    grobs.list <- list()
    for (i in 1:length(datasets)){
        dataset <- datasets[i]
        expression <- getExpression(conn, dataset, input$gene) 
	if (nrow(expression) == 0){next}
	metadata <- getMetadata(conn, dataset)
        p <- plotGeneOfInterest(dataset, expression, metadata, variable="treatment")
        grobs.list[[i]] <- p
    }
    grid.arrange(grobs=grobs.list, ncol=length(datasets))
    })

    output$gene.expression <- renderPlot({
    expression()
    })

    #####################
    # significant sets
    #####################
    significant <- eventReactive(input$significant, {

    statement <- 'SELECT dataset FROM reference'
    datasets <- dbGetQuery(conn, statement)$dataset

    dfs.list <- list()
    idx <- 1
    for (i in 1:length(datasets)){
        dataset <- datasets[i]

	# get the contrasts for the dataset
	contrasts <- getContrasts(conn, dataset)

	# iterate over contrasts
	for (j in 1:length(contrasts)){
            contrast <- contrasts[j]
            significant <- getSignificant(conn, dataset, contrast, input$lfc, input$padj, input$gene)
	    if (!(is.na(significant)) > 0){
	        significant$contrast <- contrast
	        dfs.list[[idx]] <- significant
		idx <- idx + 1
		}
            }
        }
    df.out <- bind_rows(dfs.list)
    })
    output$significant.results <- renderDataTable({
    significant()
    })

    #####################
    # plots
    #####################

    output$variable <- renderUI({
    variables <- append("none", getMetadataList(conn, input$choose.dataset))
    names(variables) <- variables
    selectInput("variable", "variable", variables)
    })
    

    PCA <- eventReactive(input$PCA, {
    mat <- getMatrix(conn, input$choose.dataset)

    metadata <- getMetadata(conn, input$choose.dataset)
    metadata <- sortMetadata(mat, metadata)
    metadata$none <- "none"

    pc <- runPCA(mat)
    plotPrincipleComponents(pc, metadata, colourby=input$variable)
    })
    output$PCA <- renderPlot({
    PCA()
    })

    output$ma.contrasts <- renderUI({
    contrasts <- getContrasts(conn, input$choose.dataset)
    names(contrasts) <- contrasts
    selectInput("ma.contrasts", "contrast", contrasts)
    })

    MA <- eventReactive(input$MA, {
    result <- getResultSet(conn, input$choose.dataset, input$ma.contrast)
    plotMA(result, lfc=input$ma.lfc, title=paste(input$choose.dataset, input$ma.contrast, sep=": "))
    })
    output$MA <- renderPlot({
    MA()
    })
}

shinyApp(ui, server)