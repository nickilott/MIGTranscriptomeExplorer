# load packages
library(shiny)
library(MIGTranscriptomeExplorer)
library(gridExtra)
library(RSQLite)
library(dplyr)
library(shinythemes)
library(shinyjs)

# Define UI for MIGTranscriptomeExplorer app

db <- system.file("data/csvdb", package="MIGTranscriptomeExplorer")
conn <- connect(db=db)

ui <- fluidPage(theme=shinytheme("flatly"),
  
    # App title
    titlePanel("MIGTranscriptomeExplorer"),
  
    # Sidebar panel for inputs
    sidebarLayout(
        sidebarPanel(
            shinyjs::useShinyjs(),
            id = "side-panel",

	    h2("Explore gene expression across datasets"),

	    h4("Search for gene in database"),
	    textInput("gene", label="Gene:", value = ""),
            actionButton("gene.search", "get expression"),

	    h4("Significant contrasts in database"),
	    h5("Specify thresholds"),
	    numericInput("lfc", label="lfc:", value = 0),
	    numericInput("padj", label="padj", value = 0.05),
	    actionButton("significant", "get significant"),

	    h2("Explore specific dataset"),
	    h4("Choose dataset"),
	    selectInput("choose.dataset", "Dataset:", choices=showDatasets(conn)$dataset),

	    h4("Principle Components Analysis"),

	    selectInput("PCs", "PC:", choices=c("PC1 vs. PC2" = "PC1_PC2",
	                                        "PC1 vs. PC3" = "PC1_PC3",
						"PC2 vs. PC3" = "PC2_PC3")),
	    h5("Colour by"),
	    uiOutput("variable"),
	    actionButton("PCA", "PCA"),

	    h4("Differential expression results"),
	    uiOutput("ma.contrast"),
	    h5("Thresholds"),
	    numericInput("ma.lfc", label="lfc", value = 1),
	    actionButton("MA", "MA plot"),
	    actionButton("heatmap", "Heatmap"),

            h5("Dispay results table"),
	    actionButton("show.results", "Tabulate"),

	    h2("Compare results across datasets/contrasts"),
	    selectInput("dataset1", "dataset 1", choices=getDatasetToContrastNames(conn)),
	    selectInput("dataset2", "dataset 2", choices=getDatasetToContrastNames(conn)),
	    actionButton("scatter.lfc", "Scatterplot lfc"),

            tags$hr(),
            actionButton("reset_input", "Reset inputs")
        ),

    # Main panel for displaying outputs
    mainPanel(dataTableOutput("dataset.table"),
              plotOutput("gene.expression"),
	      dataTableOutput("significant.results"),
	      plotOutput("PCA"),
	      plotOutput("MA"),
	      plotOutput("heatmap"),
              dataTableOutput("tabulate.results"),
	      plotOutput("scatter.lfc")
              )
    )
)



db <- system.file("data/csvdb", package="MIGTranscriptomeExplorer")
conn <- connect(db=db)

# Define server logic
server <- function(input, output) {

    ######################
    # displaying datasets
    ######################

    output$dataset.table <- renderDataTable({
        showDatasets(conn)
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
            expression <- na.omit(expression)
	    metadata <- getMetadata(conn, dataset)
            p <- plotGeneOfInterest(dataset, expression, metadata, variable="treatment")
            grobs.list[[i]] <- p
        }
	# hardcoded
        grid.arrange(grobs=grobs.list, nrow=2, ncol=2)
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
	pcs <- unlist(strsplit(input$PCs, "_"))
        pc <- runPCA(mat)
        plotPrincipleComponents(pc, metadata, colourby=input$variable, pcs=pcs)
    })
    output$PCA <- renderPlot({
        PCA()
    })

    output$ma.contrast <- renderUI({
        contrasts <- getContrasts(conn, input$choose.dataset)
        names(contrasts) <- contrasts
        selectInput("ma.contrast", "contrast", contrasts)
    })
    MA <- eventReactive(input$MA, {
        result <- getResultSet(conn, input$choose.dataset, input$ma.contrast)
        plotMA(result, lfc=input$ma.lfc, title=paste(input$choose.dataset, input$ma.contrast, sep=": "))
    })
    output$MA <- renderPlot({
        MA()
    })

    heatmap <- eventReactive(input$heatmap, {
        mat <- getDiffMatrix(conn, input$choose.dataset, input$ma.contrast, input$ma.lfc, 0.05)
    	heatmapMatrix(mat)
    })
    output$heatmap <- renderPlot({
        heatmap()
    })

    tabulate <- eventReactive(input$show.results, {
        tabulateResults(conn, input$choose.dataset, input$ma.contrast)
    })
    output$tabulate.results <- renderDataTable({
        tabulate()
    })

    scatterlfc <- eventReactive(input$scatter.lfc, {
        df <- buildComparisonSet(conn, input$dataset1, input$dataset2)
	scatterComparisons(df)
    })
    output$scatter.lfc <- renderPlot({
        scatterlfc()
    })

    observeEvent(input$gene, {
        if (is.null(input$gene)){
            shinyjs::diasable("gene.search")}
	else{
	    shinyjs::enable("gene.search")}
    })

    observeEvent(input$reset_input, {
        shinyjs::reset("side-panel")
    })

}

shinyApp(ui, server)