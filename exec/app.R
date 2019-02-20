# load packages
library(shiny)
library(MIGTranscriptomeExplorer)
library(gridExtra)
library(RSQLite)
library(dplyr)
library(shinythemes)
library(shinyjs)
library(ggplot2)
library(VennDiagram)

# Define UI for MIGTranscriptomeExplorer app

db <- system.file("data/csvdb", package="MIGTranscriptomeExplorer")
conn <- connect(db=db)

ui <- fluidPage(theme=shinytheme("flatly"),
  
    # App title
    titlePanel("MIGTranscriptomeExplorer"),
  
    # Sidebar panel for inputs
    sidebarLayout(
        sidebarPanel(
	    h2("1.Explore gene expression across datasets"),

	    h4("Search for gene in database"),
	    textInput("gene", label="Gene:", value = ""),
            actionButton("gene.search", "get expression"),

	    h4("Significant contrasts in database"),
	    h5("Specify thresholds"),
	    numericInput("lfc", label="lfc:", value = 0),
	    numericInput("padj", label="padj", value = 0.05),
	    actionButton("significant", "get significant"),

	    h2("2.Explore specific dataset"),
	    h4("Choose dataset"),
	    selectInput("choose.dataset", "Dataset:", choices=showDatasets(conn)$dataset),

	    h4("Principle Components Analysis"),

	    selectInput("PCs", "PC:", choices=c("PC1 vs. PC2" = "PC1_PC2",
	                                        "PC1 vs. PC3" = "PC1_PC3",
						"PC2 vs. PC3" = "PC2_PC3")),
	    h5("Colour by"),
	    uiOutput("variable"),
	    actionButton("PCA", "PCA"),
	    downloadButton("download.pca", "Download"),

	    h4("Differential expression results"),
	    uiOutput("ma.contrast"),
	    h5("Thresholds"),
	    numericInput("ma.lfc", label="lfc", value = 1),
	    actionButton("MA", "MA plot"),
	    downloadButton("download.ma", "Download"),

            h5("Dispay results table"),
	    actionButton("show.results", "Tabulate"),

	    h5("Export results to current directory"),
            downloadButton("download.table", "Download"),

	    h2("3.Compare results across datasets/contrasts"),
	    selectInput("dataset1", "dataset 1", choices=getDatasetToContrastNames(conn)),
	    selectInput("dataset2", "dataset 2", choices=getDatasetToContrastNames(conn)),
	    actionButton("scatter.lfc", "Scatterplot lfc"),
            downloadButton("download.scatter", "Download"),

            h5("Venn diagram"),
	    numericInput("venn.lfc", "lfc", value=1),
	    actionButton("venn", "venn diagram")
        ),

    # Main panel for displaying outputs
    mainPanel(
              tabsetPanel(
	          tabPanel("datasets", dataTableOutput("dataset.table")),
                  tabPanel("1.Expression across datasets", plotOutput("gene.expression", height=800),
	                                                 dataTableOutput("significant.results")),
	          tabPanel("2.Explore dataset", plotOutput("PCA"),
	                                      plotOutput("MA", brush = "plot_brush_ma"),
                                              verbatimTextOutput("gene.info.ma"),
                                              dataTableOutput("tabulate.results")),
                  tabPanel("3.Compare datasets", plotOutput("scatter.lfc", brush = "plot_brush"),
                                               verbatimTextOutput("gene.info"),
	                                       plotOutput("venn"))
                 )
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
        grid.arrange(grobs=grobs.list, nrow=3, ncol=2)
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

    ma.df <- eventReactive(input$MA, {
        getResultSet(conn, input$choose.dataset, input$ma.contrast)
    })

    MA <- eventReactive(input$MA, {
        plotMA(ma.df(), lfc=input$ma.lfc, title=paste(input$choose.dataset, input$ma.contrast, sep=": "))
    })

    output$MA <- renderPlot({
        MA()
    })

    output$gene.info.ma <- renderPrint({
	brushedPoints(na.omit(ma.df()), input$plot_brush_ma)          
    })

    tabulate <- eventReactive(input$show.results, {
        tabulateResults(conn, input$choose.dataset, input$ma.contrast)
    })

    output$tabulate.results <- renderDataTable({
        tabulate()
    })

    df <- eventReactive(input$scatter.lfc, {
        buildComparisonSet(conn, input$dataset1, input$dataset2)
    })

    scatterlfc <- eventReactive(input$scatter.lfc, {
      scatterComparisons(na.omit(df()))
    })

    output$scatter.lfc <- renderPlot({
        scatterlfc()
    })

    output$gene.info <- renderPrint({
	brushedPoints(na.omit(df()), input$plot_brush)          
    })

    df.venn <- eventReactive(input$venn, {
        buildComparisonSet(conn, input$dataset1, input$dataset2)
    })

    venndiagram <- eventReactive(input$venn, {
        vennComparisons(df.venn(), input$venn.lfc)
    })

    output$venn <- renderPlot({
        venndiagram()
    })

    ##############
    # downloads
    ##############

    output$download.pca <- downloadHandler(
	filename = function() {
	    paste0(input$choose.dataset, "_", "pca", ".pdf")
	},
        content = function(file){
	    ggsave(file, plot=PCA())
    })

    output$download.ma <- downloadHandler(
	filename = function() {
	    paste0(input$choose.dataset, "__", input$ma.contrast, "_", "ma", ".pdf")
	},
        content = function(file){
	    ggsave(file, plot=MA())
    })

    output$download.table <- downloadHandler(
	filename = function() {
	    paste0(input$choose.dataset, "__", input$ma.contrast, "_", "result", ".tsv")
	},
        content = function(file){
	    write.table(tabulate(), file, sep="\t", quote=F, row.names=F)
    })
    
    output$download.scatter <- downloadHandler(
	filename = function() {
	    paste0(input$choose.dataset, "_", input$dataset1, "_vs_", input$dataset2, ".pdf")
	},
        content = function(file){
	    ggsave(file, plot=output$scatter.lfc)
   })
}

shinyApp(ui, server)