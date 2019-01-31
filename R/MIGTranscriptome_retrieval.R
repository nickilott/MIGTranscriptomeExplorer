# Data accessing from within R using RSQLite3

library(RSQLite)

##################################################
##################################################
##################################################

#' Connect to the database
#'
#' Connects to the SQLite database
#' @param db database to connect to (defaults to csvdb)
#' @return an sqlite connection
#' @import RSQLite 
#' @examples
#' connect(db="csvdb")
#' @export

connect <- function(db="csvdb"){

    # connect to a database
    sqlite <- dbDriver("SQLite")
    con <- dbConnect(sqlite, db)
    return(con)
}

##################################################
##################################################
##################################################

#' Show datasets and their attributes
#'
#' Display datasets and attributes related to those
#' datasets including metadata and contacts
#' @param connection RSQLite database connection
#' @return data frame of datasets and their attributes
#' @examples
#' showDatasets(connect("csvdb"))
#' @export

showDatasets <- function(connection){

    # return data frame of dataset information
    statement = 'SELECT * FROM reference'
    datasets <- dbGetQuery(connection, statement) 
    return(datasets)
}

##################################################
##################################################
##################################################

#' Return tablename based on attributes
#'
#' Return tablename based on attributes
#' @param dataset dataset to base tablename prefix on
#' @param type one of matrix, probe2gene_map, metadata
#' @return string
#' @export
#' @examples
#' getTablename("MIGTranscriptome_0001", type="matrix")

getTablename <- function(dataset, type="matrix"){

    # return tablename as string depending on type
    tablename = paste(dataset, type, sep="_")
    return(tablename)
}

##################################################
##################################################
##################################################

#' Return matrix of expression values for a given dataset
#'
#' Return a normalised expression matrix for a given
#' dataset in the database
#' @param connection sqlite connection
#' @param dataset string (dataset to retrieve matrix for)
#' @return data frame
#' @examples
#' getMatrix(connect(), "MIGTranscriptome_0001")
#' @export

getMatrix <- function(connection, dataset){

    # return dataframe of normalised counts
    tablename <- getTablename(dataset, type="matrix")
    statement <- paste0('SELECT * from ', tablename)
    mat <- dbGetQuery(connection, statement)
    rownames(mat) <- mat$test_id
    mat <- mat[,c(2:ncol(mat))]
    return(mat)
}

##################################################
##################################################
##################################################

#' Convert gene name to uppercase
#'
#' Convert gen name to uppercase
#' @param gene string (gene name in upper or lower case)
#' @return string
#' @examples
#' convertGene("S100a9") # returns S100A9

convertGene <- function(gene){

    # return upper case version of input
    gene <- toupper(gene)
    gene <- paste0('"', gene, '"')
    return(gene)
}

##################################################
##################################################
##################################################

#' Return probe/ensembl id set 
#'
#' Get all probes/ensembl ids that match a given gene
#' from the database
#' @param connection sqlite connection
#' @param dataset string (dataset to retrieve ids from)
#' @param gene gene to search for matching probes/ensembl ids
#' @return vector of probes/ensembl gene ids
#' @examples
#' getProbes(connect(), "MIGTranscriptome_0001", "S100a9")

getProbes <- function(connection, dataset, gene){

    # return vector of probes/ensembl ids for a specified
    # gene
    tablename <- getTablename(dataset, type="probe2gene_map")
    statement <- paste0('SELECT probe FROM ', tablename, ' WHERE gene_name==', convertGene(gene))
    probes <- dbGetQuery(connection, statement)
    return(as.character(probes$probe))
}

##################################################
##################################################
##################################################

#' Return expression matrix for a given gene
#'
#' Return an expression matrix that contains values
#' for each probe/ensembl id for the given gene
#' @param connection sqlite connection
#' @param dataset string (dataset to retrieve ids from)
#' @param gene gene to get expression values for
#' @return data frame of expression values
#' @examples
#' getExpression(connect(), "MIGTranscriptome_0001", "S100a9")
#' @export

getExpression <- function(connection, dataset, gene){

    # return expression matrix for probes/ensembl ids
    # matching gene
    probes <- getProbes(connection, dataset, gene)
    mat <- getMatrix(connection, dataset)
    mat <- mat[probes,]
    return(mat)
}

##################################################
##################################################
##################################################

#' Get metadata
#'
#' Return the metadata for a given dataset
#' @param connection sqlite connection
#' @param dataset string (dataset to retrieve metadata for)
#' @examples
#' getMetadata(connect(), "MIGTranscriptome_0001")
#' @export

getMetadata <- function(connection, dataset){

    # return data frame of metadata
    tablename <- getTablename(dataset, type="metadata")
    statement <- paste0('SELECT * FROM ', tablename)
    metadata <- dbGetQuery(connection, statement)
    rownames(metadata) <- metadata$sample
    return(metadata)
}

##################################################
##################################################
##################################################

#' Get list of variables
#'
#' Return a list of metadata variables for a given dataset
#' @param connection sqlite connection
#' @param dataset string (dataset to retrieve metadata for)
#' @examples
#' getMetadataList(connect(), "MIGTranscriptome_0001")
#' @export

getMetadataList <- function(connection, dataset){

    # return a list of variables from metadata
    statement <- paste0('SELECT metadata FROM reference WHERE dataset==', '"', dataset, '"')
    metadata <- dbGetQuery(connection, statement)$metadata
    metadata <- unlist(strsplit(metadata, ","))
    return(metadata)
}

##################################################
##################################################
##################################################

#' Sort metadata
#'
#' Sort metadata by sample names given in an expression matrix
#' @param mat expression matrix (data frame)
#' @param metadata metadata returned by getMetadata (data frame)
#' @export
#' @examples
#' mat <- getMatrix(connect(), "MIGTranscriptome_0001")
#' metdata <- getMetadata(connect(), "MIGTranscriptome_0001")
#' sortMetadata(mat, metadata)

sortMetadata <- function(mat, metadata){

    # sort metadata according to columns in
    # matrix
    metadata <- metadata[colnames(mat),]
    return(metadata)
}


##################################################
##################################################
##################################################

#' Get contrasts
#'
#' Get contrasts for a specified dataset
#' @param connection database connection
#' @param dataset dataset of interest in database
#' @export
#' @examples
#' getContrasts(conn, "MIGTranscriptome_0001")

getContrasts <- function(conn, dataset){

    statement <- paste0('SELECT contrasts FROM reference WHERE dataset==', '"', dataset, '"')
    contrasts <- dbGetQuery(conn, statement)
    contrasts <- gsub("-", "_", contrasts$contrasts)
    contrasts <- unlist(strsplit(contrasts, ","))
    return(contrasts)
}

##################################################
##################################################
##################################################

#' Get significant sets
#'
#' Get significant expression differences for dataset and contrast
#' of interest
#' @param connection database connection
#' @param dataset dataset of interest in database
#' @param contrast contrast of interest
#' @param lfc log2 fold change threshold for defining significance
#' @param padj adjusted p-value threshold for defining significance
#' @param gene gene to find out if significant
#' @export
#' @examples
#' getSignificant(conn, "MIGTranscriptome_0001", "SNHh_TSB", "IL10")

getSignificant <- function(connection, dataset, contrast, lfc, padj, gene){

    # input comes as characters
    # so need to change here
    lfc <- as.numeric(lfc)
    padj <- as.numeric(padj)

    probes <- getProbes(connection, dataset, gene)
    if (length(probes) == 1){
        probes <- convertGene(probes)
    }
    
    probes <- paste0(probes, collapse=",")
    tablename <- getTablename(getTablename(dataset, contrast), type="result")

    # query
    statement <- paste0('SELECT * FROM ', tablename, ' WHERE test_id IN ', '(', probes, ')', ' AND padj < ', padj, ' AND ABS(l2fold) > ', lfc)

    significant <- dbGetQuery(connection, statement)
    if (nrow(significant) == 0){
        significant <- NA
    }
    else
    {
    significant$dataset <- dataset
    significant <- significant[, c("dataset", "test_id", "l2fold", "padj")]
    }
    return(significant)		
}

##################################################
##################################################
##################################################

#' Get differential expression results
#'
#' Get the differential expression results for
#' a given dataset and contrast
#' @param conn connection
#' @param dataset dataset to choose
#' @param contrast contrast of interest
#' @export
#' @examples
#' getResultSet(connect(), "MIGTranscriptome_0001", "SNHh_TSB")

getResultSet <- function(conn, dataset, contrast){

    tablename <- getTablename(dataset, type="result")
    tablename <- paste0(dataset, "_", contrast, "_result")
    statement <- paste0('SELECT * FROM ', '"', tablename, '"')
    dat <- dbGetQuery(conn, statement)
    probe2gene <- getTablename(conn, input$dataset, type="probe2gene_map")
    rownames(probe2gene$probe)
    dat$gene_name <- probe2gene[dat$test_id,]$gene_name
    return(dat)
}

##################################################
##################################################
##################################################

#' Run PCA
#'
#' Run PCA on matrix
#' @param df data frame of normalised counts
#' @export
#' @examples
#' runPCA(df)

runPCA <- function(df){

    pc <- prcomp(t(df))
    return (pc)
}

##################################################
##################################################
##################################################

#' Get PCA
#'
#' Get princple components
#' @param pc prcomp object
#' @export
#' @examples
#' getPCA(runPCA(df))

getPCA <- function(pc){

    pcs <- data.frame(pc$x)
    return(pcs)
}

##################################################
##################################################
##################################################

#' Get variance explained
#'
#' Get variance explained for prcomp object
#' @param pc prcomp object
#' @param component string (component to get VE for)
#' @export
#' @examples
#' getVE(runPCA(df))

getVE <- function(pc, component="PC1"){

    pve <- summary(pc)$importance[,component][[2]]
    return (pve)
}



##################################################
##################################################
##################################################

#' Tabulate results

#' Tabulate probes/genes sorted by padj
#' @param connection connection to database
#' @param dataset dataset to look at
#' @param contrast contrast of interest
#' @export
#' @examples
#' tabulateResults(connect(), "MIGTranscriptome_0001", "SNHh_TSB")

tabulateResults <- function(connection, dataset, contrast){

    result.tablename <- getTablename(dataset, type="result")
    result.tablename <- paste0(dataset, "_", contrast, "_result")
    probe2gene <- getTablename(dataset, type="probe2gene_map")
    statement <- paste0('SELECT a.*, b.gene_name FROM ', result.tablename, ' as a, ', probe2gene, ' as b ', ' WHERE a.test_id==b.probe ')
    result <- dbGetQuery(connection, statement)
    res <- result[order(result$padj, decreasing=FALSE),]
    return(res)
}

##################################################
##################################################
##################################################

#' Get expression values for differentially expressed probes/genes

#' Get expression matrix for differentially expressed probes/genes
#' @param connection connection to database
#' @param dataset dataset to see
#' @param contrast contrast of interest
#' @param lfc log fold change threshold
#' @param padj adjusted p-value threshold
#' @export
#' @examples
#' getDiffMatrix(connect(), "MIGTranscriptome_0001", "SNHh_TSB", 1, 0.05)

getDiffMatrix <- function(connection, dataset, contrast, lfc, padj){

    result.tablename <- getTablename(dataset, type="result")
    result.tablename <- paste0(dataset, "_", contrast, "_result")

    matrix.tablename <- getTablename(dataset, type="matrix")

    statement <- paste0('SELECT a.* FROM ', matrix.tablename, ' as a, ', result.tablename,  ' as b  WHERE a.test_id==b.test_id AND b.padj < ', padj, ' AND abs(b.l2fold) > ', lfc)
    mat <- dbGetQuery(connection, statement)
    rownames(mat) <- mat$test_id
    mat <- mat[,2:ncol(mat)]
    return(mat)
}

##################################################
##################################################
##################################################

#' Get choices for comparing datasets

#' These form the choices for comparisons
#' @param connection connection to database
#' @export
#' @examples
#' getDatasetToContrastNames(connect())

getDatasetToContrastNames <- function(connection){

    statement <- 'SELECT dataset, contrasts FROM reference'
    datasets <- dbGetQuery(connection, statement)
    comparison.choices <- c()
    for (dataset in datasets$dataset){
        contrasts <- getContrasts(connection, dataset)
	for (contrast in contrasts){
	    newname <- paste0(dataset, "__", contrast)
	    comparison.choices <- append(newname, comparison.choices)
        }
    }
    return(comparison.choices)
}

##################################################
##################################################
##################################################

#' Build dataframe for comparing fold changes

#' Build data frame for comparing fold changes
#' @param connection connection to database
#' @param dataset1 string (dataset1)
#' @param dataset2 string (dataset2)
#' @export
#' @examples
#' To follow

buildComparisonSet <- function(connection, dataset1, dataset2){


    d1 <- unlist(strsplit(dataset1, "__"))
    d2 <- unlist(strsplit(dataset2, "__"))

    data1 <- d1[1]
    contrast1 <- d1[2]
    data2 <- d2[1]
    contrast2 <- d2[2]

    tablename1 <- getTablename(data1, type=contrast1)
    tablename1 <- getTablename(tablename1, type="result")
    probe2gene1 <- getTablename(data1, type="probe2gene_map")

    tablename2 <- getTablename(data2, type=contrast2)
    tablename2 <- getTablename(tablename2, type="result")
    probe2gene2 <- getTablename(data2, type="probe2gene_map")

    statement1 <- paste0('SELECT a.gene_name, AVG(b.l2fold) as l2fold, b.padj as padj FROM ', probe2gene1,  ' as a, ',  tablename1, ' as b ', 'WHERE a.probe==b.test_id GROUP BY a.gene_name')
    l2fold1 <- dbGetQuery(connection, statement1)

    statement2 <- paste0('SELECT a.gene_name, AVG(b.l2fold) as l2fold, b.padj as padj FROM ', probe2gene2,  ' as a, ',  tablename2, ' as b ', 'WHERE a.probe==b.test_id GROUP BY a.gene_name')
    l2fold2 <- dbGetQuery(connection, statement2)

    df <- merge(l2fold1, l2fold2, by.x="gene_name", by.y="gene_name")

    colnames(df) <- c("gene_name", paste0(dataset1, "_l2fold"), paste0(dataset2, "_l2fold"), paste0(dataset1, "_padj"), paste0(dataset2, "_padj"))
    return(df)

}    

##################################################
##################################################
##################################################
