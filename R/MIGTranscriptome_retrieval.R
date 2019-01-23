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
    dbDisconnect(connection)
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

#' Sort metadata
#'
#' Sort metadata by sample names given in an expression matrix
#' @param mat expression matrix (data frame)
#' @param metadata metadata returned by getMetadata (data frame)
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

