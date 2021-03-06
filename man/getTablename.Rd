% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MIGTranscriptome_retrieval.R
\name{getTablename}
\alias{getTablename}
\title{Return tablename based on attributes}
\usage{
getTablename(dataset, type = "matrix")
}
\arguments{
\item{dataset}{dataset to base tablename prefix on}

\item{type}{one of matrix, probe2gene_map, metadata}
}
\value{
string
}
\description{
Return tablename based on attributes
}
\examples{
getTablename("MIGTranscriptome_0001", type="matrix")
}
