# plotting functions for MITTranscriptome db


library(ggplot2)
library(reshape)

#' Plot expression values for a gene of interest
#'
#' Return a plot of expression values
#' @param dataset string (dataset to retrieve data from)
#' @param mat expression matrix(subsetted)
#' @param metadata dataframe of metadata
#' @param variable variable present in metadata to colour by
#' @import ggplot2
#' @import reshape
#' @examples
#' @export

plotGOI <- function(dataset, mat, metadata, variable="treatment"){

    # sort the metadata
    sortMetadata(mat, metadata)    

    # add the test_id as variable
    mat$test_id <- rownames(mat)

    # reshape
    mat.m <- melt(mat)

    # add metadata
    mat.m$covariate <- rep(metadata[,variable], nrow(mat))

    colours <- rainbow(length(unique(metadata[,variable])), s=0.7, v=0.6)

    # plotting
    plot1 <- ggplot(mat.m, aes(x="covariate", y=value, colour=covariate))
    plot2 <- plot1 + geom_boxplot(outlier.alpha=0)
    plot3 <- plot2 + geom_jitter(height=0, width=0.15)
    plot4 <- plot3 + theme_bw()
    plot5 <- plot4 + ggtitle(dataset)
    plot6 <- plot5 + facet_wrap(~test_id, nrow=1)
    plot7 <- plot6 + ylab("Expression level")
    plot8 <- plot7 + scale_colour_manual(values=colours) + xlab("")
    return(plot8)
}