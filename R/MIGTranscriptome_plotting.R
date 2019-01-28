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

plotGeneOfInterest <- function(dataset, mat, metadata, variable="treatment"){

    # plot for if there's no data
    if (nrow(mat) == 0){
        df <- data.frame(x=1, y=1, text.output="No data")
        plot1 <- ggplot(df, aes(x=x, y=y)) + geom_text(label=text.output)
	plot2 <- plot1 + theme(panel.background=element_rect(fill="white", colour="white"))
	plot3 <- plot2 + theme(panel.grid.major=element_line(colour="white"))
        plot4 <- plot3 + theme(panel.grid.minor=element_line(colour="white"))
        plot5 <- plot4 + theme(legend.position="none")
        plot6 <- plot5 + theme(axis.line=element_blank())
        plot7 <- plot6 + theme(axis.text.x=element_blank())
        plot8 <- plot7 + theme(axis.text.y=element_blank())
        plot9 <- plot8 + theme(axis.ticks=element_blank())
        plot10 <- plot9 + theme(axis.title.x=element_blank())
        plot10 <- plot10 + theme(axis.title.y=element_blank()) + ggtitle(dataset)
        return(plot10)
    }else{

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
        plot1 <- ggplot(mat.m, aes(x=factor(covariate, levels=mixedsort(unique(mat.m$covariate))), y=value, colour=covariate))
        plot2 <- plot1 + geom_boxplot(outlier.alpha=0)
        plot3 <- plot2 + geom_jitter(height=0, width=0.15)
        plot4 <- plot3 + theme_bw()
        plot5 <- plot4 + ggtitle(dataset)
        plot6 <- plot5 + facet_wrap(~test_id, nrow=1)
        plot7 <- plot6 + ylab("Expression level")
        plot8 <- plot7 + scale_colour_manual(values=colours) + xlab("")
        return(plot8)
    }
}

##################################################
##################################################
##################################################

#' Plot PCA
#'
#' Plot PCA of dataset of interest
#' @param pc prcomp object
#' @param colourby vector of conditions
#' @param pcs vector of pcs to plot
#' @export
#' @examples
#' plotPCA(PCA(mat))

plotPrincipleComponents <- function(pc, metadata, colourby="none", pcs=c("PC1", "PC2")){

    # covariate must be in same order as pc rownames

    # get variance explained for each component
    ve1 <- getVE(pc, component=pcs[1])
    ve2 <- getVE(pc, component=pcs[2])

    ve1 <- round(ve1, 2)*100
    ve2 <- round(ve2, 2)*100

    # get data frame of components
    pca <- data.frame(pc$x)

    # add conditions
    pca$condition <- metadata[,colourby]
    pca$condition <- factor(pca$condition, levels=unique(pca$condition))

    # plot
    pc1 <- pcs[1]
    pc2 <- pcs[2]

    # labels
    xlabel <- paste(pc1, ve1, sep=" (")
    xlabel <- paste(xlabel, "%", sep="")
    xlabel <- paste(xlabel, ")", sep="")
    ylabel <- paste(pc2, ve2, sep=" (")
    ylabel <- paste(ylabel, "%", sep="")	
    ylabel <- paste(ylabel, ")", sep="")

    n <- length(unique(pca$condition))
    colours <- rainbow(n, s=0.7, v=0.6)
    
    plot1 <- ggplot(pca, aes_string(x=pc1, y=pc2, group="condition", colour="condition"))
    plot2 <- plot1 + geom_point(size=3)
    plot3 <- plot2 + theme_bw() 
    plot4 <- plot3 + xlab(xlabel) + ylab(ylabel) + scale_colour_manual(values=colours)
    return(plot4) 
}

##################################################
##################################################
##################################################

#' Plot MA
#'
#' Plot expression level against fold change
#' @param dat data frame of results
#' @param lfc log2 fold change threshold
#' @param title title for the plot
#' @export
#' @examples
#' plotMA(dat)

plotMA <- function(dat, lfc=1, title="default title"){

       nup <- nrow(dat[dat$padj < 0.05 & dat$l2fold > lfc & !(is.na(dat$padj)) & !(is.na(dat$l2fold)),])
       ndown <- nrow(dat[dat$padj < 0.05 & dat$l2fold < (-lfc) & !(is.na(dat$padj)) & !(is.na(dat$l2fold)),])

       nup <- paste("Upregulated = ", nup, sep="")
       ndown <- paste("Downregulated = ", ndown, sep="")

       dat$significant <- ifelse(dat$padj < 0.05 & abs(dat$l2fold) > lfc & !(is.na(dat$padj)) & !(is.na(dat$l2fold)), "Yes", "No")

       plot1 <- ggplot(na.omit(dat), aes(x=aveexprs, y=l2fold, colour=significant))
       plot2 <- plot1 + geom_point(pch=18, size=1)
       plot3 <- plot2 + scale_colour_manual(values=c("grey", "darkRed"))
       plot4 <- plot3 + theme_bw()
       plot5 <- plot4 + theme(text=element_text(size=10))
       plot6 <- plot5 + geom_hline(yintercept=c(-(lfc),lfc,0), linetype="dashed")
       plot7 <- plot6 + xlab("Mean expression") + ylab("Log2 fold change")
       plot8 <- plot7 + annotate("text", x=9, y=6, label=nup, size=3)
       plot9 <- plot8 + annotate("text", x=9, y=-6, label=ndown, size=3)
       plot10 <- plot9 + ggtitle(title)
       return (plot10)
       }


##################################################
##################################################
##################################################

#' Heatmap
#'
#' Heatmap of differential expression results
#' @param mat expression matrix
#' @param distfun distance function (dist)
#' @param clustfun cluster function (hclust)
#' @import gplots
#' @import gtools
#' @export
#' @examples
#' heatmap(mat)

heatmapMatrix <- function(mat, distfun="euclidean", clustfun="ward.D2"){

    distf <- function(x) dist(x, method=distfun)
    clustf <- function(x) hclust(x, method=clustfun)

    colours <- colorRampPalette(c("blue", "white", "red"))(75)
    mat.s <- data.frame(t(apply(mat, 1, scale)))
    rownames(mat.s) <- rownames(mat)
    colnames(mat.s) <- colnames(mat)
    mat.s <- na.omit(mat.s)
    mat.s <- mat.s[,mixedsort(colnames(mat.s))]
    heatmap.2(as.matrix(mat.s),
              labRow=NA,
              scale="none",
              col=colours,
	      Rowv=T,
	      Colv=T,
	      trace="none",
	      margins=c(15,15),
	      hclustfun=clustf,
	      distfun=distf)
}

##################################################
##################################################
##################################################

#' Scatterplot of fold changes
#'
#' Scatterplot fold changes between datasets/contrasts
#' @param df data frame 
#' @export
#' @examples
#' scatterComparisons(df)

scatterComparisons <- function(df){

    x <- df[,2]
    y <- df[,3]

    namex <- colnames(df)[2]
    namey <- colnames(df)[3]

    df$d <- densCols(x, y, colramp = colorRampPalette(rev(rainbow(10, v = 0.5, s = 0.5, end = 4/6))))
    p <- ggplot(df) +
    geom_point(aes_string(namex, namey, col = "d"), size = 1) +
    scale_color_identity() + theme_bw()
    p <- p + geom_hline(yintercept=c(-1, 1), colour="darkGrey", linetype="dashed")
    p <- p + geom_vline(xintercept=c(-1, 1), colour="darkGrey", linetype="dashed")
    
    return(p)
}