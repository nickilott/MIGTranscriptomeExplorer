# make metadata

makeMetadataMicroarray <- function(matrix.file, output.dataset.name){

    dat <- read.csv(matrix.file, header=T, stringsAsFactors=F, sep="\t")
    samples <- colnames(dat)[1:ncol(dat)-1]
    treatment <- unlist(strsplit(samples, "\\.R[0-9]*"))
    metadata <- data.frame("sample"=samples, "treatment"=treatment)
    write.table(metadata, file=paste0(output.dataset.name, ".metadata"), sep="\t", row.names=F, quote=F)
}

convertlfc <- function(result.file, outfile, platform="microarray"){

    dat <- read.csv(result.file, header=T, stringsAsFactors=F, sep="\t")

    if (paltform == "microarray"){
        dat$logFC <- dat$logFC*-1
        dat$ID <- as.character(dat$ID)
        write.table(dat, file=outfile, sep="\t", row.names=F)}
    else{
        dat$log2FoldChange <- dat$log2FoldChange*-1
        write.table(dat, file=outfile, sep="\t", row.names=F)}
        
}