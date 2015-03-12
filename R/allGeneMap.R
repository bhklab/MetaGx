###############
## Natchar Ratanasirigulchai
## March 12, 2015
###############

`allGeneMap ` <-
  function(){
    ## allGeneMap loads the expression sets, and generates a heatmap showing genes present in each dataset
    
    
    library(MetaGx2)
    OvarianEsets <- getOvCaData()
    ## mergedList is all the gene symbols
    mergedList <- fData(OvarianEsets$merged)[,"SYMBOL"]

    # geneMatrix will show 1's where there are genes present in the eset
    geneMatrix <- matrix(data=0, nrow= length(mergedList), ncol= length(OvarianEsets$each))
    dimnames(geneMatrix) <- list(mergedList, c(names(OvarianEsets$each)))


    for (n in 1:length(OvarianEsets$each)){
        print(paste("n ",n))
          genes <- intersect(mergedList,fData(OvarianEsets$each[[n]])[,"SYMBOL"])
          geneMatrix[genes,n] <- 1
      }

    # pdf(file.path(getwd(), "allGeneMap.pdf"), width=5, height=100)
    pdf(file.path(getwd(), "allGeneMap.pdf") )

    heatmap.2(geneMatrix, 
              col=c("#edf8b1", "#2c7fb8"), 
              Rowv = FALSE,
              Colv = FALSE,
              key = FALSE, 
              symkey = FALSE, 
              density.info="none", 
              trace="none",
              labRow = rownames(geneMatrix),
              cexRow = 0.1,
              cexCol = 0.5)
    # heatmap(geneMatrix, Rowv=NA, Colv=NA, col=c("#edf8b1", "#2c7fb8"), scale="none")
    dev.off()
}
