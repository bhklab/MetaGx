########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

`platformMerging` <- 
function (esets, probe=c("intersect", "union")) {
  # Args:
  #     esets: list of ExpressionSet objects
  #   
  # Returns:     
  #     The eset with all platforms merged
  
  probe <- match.arg(probe)
  
  if (class(esets) == "ExpressionSet") {
    ## a single eset, no platform to merge
    return(esets)
  }
  if (is.list(esets)) {
    esets.check <- sapply(esets, function(x) { return(class(x) == "ExpressionSet") })
    if (any(!esets.check)) { stop("Some esets in the list are not ExpressionSet") }
    pp <- lapply(esets, function (x) { return (rownames(exprs(x))) } )
    switch (probe,
      "intersect" = {
        probex <- intersectList(pp)
      },
      "union" = {
        probex <- sort(unique(unlist(pp)))
      }
    )
    if (length(probex) == 0) { stop("Not enough probes") }
    eset <- esets[[1]]
    ## update expression values
    Biobase::exprs(eset) <- matrix(NA, nrow=length(probex), ncol=nrow(Biobase::pData(eset)), dimnames=list(probex, rownames(Biobase::pData(eset))))
    probex2 <- intersect(probex, rownames(Biobase::fData(eset)))
    Biobase::exprs(eset)[probex2, ] <- Biobase::exprs(esets[[1]])[probex2, , drop=FALSE]
    ## update feature data
    fData(eset) <- data.frame(matrix(NA, nrow=length(probex), ncol=ncol(Biobase::fData(eset)), dimnames=list(probex, colnames(Biobase::fData(eset)))))
    ff <- apply(Biobase::fData(esets[[1]]), 2, function (x) { return(stripWhiteSpace(as.character(x))) })
    dimnames(ff) <- dimnames(Biobase::fData(esets[[1]]))
    Biobase::fData(eset)[probex2, ] <- ff[probex2, , drop=FALSE]
    annotation(eset) <- "merged"
    if (length(esets) > 1) {
      for (j in 2:length(esets)) {
        eset2 <- esets[[j]]
        ## merge expression values
        Biobase::exprs(eset2) <- matrix(NA, nrow=length(probex), ncol=nrow(Biobase::pData(eset2)), dimnames=list(probex, rownames(Biobase::pData(eset2))))
        probex2 <- intersect(probex, rownames(Biobase::fData(eset2)))
        Biobase::exprs(eset2)[probex2, ] <- Biobase::exprs(esets[[j]])[probex2, , drop=FALSE]
        ## update feature data
        fData(eset2) <- data.frame(matrix(NA, nrow=length(probex), ncol=ncol(Biobase::fData(eset2)), dimnames=list(probex, colnames(Biobase::fData(eset2)))))
        ff <- apply(Biobase::fData(esets[[j]]), 2, function (x) { return (stripWhiteSpace(as.character(x))) })
        dimnames(ff) <- dimnames(Biobase::fData(esets[[j]]))
        Biobase::fData(eset2)[probex2, ] <- ff[probex2, , drop=FALSE]
        annotation(eset2) <- "merged"
        eset <- BiocGenerics::combine(eset, eset2)
      }
    }
    return(eset)
  }
}

