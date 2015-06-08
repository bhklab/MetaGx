########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

########################
## Natchar Ratanasirigulchai
## changes made on March 2, 2015
########################


`datasetMerging` <- 
function (esets, method=c("union", "intersect"), standardization=c("quantile", "robust.scaling", "scaling", "none"), nthread=1) {
  ## function merging all individual esets and merging them into a big eset
  #
  # Args:
  #   esets: The list containing all GSE file that need to be merged.      
  #   duplication.checker: A marker, either TRUE or FALSE if you you want to verify
  #                        wheter or not you have duplicate samples into your master 
  #                        gene expression matrix.
  #   survdata: For t.fs and e.fs
  #   time.cens: maximum follow up (years)
  #   Method: either "unique" or "intersect" is use to for selecting geneid
  #
  # Returns:     
  #     The merging eset
  
  method <- match.arg(method)
  standardization <- match.arg(standardization)
  
  ## all unique Entrez gene ids
  ## gene ids
  ugid <- lapply(esets, function(x) { return(Biobase::fData(x)) })
  ugid <- do.call(rbind, ugid)
  ugid <- ugid[!is.na(ugid[ , "EntrezGene.ID"]) & !duplicated(ugid[ , "EntrezGene.ID"]), , drop=FALSE]
  rownames(ugid) <- gsub(sprintf("(%s).", paste(names(esets), collapse="|")), "", rownames(ugid))
  switch (method,
    "union" = {
      feature.merged <- ugid
    },
    "intersect" = {
      feature.merged <- lapply(esets, function(x) { return(stripWhiteSpace(as.character(Biobase::fData(x)[ , "EntrezGene.ID"]))) })
      feature.merged <- table(unlist(feature.merged))
      feature.merged <- names(feature.merged)[feature.merged == length(esets)]
      feature.merged <- ugid[match(feature.merged, stripWhiteSpace(as.character(ugid[ , "EntrezGene.ID"]))), , drop=FALSE]
    },
    {
      stop("Unknown method")
    }
  )
  ## expression data
  exprs.merged <- lapply(esets, function (x, y) {
    ee <- Biobase::exprs(x)[is.element(rownames(exprs(x)),rownames(feature.merged)),]
    #print(dim(ee))
    eem <- matrix(NA, nrow=length(y), ncol=ncol(ee), dimnames=list(y, colnames(ee)))
    #print(dim(eem))
    #print(length(intersect(rownames(ee),rownames(eem))))
    eem[rownames(ee), colnames(ee)] <- ee
    return (eem)
  }, y=rownames(feature.merged))
  exprs.merged <- do.call(cbind, exprs.merged)
  ## clinical info
  ucid <- lapply(esets, function(x) { return(colnames(Biobase::pData(x))) })
  ucid <- table(unlist(ucid))
  ucid <- names(ucid)[ucid == length(esets)]
  clinicinfo.merged <- lapply(esets, function (x , y) {
    ee <- Biobase::pData(x)[ , y, drop=FALSE]
  }, y=ucid)
  clinicinfo.merged <- do.call(gdata::combine, clinicinfo.merged)
  rownames(clinicinfo.merged) <- colnames(exprs.merged)
#   rownames(clinicinfo.merged) <- gsub(   sprintf("(%s).", paste(names(esets), collapse="|")), "", rownames(clinicinfo.merged)         )
#   ## create a merged expressionSet object
  eset.merged <- ExpressionSet(assayData=exprs.merged, phenoData=AnnotatedDataFrame(data=clinicinfo.merged), featureData=AnnotatedDataFrame(data=feature.merged))
  experimentData(eset.merged)@preprocessing <- list("normalization"="mixed", package="unspecified", version="0")
  annotation(eset.merged) <- "mixed"
  ## subtyping
  sbtn <- lapply(esets, function (x) {
    return (colnames(getSubtype(eset=x, method="fuzzy")))
  })
  if (!all(sapply(sbtn, is.null))) {
    sbtn <- table(unlist(sbtn))
    if (!all(sbtn == length(esets))) { stop("Different subtyping across esets") }
    sclass <- lapply(esets, getSubtype, method="class")
#     nn <- unlist(sapply(sclass, names))
#     sclass <- unlist(sclass)
#     names(sclass) <- nn
#     names(sclass) <- names(esets)
    sclass <- unlist(sclass)
    names(sclass) <- unlist(sapply(esets, sampleNames))
    sfuzzy <- do.call(rbind, lapply(esets, getSubtype, method="fuzzy"))
    rownames(sfuzzy) <- sampleNames(eset.merged)
    scrisp <- do.call(rbind, lapply(esets, getSubtype, method="crisp"))
    message("going to set the subtype")
    eset.merged <- setSubtype(eset=eset.merged, subtype.class=sclass, subtype.fuzzy=sfuzzy, subtype.crisp=scrisp)
    rownames(fData(eset.merged)) <- fData(eset.merged)[,"EntrezGene.ID"]
    message("set the subtype complete for merged")
  }
    
  ## standardization
  switch(standardization,
    "none" = {
      ## do nothing
    },
    "quantile" = {
      ## robust scaling followed by quantile normalization
      ee <- exprs(eset.merged)
      # ee <- apply(ee, 2, genefu::rescale)
      splitix <- parallel::splitIndices(nx=ncol(ee), ncl=nthread)
      mcres <- parallel::mclapply(splitix, function(x, data) {
        res <- apply(data[ , x, drop=FALSE], 2, function (dx) {
          return ((genefu::rescale(dx, q=0.05, na.rm=TRUE) - 0.5) * 2)
        })
        return (res)
      }, data=ee, mc.cores=nthread)
      ee <- do.call(cbind, mcres)
      ## quantile normalization
      ee <- limma::normalizeBetweenArrays(object=ee, method="quantile")
      exprs(eset.merged) <- ee
    },
    "robust.scaling" = {
      ## robust scaling
      ee <- exprs(eset.merged)
      # ee <- apply(ee, 2, genefu::rescale)
      splitix <- parallel::splitIndices(nx=ncol(ee), ncl=nthread)
      mcres <- parallel::mclapply(splitix, function(x, data) {
        res <- apply(data[ , x, drop=FALSE], 2, function (dx) {
          return ((genefu::rescale(dx, q=0.05, na.rm=TRUE) - 0.5) * 2)
        })
        return (res)
      }, data=ee, mc.cores=nthread)
      ee <- do.call(cbind, mcres)
      exprs(eset.merged) <- ee
    },
    "scaling" = {
      ## traditional scaling
      # robust scaling
      ee <- exprs(eset.merged)
      # ee <- apply(ee, 2, genefu::rescale)
      splitix <- parallel::splitIndices(nx=ncol(ee), ncl=nthread)
      mcres <- parallel::mclapply(splitix, function(x, data) {
        return(apply(data[ , x, drop=FALSE], 2, scale))
      }, data=ee, mc.cores=nthread)
      ee <- do.call(cbind, mcres)
      exprs(eset.merged) <- ee
    },
    {
      stop("Unknown data standardization method")
    }
  )
    
  return (eset.merged)
}

