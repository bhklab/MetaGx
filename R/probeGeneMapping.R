########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################


`probeGeneMapping` <- 
function (eset, platform=c("MISC", "GPL8300", "GPL96", "GPL3921", "GPL97", "GPL570", "GPL571", "GPL1352"), method=c("variance", "jetset")){
  ## probe-gene mapping: which package to use for which platform
  #
  # Args:
  #    eset: ExpressionSet object
  #    platform: identifier of the microarray platform (usually a GPL id from GEO)
  #    method: either the most variant probes or jetset for Affymetrix platform
  #   
  # Returns:     
  #     updated ExpressionSet object with single probe per Entrez gene id

  
  platform <- match.arg(platform)
  method <- match.arg(method)
  
  platf.map <- rbind(c("MISC", "variance", ""),
    c("GPL8300", "jetset", "hgu95av2"),
    c("GPL96", "jetset", "hgu133a"),
    c("GPL3921", "jetset", "hgu133a"),
    c("GPL97", "jetset", "hgu133plus2"),
    c("GPL570", "jetset", "hgu133plus2"),
    c("GPL571", "jetset", "hgu133a"),
    c("GPL1352", "jetset", "u133x3p"))
  dimnames(platf.map) <- list(platf.map[ , 1], c("platform", "method", "parameters"))
  if (!is.element(method, platf.map[platf.map[ , "platform"], "method"])) {
    stop(sprintf("Method %s cannot be applied on platform %s\nUse the following method(s) instead: %s", method, platform, paste(x=platf.map[platf.map[ , "platform"] == platform, "method"], collapse=", ")))
  }
  params <- platf.map[which(platform == platf.map[ , "platform"]), "parameters"]
  
  ## keep only ENTREZID and SYMBOL in feature annotation
  Biobase::fData(eset) <- Biobase::fData(eset)[ , c("ENTREZID", "SYMBOL"), drop=FALSE]
  Biobase::fData(eset)[ , "ENTREZID"] <- stripWhiteSpace(as.character(Biobase::fData(eset)[ , "ENTREZID"]))
  Biobase::fData(eset)[ , "SYMBOL"] <- stripWhiteSpace(as.character(Biobase::fData(eset)[ , "SYMBOL"]))
  
  switch (method,
    "jetset" = {
      js <- jetset::jscores(chip=params, probeset=rownames(Biobase::exprs(eset)))
      js <- js[rownames(Biobase::exprs(eset)), , drop=FALSE]
      ## identify the best probeset for each Entrez Gene ID
      geneid1 <- stripWhiteSpace(as.character(js[ ,"EntrezID"]))
      names(geneid1) <- rownames(js)
      geneid2 <- sort(unique(geneid1))
      names(geneid2) <- paste("geneid", geneid2, sep=".")
      gix1 <- !is.na(geneid1)
      gix2 <- !is.na(geneid2)
      geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
      ## probes corresponding to common gene ids
      gg <- names(geneid1)[is.element(geneid1, geneid.common)]
      gid <- geneid1[is.element(geneid1, geneid.common)]
      ## duplicated gene ids
      gid.dupl <- unique(gid[duplicated(gid)])
      gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
      ## unique gene ids
      gid.uniq <- gid[!is.element(gid, gid.dupl)]
      gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
      ## which are the best probe for each gene
      js <- data.frame(js, "best"=FALSE, stringsAsFactors=FALSE)
      js[gg.uniq, "best"] <- TRUE
      ## data for duplicated gene ids
      if(length(gid.dupl) > 0) {	
      	## use jetset oevrall score to select the best probesets
      	myscore <- js[gg.dupl,"overall"]
      	myscore <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "score"=myscore)
      	myscore <- myscore[order(as.numeric(myscore[ , "score"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
      	myscore <- myscore[!duplicated(myscore[ , "gid"]), , drop=FALSE]
      	js[myscore[ ,"probe"], "best"] <- TRUE
      }
      ## update the esets
      probes <- rownames(Biobase::exprs(eset))[js[ , "best"]]
      names(probes) <- paste("geneid", js[js[ , "best"], "EntrezID"], sep=".")
      gid <- js[js[ , "best"], "EntrezID"]
      gsymb <- js[js[ , "best"], "symbol"]
      Biobase::exprs(eset) <- Biobase::exprs(eset)[probes, , drop=FALSE]
      rownames(Biobase::exprs(eset)) <- names(probes)
      Biobase::fData(eset) <- Biobase::fData(eset)[probes, , drop=FALSE]
      rownames(Biobase::fData(eset)) <- names(probes)
      Biobase::fData(eset)[ , "ENTREZID"] <- gid
      Biobase::fData(eset)[ , "SYMBOL"] <- gsymb
      Biobase::fData(eset) <- cbind(Biobase::fData(eset), "PROBEID"=probes, stringsAsFactors=FALSE)
    },
    "variance" = {
      ## other platform, select the most variant probe per Entrez Gene ID
      gid <- stripWhiteSpace(as.character(Biobase::fData(eset)[ , "ENTREZID"]))
      names(gid) <- rownames(Biobase::exprs(eset))
      ugid <- sort(unique(gid))
      rr <- genefu::geneid.map(geneid1=gid, data1=t(Biobase::exprs(eset)), geneid2=ugid)
      probes <- colnames(rr$data1)
      names(probes) <- paste("geneid", rr$geneid1, sep=".")
      Biobase::exprs(eset) <- Biobase::exprs(eset)[probes, , drop=FALSE]
      rownames(Biobase::exprs(eset)) <- names(probes)
      Biobase::fData(eset) <- Biobase::fData(eset)[probes, , drop=FALSE]
      rownames(Biobase::fData(eset)) <- names(probes)
      Biobase::fData(eset) <- cbind(Biobase::fData(eset), "PROBEID"=probes, stringsAsFactors=FALSE)
      ## get the gene symbols from entrez gene id using org.Hs.eg.db
      gs <- toTable(org.Hs.egSYMBOL)
      gs <- gs[!duplicated(gs[ , "gene_id"]), , drop=FALSE]
      rownames(gs) <- gs[ , "gene_id"]
      gs <- gs[stripWhiteSpace(as.character(Biobase::fData(eset)[ , "ENTREZID"])), "symbol"]
      Biobase::fData(eset)[ , "SYMBOL"] <- stripWhiteSpace(as.character(gs))
    },
    {
      stop(sprintf("Unknown method for probe-gene mapping for platform %s", platform))
    }
  )
  return (eset)
}

