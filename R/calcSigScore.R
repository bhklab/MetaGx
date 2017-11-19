#' Function to determine the risk prediction scores of patients
#'
#' This function uses the gene expression values for the genes of interest along with the direction of association for the genes of interest supplied by the user to determine a score for the patients. Presumably the genes are believed to be 
#' associated with an event, such as survival, and the scores are used to determine the likelihood of the event for a patient.
#' @param x Matrix containing the gene(s) in the gene list in rows and at least three columns: "probe", "EntrezGene.ID" and "coefficient" standing for the name of the probe, the NCBI Entrez Gene id and the coefficient giving the direction and the strength of the association of each gene in the gene list.
#' @param data Matrix of gene expressions with samples in rows and probes in columns, dimnames being properly defined.
#' @param annot Matrix of annotations with at least one column named "EntrezGene.ID" (for ssp, scm, AIMS, and claudinLow models) or "Gene.Symbol" (for the intClust model), dimnames being properly defined.
#' @param doMapping TRUE if the mapping through Entrez Gene ids must be performed (in case of ambiguities, the most variant probe is kept for each gene), FALSE otherwise.
#' @param mapping Matrix with columns "EntrezGene.ID" and "probe" used to force the mapping such that the probes are not selected based on their variance
#' @param size Integer specifying the number of probes to be considered in  signature computation. The probes will be sorted by absolute value of coefficients.
#' @param cutoff Only the probes with coefficient greater than cutoff will be considered in signature computation.
#' @param signed TRUE if only the sign of the coefficient must be considered in signature computation, FALSE otherwise.
#' @param verbose TRUE to print informative messages, FALSE otherwise.
#' @return a list with 3 elements. "score", which has the risk predictions/scores for each patient, "mapping", which contains info regarding the number of genes in x that were mapped to the genes in data, and "probe", which has the IDs of the probes that were successfully mapped.
#' @export
#' @examples
#'
#' eset = loadMetaData(cancerType = "ovarian", survivalMetric = "overall")[1]
#' esetData = t(eset$E.MTAB.386@assayData$exprs)
#' colnames(esetData) = eset$E.MTAB.386@featureData@data$EntrezGene.ID
#' geneEntrezIds = eset$E.MTAB.386@featureData@data$EntrezGene.ID[1:5]
#' probes = eset$E.MTAB.386@featureData@data$probeset[1:5]
#' geneDirecs = rep(1, length(geneEntrezIds))
#' sig = cbind("probe"=probes, "EntrezGene.ID"=geneEntrezIds, "coefficient"= geneDirecs)
#' rownames(sig) = probes
#' scoreList = calcSigScore(x=sig, esetData, annot=NULL, doMapping=FALSE, signed=TRUE)
#'

calcSigScore = function(x, data, annot, doMapping=FALSE, mapping = NULL, size=0, cutoff=NA, signed=TRUE, verbose=FALSE)
{
  `geneid.map` <-
    function(geneid1, data1, geneid2, data2, verbose=FALSE) {

      nn <- names(geneid1)
      geneid1 <- as.character(geneid1)
      names(geneid1) <- nn
      nn <- names(geneid2)
      geneid2 <- as.character(geneid2)
      names(geneid2) <- nn
      if(is.null(names(geneid1))) { names(geneid1) <- dimnames(data1)[[2]] }
      if(!missing(data2) && is.null(names(geneid2))) { names(geneid2) <- dimnames(data2)[[2]] }
      if(!missing(data1) && !missing(geneid1) && !missing(geneid2)) {
        ## remove probes without any measurements
        na.ix <- apply(data1, 2, function(x) { return(all(is.na(x))) })
        data1 <- data1[ , !na.ix, drop=FALSE]
        geneid1 <- geneid1[!na.ix]
      } else { stop("data1, geneid1 and geneid2 parameters are mandatory!") }
      if(!missing(data2)) {
        ## remove probes without any measurements
        na.ix <- apply(data2, 2, function(x) { return(all(is.na(x))) })
        data2 <- data2[ , !na.ix, drop=FALSE]
        geneid2 <- geneid2[!na.ix]
      } else { data2 <- NULL }

      gix1 <- !is.na(geneid1)
      gix2 <- !is.na(geneid2)

      geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
      if(length(geneid.common) == 0) {
        warning("no gene ids in common!")
        return(list("geneid1"=NA, "data1"=NA, "geneid2"=NA, "data2"=NA))
      }

      ## dataset1
      ## probes corresponding to common gene ids
      gg <- names(geneid1)[is.element(geneid1, geneid.common)]
      gid <- geneid1[is.element(geneid1, geneid.common)]
      ## duplicated gene ids
      gid.dupl <- unique(gid[duplicated(gid)])
      gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
      ## unique gene ids
      gid.uniq <- gid[!is.element(gid, gid.dupl)]
      gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
      ## data corresponding to unique gene ids
      datat <- data1[ ,gg.uniq,drop=FALSE]
      ## data for duplicated gene ids
      if(length(gid.dupl) > 0) {
        if(verbose) { message("\ndataset1 duplicates...") }
        ## compute the standard deviation with a penalization on the number of missing values
        ## this should avoid selecting the most variant probe with a lot of missing values
        pena <- apply(X=data1[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=function(x) { return(sum(is.na(x))) })
        pena <- log((nrow(data1) + 1) / (pena + 1)) + 1
        #pena <- 1
        sdr <- drop(apply(X=data1[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=sd, na.rm=TRUE)) * pena
        mysd <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "sd"=sdr)
        mysd <- mysd[order(as.numeric(mysd[ , "sd"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
        mysd <- mysd[!duplicated(mysd[ , "gid"]), , drop=FALSE]
        datat <- cbind(datat, data1[ , mysd[ , "probe"], drop=FALSE])
      }
      data1 <- datat
      geneid1 <- geneid1[dimnames(data1)[[2]]]

      #dataset2
      if(is.null(data2)) {
        #keep arbitrarily the first occurence of each duplicated geneid
        geneid2 <- geneid2[!duplicated(geneid2) & is.element(geneid2, geneid.common)]
      }
      else {
        ## probes corresponding to common gene ids
        gg <- names(geneid2)[is.element(geneid2, geneid.common)]
        gid <- geneid2[is.element(geneid2, geneid.common)]
        ## duplicated gene ids
        gid.dupl <- unique(gid[duplicated(gid)])
        gg.dupl <- names(geneid2)[is.element(geneid2, gid.dupl)]
        ## unique gene ids
        gid.uniq <- gid[!is.element(gid, gid.dupl)]
        gg.uniq <- names(geneid2)[is.element(geneid2, gid.uniq)]
        ## data corresponding to unique gene ids
        datat <- data2[ ,gg.uniq,drop=FALSE]
        ## data for duplicated gene ids
        if(length(gid.dupl) > 0) {
          if(verbose) { message("\ndataset2 duplicates...") }
          ## compute the standard deviation with a penalization on the number of missing values
          ## this should avoid selecting the most variant probe with a lotof missing values
          pena <- apply(X=data2[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=function(x) { return(sum(is.na(x))) })
          pena <- log((nrow(data2) + 1) / (pena + 1)) + 1
          #pena <- 1
          sdr <- drop(apply(X=data2[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=sd, na.rm=TRUE)) * pena
          mysd <- cbind("probe"=gg.dupl, "gid"=geneid2[gg.dupl], "sd"=sdr)
          mysd <- mysd[order(as.numeric(mysd[ , "sd"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
          mysd <- mysd[!duplicated(mysd[ , "gid"]), , drop=FALSE]
          datat <- cbind(datat, data2[ , mysd[ , "probe"], drop=FALSE])
        }
        data2 <- datat
        geneid2 <- geneid2[dimnames(data2)[[2]]]
      }

      #same order for the two datasets
      rix <- match(geneid2, geneid1)
      geneid1 <- geneid1[rix]
      data1 <- data1[ ,rix,drop=FALSE]
      return(list("geneid1"=geneid1, "data1"=data1, "geneid2"=geneid2, "data2"=data2))
    }


    if(missing(data) || missing(annot)) { stop("data and annot parameters must be specified") }
    x <- as.data.frame(x, stringsAsFactors=FALSE)
    if(nrow(x) == 0) { stop("empty gene list!"); }

    myprobe <- as.character(x[ ,"probe"])
    mygid <- as.character(x[ ,"EntrezGene.ID"])
    mycoef <- as.numeric(x[ ,"coefficient"])
    names(mycoef) <- names(mygid) <- names(myprobe) <- myprobe

    nix <- order(abs(mycoef), decreasing=TRUE, na.last=NA)
    myprobe <- myprobe[nix]
    mygid <- mygid[nix]
    mycoef <- mycoef[nix]

    if(doMapping) { ## mapping is requested
      gid1 <- mygid
      gid2 <- as.character(annot[ ,"EntrezGene.ID"])
      names(gid2) <- dimnames(annot)[[1]]
      ## remove missing and duplicated geneids from the gene list
      rm.ix <- is.na(gid1) | duplicated(gid1)
      gid1 <- gid1[!rm.ix]

      rr <- geneid.map(geneid1=gid2, data1=data, geneid2=gid1, verbose=FALSE)
      if(is.na(rr$geneid1[1])) {
        #no gene ids in common
        res <- rep(NA, nrow(data))
        names(res) <- dimnames(data)[[1]]
        gf <- c("mapped"=0, "total"=nrow(x))
        if(verbose) { message(sprintf("probe candidates: 0/%i", nrow(x))) }
        return(list("score"=res, "mapping"=gf, "probe"=cbind("probe"=NA, "EntrezGene.ID"=NA, "new.probe"=NA)))
      }
      nix <- match(rr$geneid2, mygid)
      myprobe <- myprobe[nix]
      mygid <- mygid[nix]
      mycoef <- mycoef[nix]
      gid1 <- rr$geneid2
      if(is.null(names(gid1))) { stop("problem with annotations!") }
      gid2 <- rr$geneid1
      if(is.null(names(gid2))) { stop("problem with annotations!") }
      data <- rr$data1

      #change the names of probes in x and data
      names(mycoef) <- names(mygid) <- mygid <- names(myprobe) <- myprobe <- as.character(gid1)
      dimnames(data)[[2]] <- as.character(gid2)
    } else { ## no mapping
      nix <- is.element(myprobe, dimnames(data)[[2]])
      myprobe <- myprobe[nix]
      mygid <- mygid[nix]
      mycoef <- mycoef[nix]
      gid1 <- gid2 <- mygid
      data <- data[ ,myprobe,drop=FALSE]
    }
    if(length(myprobe) == 0) {
      if(verbose) { message(sprintf("probe candidates: 0/%i", size)) }
      tt <- rep(NA, nrow(data))
      names(tt) <- dimnames(data)[[1]]
      return(list("score"=tt, "mapping"=c("mapped"=0, "total"=nrow(x)), "probe"=cbind("probe"=names(gid1), "EntrezGene.ID"=gid1, "new.probe"=names(gid2))))
    }

    if(size == 0 || size > nrow(x)) { size <- length(myprobe) }
    nix <- 1:size
    myprobe <- myprobe[nix]
    mygid <- mygid[nix]
    mycoef <- mycoef[nix]
    gid1 <- gid1[nix]
    gid2 <- gid2[nix]
    if(!is.na(cutoff)) {
      nix <- abs(mycoef) > cutoff
      myprobe <- myprobe[nix]
      mygid <- mygid[nix]
      mycoef <- mycoef[nix]
      gid1 <- gid1[nix]
      gid2 <- gid2[nix]
    }
    probe.candp <- myprobe[mycoef >= 0]
    probe.candn <- myprobe[mycoef < 0]
    gf <- length(myprobe)

    gf <- c("mapped"=gf, "total"=nrow(x))
    if(verbose) { message(sprintf("probe candidates: %i/%i",gf[1], gf[2])) }

    nprobe <- c(probe.candp, probe.candn)
    myw <- c("p"=length(probe.candp) / length(nprobe), "n"=length(probe.candn) / length(nprobe))
    res <- rep(0, nrow(data))

    if(signed) {
      ## consider only the sign of the coefficients
      if(length(probe.candp) > 0) { res <- myw["p"] * (apply(X=data[ ,probe.candp,drop=FALSE], MARGIN=1, FUN=sum, na.rm=TRUE) / apply(X=data[ ,probe.candp,drop=FALSE], MARGIN=1, FUN=function(x) { return(sum(!is.na(x))) })) }
      if(length(probe.candn) > 0) { res <- res - myw["n"] * (apply(X=data[ ,probe.candn,drop=FALSE], MARGIN=1, FUN=sum, na.rm=TRUE) / apply(X=data[ ,probe.candn,drop=FALSE], MARGIN=1, FUN=function(x) { return(sum(!is.na(x))) })) }
    } else {
      ## consider the exact value of the coefficients
      if(length(probe.candp) > 0) { res <- myw["p"] * (apply(X=data[ ,probe.candp,drop=FALSE], MARGIN=1, FUN=function(x, y) { nix <- is.na(x); return(sum(x * y, na.rm=TRUE) / sum(y[!nix])) }, y=abs(mycoef[probe.candp]))) }
      if(length(probe.candn) > 0) { res <- res - myw["n"] * (apply(X=data[ ,probe.candn,drop=FALSE], MARGIN=1, FUN=function(x, y) { nix <- is.na(x); return(sum(x * y, na.rm=TRUE) / sum(y[!nix])) }, y=abs(mycoef[probe.candn]))) }
    }
    return(list("score"=res, "mapping"=gf, "probe"=cbind("probe"=names(gid1), "EntrezGene.ID"=gid1, "new.probe"=names(gid2))))
}


