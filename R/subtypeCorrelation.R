########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

########################
## Natchar Ratanasirigulchai
## Changes Made on: March 12, 2015
########################


`subtypeCorrelation` <- 
function (eset, sig, method=c("pearson", "spearman", "kendall"), weighted=FALSE, condensed=TRUE, plot=TRUE, resdir="cache", nthread=1, sig.method, sig.scaling, label=c("symbol", "entrez")) {
  ## assess (weighted) correlation between gene expression with respect to subtypes
  #
  # Arga:
  #   eset: an expressionSet object
  #   geneid: vector of Entrez Gene IDs. If missing, all genes will be considered.
  #   plot: should the correlation heatmap be plotted?
  #   method: method for correlation
  #   label: labels by symbol or entrezgene ID
  #   resdir
  #   nthread:
  #   ...: parameters to be passed to sigScore
  #
  # Returns
  #   list containing p-values for comparisons
  #   Kruskal-Wallist to test whether the expression of the genes(s) of interest is dependent on the molecular subtypes
  #   pairwise wilcoxon rank sum test p-values, is the expression of the gene(s) of interest higher in the subtype in rows compared to the subtype in column?
  
  if (class(eset) != "ExpressionSet") {
    stop("Handling list of expressionSet objects is not implemented yet")
  }
  
  if (missing(sig)) {
    sig <- as.list(rownames(Biobase::fData(eset)))
    #     sig <- as.list(rownames(Biobase::fData(eset)[ , "ENTREZID"]))
    ## assign gene symbol as signature names
    gsymb <- Biobase::fData(eset)[ , "SYMBOL"]
    gsymb[is.na(gsymb)] <- paste("ENTREZID", Biobase::fData(eset)[is.na(gsymb), "ENTREZID"], sep=".")
    names(sig) <- gsymb
  }
  if (!is.list(sig) || (is.list(sig) && is.data.frame(sig))) {
    sig <- list("SIG"=sig)
  }
  if (is.null(names(sig))) {
    names(sig) <- paste("SIG", 1:length(sig), sep=".")
  }
  
  label <- match.arg(label)
  
  if (!file.exists(file.path(resdir))) { dir.create(file.path(resdir), showWarnings=FALSE, recursive=TRUE) }
  
  ## for a single expressionSet object
  
  ## extract subtypes
  sbts <- getSubtype(eset=eset, method="class")
  if (sum(table(sbts) > 3) < 2) {
    warning("Not enough tumors in each subtype")
    return(NULL)
  }
  if (!weighted) {
    sbts.proba <- getSubtype(eset=eset, method="crisp")
  } else {
    sbts.proba <- getSubtype(eset=eset, method="fuzzy")  
  }
  sbts.proba <- cbind("Global"=1, sbts.proba)
  sbts.crisp <- getSubtype(eset=eset, method="crisp")
  sbts.crisp <- cbind("Global"=1, sbts.crisp)
  sbtu <- colnames(sbts.proba)
  
  ## build matrix of signature scores in parallel
  splitix <- parallel::splitIndices(nx=length(sig), ncl=nthread)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(x, sig, eset, sigm, sigs) {
    res <- lapply(sig[x], function (x, eset, sigm, sigs) {
      res <- sigScore(eset=eset, sig=x, method=sigm, scaling=sigs)
      return (res)
    }, eset=eset, sigm=sigm, sigs=sigs)
    return (res)
  }, sig=sig, eset=eset, sigm=sig.method, sigs=sig.scaling)
  expr <- t(do.call(cbind, do.call(c, mcres)))
  
 if(label == "symbol"){
   rownames(expr) <- fData(eset)[names(sig),"SYMBOL"]
 } else if(label == "entrez"){
   rownames(expr) <- names(sig)
 } else{
  rownames(expr) <- names(sig)
 }  
    
  ## compute subtype-specific pairwise correlation across query genes
  
  ## slow code
  # if (method == "spearman") {
  #   expr <- t(apply(expr, 1, rank))
  # }
  # pairs <- t(combn(1:length(gid), 2, simplify=TRUE))
  # splitix <- parallel::splitIndices(nx=nrow(pairs), ncl=nthread)
  # splitix <- splitix[sapply(splitix, length) > 0]
  # mcres <- parallel::mclapply(splitix, function(x, ...) {    
  #   res <- apply(pairs[x, , drop=FALSE], 1, function (x, ...) {
  #     res <- apply(sbts.proba, 2, function (w, x, expr) {
  #       return (wcor(d=t(expr[x, , drop=FALSE]), w=w))
  #     }, x=x, expr=expr)
  #     return (res)
  #   })
  #   return (res)
  # }, expr=expr, sbts.proba=sbts.proba)
  # res <- t(do.call(cbind, mcres))
  # rr <- unlist(apply(res, 2, function (x, y, gid) {
  #   rr <- matrix(NA, nrow=length(gid), ncol=length(gid), dimnames=list(gid, gid))
  #   rr[y] <- x
  #   rr[y[ , 2:1]] <- x
  #   diag(rr) <- 1
  #   return (list(rr))
  # }, y=pairs, gid), recursive=FALSE)

  ## using mRMRe
  nn <- mRMRe::get.thread.count()
  mRMRe::set.thread.count(nthread)
  rr <- unlist(apply(sbts.proba, 2, function (w, expr, method) {
    expr2 <- mRMRe::mRMR.data(data=data.frame(t(expr)), weights=w)
    cor.genes <- mRMRe::mim(object=expr2, continuous_estimator=method, method="cor")
    return (list(cor.genes))
  }, expr=expr, method=method), recursive=FALSE)
  mRMRe::set.thread.count(nn)
  
  if (plot) {
    
    # Function to plot color bar
    color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
      scale = (length(lut)-1)/(max-min)
      
#       dev.new(width=1.75, height=5)
      plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='correlation', main=title)
      axis(2, ticks, las=1)
      for (i in 1:(length(lut)-1)) {
        y = (i-1)/scale + min
        rect(0,y,10,y+1/scale, col=lut[i], border=NA)
      }	
    }
    
    if (condensed) {
      nc <- 4
      nr <- ceiling(length(rr) / nc)
      pdf(file.path(resdir, "legend.pdf"), width= 1.75, height=5)
      color.bar(colorRampPalette(c("blue","white","red"))(256), -1)
      
      dev.off()
      pdf(file.path(resdir, "subtype_correlation.pdf"), width=nc * 5, height=nr * 5)
      par(mfrow=c(nr, nc), cex=0.8, mar=c(7, 7, 3, 1) + 0.1, xaxt="n", yaxt="n")
      for (i in 1:length(rr)) {
        xx <- rr[[i]]
        image(z=xx, zlim=c(-1, 1), col=colorRampPalette(c("blue","white","red"))(256))
        # axis(1, at=seq(1, nrow(xx), by=1), labels=FALSE)
        text(x=seq(par("usr")[1] + par("usr")[2] * 0.075, par("usr")[2] - par("usr")[2] * 0.025, length.out=nrow(xx)), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=rownames(xx), srt=45, xpd=NA, cex=0.8, font=2)
        text(x=par("usr")[1], y=seq(par("usr")[3] + par("usr")[4] * 0.05, par("usr")[4] - par("usr")[4] * 0.05, length.out=nrow(xx)), pos=2, labels=rownames(xx), srt=0, xpd=NA, cex=0.8, font=2)
        if (names(rr)[i] == "Global") {
          ss <- sprintf("Co-expression in %s population", names(rr)[i])
        } else {
          ss <- sprintf("Co-expression in %s subtype", names(rr)[i])
        }
        title(main=ss)
      }
      
      dev.off()
    } else {
      for (i in 1:length(rr)) {
        pdf(file.path(resdir, sprintf("subtype_correlation_%s.pdf", names(rr)[i])), width=5, height=5)
        par(cex=0.8, mar=c(7, 7, 3, 1) + 0.1, xaxt="n", yaxt="n")
        xx <- rr[[i]]
        image(z=xx, zlim=c(-1, 1), col=blueYellow(256))
        # axis(1, at=seq(1, nrow(xx), by=1), labels=FALSE)
        text(x=seq(par("usr")[1] + par("usr")[2] * 0.075, par("usr")[2] - par("usr")[2] * 0.025, length.out=nrow(xx)), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=rownames(xx), srt=45, xpd=NA, cex=0.8, font=2)
        text(x=par("usr")[1], y=seq(par("usr")[3] + par("usr")[4] * 0.05, par("usr")[4] - par("usr")[4] * 0.05, length.out=nrow(xx)), pos=2, labels=rownames(xx), srt=0, xpd=NA, cex=0.8, font=2)
        if (names(rr)[i] == "Global") {
          ss <- sprintf("Co-expression for %s population", names(rr)[i])
        } else {
          ss <- sprintf("Co-expression for %s subtype", names(rr)[i])
        }
        title(main=ss)
        color.bar(colorRampPalette(c("blue", "yellow"))(256), -1)
        dev.off()
      }
    }
  }
  
  dd <- lapply(rr, function (x, glabel) {
    dd <- data.frame(x)
    return (dd)
  }, glabel=colnames(expr))
  if (condensed) {
    WriteXLS::WriteXLS(x="dd", ExcelFileName=file.path(resdir, sprintf("subtype_correlation_%s.xls", method)), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
  } else {
    mapply(function(x, y, method, resdir) {
       WriteXLS::WriteXLS("x", ExcelFileName=file.path(resdir, sprintf("subtype_correlation_%s_%s.xls", method, y)), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
     }, x=dd, y=names(dd), method=method, resdir=resdir)
  }
  return (rr)
}


