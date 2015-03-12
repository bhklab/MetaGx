########################/
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

########################
## Natchar Ratanasirigulchai
## Changes made on: March 12, 2015
########################

`subtypeAssociation` <- 
function (eset, sig, plot=TRUE, subtype.col, weighted=FALSE, condensed=TRUE, resdir="cache", nthread=1, label=c("symbol", "entrez"), ...) {
  ## assess association between gene expression and subtypes
  #
  # Arga:
  #   eset: an expressionSet object
  #   sig: signature (vector of gene identifiers and optionnally coefficients; see sigScore) gene or list of signatures for signatures. If missing, all individual genes will be considered.
  #   subtype.col: color for each molecular subtype
  #
  # Returns
  #   list containing p-values for comparisons
  #   Kruskal-Wallist to test whether the expression of the genes(s) of interest is dependent on the molecular subtypes
  #   pairwise wilcoxon rank sum test p-values, is the expression of the gene(s) of interest higher in the subtype in rows compared to the subtype in column?
  
  label <- match.arg(label)
  
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
  
  if (!file.exists(file.path(resdir))) { dir.create(file.path(resdir), showWarnings=FALSE, recursive=TRUE) }
  
  ## extract subtypes
  sbts <- getSubtype(eset=eset, method="class")
  sbtu <- levels(sbts)
  if (sum(table(sbts) > 3) < 2) {
    warning("Not enough tumors in each subtype")
    return(NULL)
  }
  
  ## extract subtypes
  sbts <- getSubtype(eset=eset, method="class")
  if (sum(table(sbts) > 3) < 2) {
    warning("Not enough tumors in each subtype")
    return(NULL)
  }
  if (!weighted) {
    weights <- array(1, dim=length(sbts), dimnames=names(sbts))
  } else {
    weights <- getSubtype(eset=eset, method="fuzzy")
    stop("Weighted association analysis is not implemented yet")
    ## use the kruskal_tes and wilcox_test function in library coin
  }
  if (missing(subtype.col)) {
    subtype.col <- rainbow(length(sbtu), alpha=0.6)
  } else {
    if (length(subtype.col) < length(sbtu)) {
      stop (sprintf("Not enough color for %i subtypes", length(sbtu)))
    }  
  }

  if (length(sig) > 100 && condensed) {
    warning("Condensed output files are not suitable for large queries (>100 genes)")
  }
  
  splitix <- parallel::splitIndices(nx=length(sig), ncl=nthread)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(x, sig, eset, sbts, sbtu, ...) {    
    pp <- lapply(sig[x], function (x, eset, sbts, sbtu, ...) {
      xx <- sigScore(eset=eset, sig=x, ...)
      ccix <- complete.cases(xx, sbts)
      if (length(unique(sbts[ccix])) < 2) {
        res <- NULL
      } else {
        ## kruskal-wallis test
        kt <- kruskal.test(x=xx, g=sbts)$p.value
        ## pairwise wilcoxon test
        wt <- matrix(NA, nrow=length(sbtu), ncol=length(sbtu), dimnames=list(sbtu, sbtu))
        wt1 <- pairwise.wilcox.test(x=xx, g=sbts, p.adjust.method="none", paired=FALSE, alternative="greater")$p.value
        wt2 <- pairwise.wilcox.test(x=xx, g=sbts, p.adjust.method="none", paired=FALSE, alternative="less")$p.value
        nix <- !is.na(wt1)
        wt[rownames(wt1), colnames(wt1)][nix] <- wt1[nix]
        nix <- !is.na(t(wt2))
        wt[colnames(wt2), rownames(wt2)][nix] <- t(wt2)[nix]
        diag(wt) <- 1
        res <- list("kruskal.pvalue"=kt, "wilcoxon.pvalue"=wt, "x"=xx)
      }
      return (res)
    }, eset=eset, sbts=sbts, sbtu=sbtu, ...)
    return (pp)
  }, sig=sig, eset=eset, sbts=sbts, sbtu=sbtu, ...)
  pp <- do.call(c, mcres)
  pp <- pp[!sapply(pp, is.null)]
  
  if(label == "symbol"){
    names(pp) <- fData(eset)[names(pp), "SYMBOL"]
  }

    for (i in 1:length(pp)) {
      pp[[i]]$name <- names(pp)[i]
    }

  # names(pp) <- names(sig)
  
  if (plot) {
    ## determine y-axis limits for all the boxplot
    mylim <- range(lapply(pp, function (x, sbts) {
      sbtsu <- sort(unique(sbts))
      res <- range(sapply(sbtsu, function (s, sbts, x) {
        return (range(boxplot.stats(x=x[!is.na(sbts) & sbts == s])$stats))
      }, sbts=sbts, x=x$x))
      return (res)
    }, sbts=sbts))
    if (condensed) { pdf(file.path(resdir, "subtype_association_boxplot.pdf")) }
    lapply(pp, function (x, condensed, sbts, subtype.col, resdir) {
      if (!condensed) { pdf(file.path(resdir, sprintf("subtype_association_boxplot_%s.pdf", x$symbol))) }
      par(las=2, mar=c(5, 4, 4, 2) + 0.1, xaxt="n")
      # dd <- c(list(" "=NA), list("  "=NA), dd2)
      # mylim <- round(range(x$x, na.rm=TRUE))
      # mylim[abs(mylim) < 2] <- c(-2, 2)[abs(mylim) < 2]
      mp <- boxplot(x$x ~ sbts, las=3, outline=FALSE, ylim=mylim, main=sprintf("%s", x$name), col=subtype.col)
      axis(1, at=1:length(mp$names), tick=TRUE, labels=TRUE)
      text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black"))
      text(x=1:length(mp$names), y=par("usr")[3], pos=3, labels=sprintf("n=%i", table(sbts)[mp$names]), col=c("black"))
      # title(sub=sprintf("Kruskal-Wallis p-value = %.1E", kt))
      # legend("topleft", legend=c(paste(mp$names, sprintf("(%i)", table(sbts)[mp$names]), sep=" "), "", sprintf("K-W p = %.1E", kt)), bty="n", cex=0.8)
      legend("topleft", legend=sprintf("Kruskal-Wallis p-value = %.1E", x$kruskal.pvalue), bty="n")
      if (!condensed) { dev.off() }
    }, condensed=condensed, sbts=sbts, subtype.col=subtype.col, resdir=resdir)
    if (condensed) { dev.off() }
  }
  
  ## write spreadsheets
  ## kruskal-wallis p-values
  dd <- sapply(pp, function(x) { return(x$kruskal.pvalue)})
  med <- t(sapply(pp, function (x, sbts, sbtu) {
    res <- sapply(sbtu, function (s, x, sbts) {
      return (median(x[sbts == s], na.rm=TRUE))
    }, x=x$x, sbts=sbts)
    return (res)
  }, sbts=sbts, sbtu=sbtu))
  colnames(med) <- paste("median.expression", sbtu, sep=".")
  dd <- data.frame("Kruskal.Wallis.pvalue"=dd, "Kruskal.Wallis.fdr"=p.adjust(dd, method="fdr"), med, stringsAsFactors=FALSE)
  WriteXLS::WriteXLS(x="dd", SheetNames="Kruska-Wallis", ExcelFileName=file.path(resdir, "subtype_association_kruskal.xls"), AdjWidth=FALSE, BoldHeaderRow=TRUE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
  ## wilcoxon p-values
  dd <- lapply(pp, function (x) { return (data.frame(x$wilcoxon.pvalue)) })
  if (condensed) {
    WriteXLS::WriteXLS(x="dd", ExcelFileName=file.path(resdir, "subtype_association_wilcoxon.xls"), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
  } else {
    mapply(function(x, y, resdir) {
       WriteXLS::WriteXLS("x", ExcelFileName=file.path(resdir, sprintf("subtype_association_wilcoxon_%s.xls", y)), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
     }, x=dd, y=names(dd), resdir=resdir)
  }
  
  return (pp)
}


