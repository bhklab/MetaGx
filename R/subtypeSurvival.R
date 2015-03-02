########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

########################
## Natchar Ratanasirigulchai
## changes made on March 2, 2015
########################


`subtypeSurvival` <- 
function (eset, sig, plot=FALSE, weighted=FALSE, surv.type=c("dfs", "rfs", "dmfs", "tdm", "os"), time.cens, condensed=TRUE, resdir="cache", nthread=1, sig.method, sig.scaling) {

  ######################
  
  concIndex <- function (x, stime, sevent, strat, weights, tau, alpha=0.05, alternative=c("two.sided", "less", "greater")) {
    if (missing(strat)) { strat <- array(1, dim=length(stime), dimnames=list(names(stime))) }
    if (missing(weights)) { weights <- array(1, dim=length(stime), dimnames=list(names(stime))) }
    resnn <- c("Dxy", "cindex", "se", "lower", "upper", "p.value")
    ccix <- complete.cases(x, stime, sevent, strat, weights)
    if (sum(ccix) < 3) {
      rr <- as.list(array(NA, dim=length(resnn), dimnames=list(resnn)))
    } else {      
      if (missing(tau)) {
        tau <- max(stime, na.rm=TRUE)
      } else {
        ss <- survcomp::censor.time(surv.time=stime, surv.event=sevent, time.cens=tau)  
        stime <- ss[[1]]
        sevent <- ss[[2]]
      }
      if (missing(strat)) { strat <- array(1, dim=length(stime), dimnames=list(names(stime))) }
      if (length(stime) != length(sevent) || length(stime) != length(strat) || length(stime) != length(x)) { stop("stime, sevent, strat and x must have the same length") }
      # rr <-  mRMRe::correlate(X=x, Y=Surv(stime, sevent), method="cindex", strata=strat, weights=weights)
      dd <- data.frame("stime"=stime, "sevent"=sevent, "x"=x, "weights"=weights, "strat"=strat, stringsAsFactors=FALSE)
      ## weights should be > 0
      dd <- dd[!is.na(dd$weights) & dd$weights > 0, , drop=FALSE]
      rr <- summary(survival::coxph(Surv(stime, sevent) ~ strata(strat) + x, data=dd, weights=dd$weights))
      cindex <- abs(rr$concordance[1] - 0.5)
      if (rr$coefficients["x", 1] < 0) { cindex <- - cindex }
      cindex <- 0.5 + cindex
      se <- rr$concordance[2]
      ci <- qnorm(p=alpha / 2, lower.tail=FALSE) * se
      lower <- cindex - ci
      upper <- cindex + ci
      switch(alternative, 
        "two.sided"={ p <- pnorm((cindex - 0.5) / se, lower.tail=cindex < 0.5) * 2 }, 
        "less"={ p <- pnorm((cindex - 0.5) / se, lower.tail=TRUE) }, 
        "greater"={  p <- pnorm((cindex - 0.5) / se, lower.tail=FALSE) }
      )
      rr <- c(list(2 * (cindex - 0.5)), cindex, se, lower, upper, p)
      names(rr) <- resnn
    }
    return (rr)
  }
  
  dIndex <- function (x, stime, sevent, strat, weights, tau, alternative=c("two.sided", "less", "greater")) {
    if (missing(strat)) { strat <- array(1, dim=length(stime), dimnames=list(names(stime))) }
    if (missing(weights)) { weights <- array(1, dim=length(stime), dimnames=list(names(stime))) }
    resnn <- c("Dindex", "coef", "se", "lower", "upper", "p.value")
    ccix <- complete.cases(x, stime, sevent, strat, weights)
    if (sum(ccix) < 3) {
      rr <- as.list(array(NA, dim=length(resnn), dimnames=list(resnn)))
    } else {      
      if (missing(tau)) {
        tau <- max(stime, na.rm=TRUE)
      } else {
        ss <- survcomp::censor.time(surv.time=stime, surv.event=sevent, time.cens=tau)  
        stime <- ss[[1]]
        sevent <- ss[[2]]
      }
      if (length(stime) != length(sevent) || length(stime) != length(strat) || length(stime) != length(x)) { stop("stime, sevent, strat and x must have the same length") }
      rr <-  survcomp::D.index(x=x, surv.time=stime, surv.event=sevent, strat=strat, weights=weights, method.test="logrank", na.rm=TRUE)
      rr <- rr[c("d.index", "coef", "se", "lower", "upper", "p.value")]
      names(rr) <- resnn
    }
    return (rr)
  }
  
  ######################
 
  surv.type <- match.arg(surv.type)
  
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
  
  if (missing(time.cens)) {
    time.cens <- max(stime, na.rm=TRUE)
  }
  
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
  message("build matrix sig scores")
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
  rownames(expr) <- names(sig)
  ## extract survival data
#   stime <- Biobase::pData(eset)[ , sprintf("t.%s", surv.type)] / 365
  switch(surv.type,
         "os" = {
           stime <- pData(eset)[, "days_to_death"]/365
         },
         "rfs" = {
           stime <- pData(eset) [,"days_to_tumor_recurrence"]/365
         })
  time.cens <- time.cens / 365
#   sevent <- Biobase::pData(eset)[ , sprintf("e.%s", surv.type)]
  switch(surv.type,
       "os" = {
         sevent <-as.numeric(pData(eset)[, "vital_status"] == "living")
       },
       "rfs" = {
         sevent <- pData(eset) [,"recurrence_status"]
       })
  ss <- survcomp::censor.time(surv.time=stime, surv.event=sevent, time.cens=time.cens)  
  stime <- ss[[1]]
  sevent <- ss[[2]]
  ## strata
  strat <- as.factor(Biobase::pData(eset)[ , "dataset"])
  message("cindex")
  ## concordance index
  splitix <- parallel::splitIndices(nx=nrow(expr), ncl=nthread)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(x, expr, stime, sevent, strat, sbts.proba) {
    ci <- lapply(rownames(expr)[x], function (x, expr, stime, sevent, strat, sbts.proba) {
      res <- t(apply(sbts.proba, 2, function (w, xx, stime, sevent, strat) {
        return (unlist(concIndex(x=xx, stime=stime, sevent=sevent, strat=strat, weights=w, alternative="two.sided")))
      }, xx=expr[x, ], stime=stime, sevent=sevent, strat=strat))
      return (res)
    }, expr=expr, stime=stime, sevent=sevent, strat=strat, sbts.proba=sbts.proba)
  }, expr=expr, stime=stime, sevent=sevent, strat=strat, sbts.proba=sbts.proba)
  rr <- unlist(mcres, recursive=FALSE)
  names(rr) <- names(sig)
  ## save results
  dd <- lapply(rr, data.frame)
  if (condensed) {
    WriteXLS::WriteXLS(x="dd", ExcelFileName=file.path(resdir, sprintf("subtype_gene_cindex.xls")), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
    ## per subtype
    dd2 <- lapply(sbtu, function (x, y) {
      res <- t(sapply(y, function (y, x) {
        return (y[x, ])
      }, x=x))
      return (data.frame(res))
    }, y=rr)
    names(dd2) <- sbtu
    WriteXLS::WriteXLS(x="dd2", ExcelFileName=file.path(resdir, sprintf("subtype_cindex.xls")), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
  } else {
    mapply(function(x, y, method, resdir) {
       WriteXLS::WriteXLS("x", ExcelFileName=file.path(resdir, sprintf("subtype_cindex_%s.xls", y)), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
     }, x=dd, y=names(dd), method=method, resdir=resdir)
  }
  cindices <- rr
  
  ## D index (hazard ratio)
message("dindex")
  splitix <- parallel::splitIndices(nx=nrow(expr), ncl=nthread)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(x, expr, stime, sevent, strat, sbts.proba) {
    ci <- lapply(rownames(expr)[x], function (x, expr, stime, sevent, strat, sbts.proba) {
      res <- t(apply(sbts.proba, 2, function (w, xx, stime, sevent, strat) {
         return (unlist(dIndex(x=xx, stime=stime, sevent=sevent, strat=strat, weights=w, alternative="two.sided")))
      }, xx=expr[x, ], stime=stime, sevent=sevent, strat=strat))
      return (res)
    }, expr=expr, stime=stime, sevent=sevent, strat=strat, sbts.proba=sbts.proba)
  }, expr=expr, stime=stime, sevent=sevent, strat=strat, sbts.proba=sbts.proba)
  rr <- unlist(mcres, recursive=FALSE)
  names(rr) <- names(sig)
  ## save results
  dd <- lapply(rr, data.frame)
  if (condensed) {
    WriteXLS::WriteXLS(x="dd", ExcelFileName=file.path(resdir, sprintf("subtype_gene_dindex.xls")), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
    ## per subtype
    dd2 <- lapply(sbtu, function (x, y) {
      res <- t(sapply(y, function (y, x) {
        return (y[x, ])
      }, x=x))
      return (data.frame(res))
    }, y=rr)
    names(dd2) <- sbtu
     WriteXLS::WriteXLS(x="dd2", ExcelFileName=file.path(resdir, sprintf("subtype_dindex.xls")), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
  } else {
    mapply(function(x, y, method, resdir) {
       WriteXLS::WriteXLS("x", ExcelFileName=file.path(resdir, sprintf("subtype_dindex_%s.xls", y)), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
     }, x=dd, y=names(dd), method=method, resdir=resdir)
  }
  dindices <- rr
  
  ## kaplan-meier survival curves
  if (plot) {
    figsize <- 6
    nc <- 3
    nr <- ceiling(ncol(sbts.proba) / nc)
    if (condensed) { pdf(file.path(resdir, "subtype_surv_curves.pdf"), height=nr * figsize, width=nc * figsize) }
    lapply(rownames(expr), function (x, expr, stime, sevent, strat, condensed, sbts.proba, sbts.crisp, subtype.col, resdir, nc, nr) {
      if (!condensed) { pdf(file.path(resdir, sprintf("subtype_surv_curves_%s.pdf", x$symbol)), height=nr * figsize, width=nc * figsize) }
      par(mfrow=c(nr, nc))
      for (i in 1:ncol(sbts.proba)) {
        w <- sbts.proba[ , i]
        ## discretize scores
        cc <- quantile(expr[x, !is.na(sbts.crisp[ , i]) & sbts.crisp[ , i] == 1], probs=c(0, 0.33, 0.66, 1), na.rm=TRUE)
        xx <- factor(cut(x=expr[x, ], breaks=cc, labels=FALSE))
        xx2 <- (as.numeric(xx) - 1) / length(levels(xx))
        dd <- data.frame("stime"=stime, "sevent"=sevent, "risk"=xx, "score"=xx2, "weights"=w, "strat"=strat, stringsAsFactors=FALSE)
        ## weights should be > 0
        dd <- dd[!is.na(dd$weights) & dd$weights > 0, , drop=FALSE]
        statn.risk <- summary(survival::coxph(formula=Surv(stime, sevent) ~ risk + strata(strat) , data=dd, weights=dd$weights))
        statn.score <- summary(survival::coxph(formula=Surv(stime, sevent) ~ score + strata(strat) , data=dd, weights=dd$weights))
        statn <- sprintf("Logrank P = %.1E\nHR = %.2g [%.2g,%2.g], P = %.1E", statn.risk$sctest["pvalue"], statn.score$conf.int[1], statn.score$conf.int[3], statn.score$conf.int[4], statn.score$coefficients[ , 5])
        survcomp::km.coxph.plot(formula.s=Surv(stime, sevent) ~ risk, data.s=dd, weight.s=dd$weights, x.label="time (years)", y.label="probability of disease-free survival", main.title=sprintf("%s\n%s", rownames(expr)[x], colnames(sbts.proba)[i]), leg.text=paste(c("Low", "Intermediate", "High"), "     ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen", "darkred"), .lty=c(1,1,1), show.n.risk=TRUE, n.risk.step=2, n.risk.cex=0.85, bty="n", leg.bty="n", o.text=statn, verbose=FALSE)
      }     
      if (!condensed) { dev.off() }
    }, expr=expr, stime=stime, sevent=sevent, strat=strat, condensed=condensed, sbts.proba=sbts.proba, sbts.crisp, subtype.col=subtype.col, resdir=resdir, nc=nc, nr=nr)
    if (condensed) { dev.off() }
  }

}


