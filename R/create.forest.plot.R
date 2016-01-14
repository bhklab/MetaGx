
create.forest.plot <- function(
                          survival.data, # a list of data frames
                          surv.time.colname,
                          surv.event.colname,
                          risk.val.colname,
                          stat=c("concordance.index","d.index", "hazard.ratio"), #returns log HR and D.index
                          dataset.names = names(survival.data),
                          pooling.method = c("both", "fixed", "random"),
                          random.pooled.name = "Pooled estimate (random)",
                          fixed.pooled.name = "Pooled estimate (fixed)",
                          x.ticks=NULL,
                          main="Forest Plot",
                          ...
                          ) {
  
  stat <- match.arg(stat)
  pooling.method <- match.arg(pooling.method)
  for(colname in c(surv.time.colname, surv.event.colname, risk.val.colname)) {
    if(!all(sapply(survival.data, function(x) colname %in% colnames(x)))) {
      stop(paste("Column name", colname, "is not present in all dataframes of survival.data"))
    }
  }
  
  if(stat == "concordance.index") {
  stat.objects <- lapply(survival.data, function(x) {
    survcomp::concordance.index(x=x[[risk.val.colname]], surv.time=x[[surv.time.colname]], surv.event=x[[surv.event.colname]], method='noether')
    })
  stat.vals <- sapply(stat.objects, function(x) x$c.index)
  }
  if(stat == "d.index") {
  stat.objects <- lapply(survival.data, function(x) {
    survcomp::D.index(x=x[[risk.val.colname]], surv.time=x[[surv.time.colname]], surv.event=x[[surv.event.colname]])
    })
  stat.vals <- sapply(stat.objects, function(x) x$dicoef)
  }
  if(stat == "hazard.ratio") {
  stat.objects <- lapply(survival.data, function(x) {
    survcomp::hazard.ratio(x=x[[risk.val.colname]], surv.time=x[[surv.time.colname]], surv.event=x[[surv.event.colname]])
    })
  stat.vals <- sapply(stat.objects, function(x) x$coef)
  }
  
  stat.se <- sapply(stat.objects, function(x) x$se)
  stat.lower <- sapply(stat.objects, function(x) x$lower)
  stat.upper <- sapply(stat.objects, function(x) x$upper)
  # If D-index or hazard ratio, log-transform upper and lower bounds. Note that
  # survcomp outputs stanard error and coef for log(D-index) and log(HR)
  if(stat == "d.index" || stat == "hazard.ratio") {
    stat.lower <- sapply(stat.lower, log)
    stat.upper <- sapply(stat.upper, log)
  }
  rma.random <- NULL
  rma.fixed <- NULL
  if(pooling.method=="random" || pooling.method=="both") {
    rma.random <- rma(stat.vals, sei=stat.se, method="DL", slab=names(survival.df.list))
  }
  if(pooling.method=="fixed" || pooling.method=="both") {
    rma.fixed <- rma(stat.vals, sei=stat.se, method="FE", slab=names(survival.df.list))
  }
  
  num.summary.stats <- 0
  if(pooling.method=="random" || pooling.method=="fixed") {
    num.summary.stats <- 1
  } else if(pooling.method=="both") {
    num.summary.stats <- 2
  }
  
  if(stat=="concordance.index") {
    zero <- 0.5
    xlab <- "Concordance Index"
  } else if(stat=="d.index") {
    zero <- 0
    xlab <- "D-index"
  } else if(stat=="hazard.ratio") {
    # Convert from natural log to log10
    zero <- 0
    xlab <- "Hazard Ratio"
  }
  
  if(pooling.method=="both") {
    forest(rma.fixed, mlab="Fixed Effects", xlab="Hazard Ratio", atransf=exp, refline = zero, ylim=c(-2.5,length(survival.df.list) + 3), main=main, annotate=FALSE, addfit=FALSE, digits=c(2,2), ...)
    abline(h=0, lwd=1)
    addpoly(rma.fixed, mlab="Fixed Effects", row=-1, atransf=exp, annotate=FALSE)
    addpoly(rma.random, mlab="Random Effects", row=-2, atransf=exp, annotate=FALSE)
  } else if(pooling.method=="fixed") {
    forest(rma.fixed, mlab="Fixed Effects", xlab="Hazard Ratio", atransf=exp, refline = zero, annotate=FALSE)
  } else if(pooling.method=="random") {
    forest(rma.random, mlab="Fixed Effects", xlab="Hazard Ratio", atransf=exp, refline = zero, annotate=FALSE)
  } 
  return(list(stat.vals=stat.vals, stat.se=stat.se, stat.lower=stat.lower, stat.upper=stat.upper, rma.random=rma.random, rma.fixed=rma.fixed))
}