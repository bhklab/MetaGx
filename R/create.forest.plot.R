
create.forest.plot <- function(
                          survival.data, # a list of data frames
                          surv.time.colname,
                          surv.event.colname,
                          risk.val.colname,
                          stat=c("concordance.index","d.index"),
                          dataset.names = names(survival.data),
                          pooling.method = c("both", "fixed", "random"),
                          random.pooled.name = "Pooled estimate (random)",
                          fixed.pooled.name = "Pooled estimate (fixed)",
                          just.meta=FALSE,
                          x.ticks=NULL
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
  stat.vals <- sapply(stat.objects, function(x) x$d.index)
  }
  
  myspace <- "    "
  
  stat.se <- sapply(stat.objects, function(x) x$se)
  stat.lower <- sapply(stat.objects, function(x) x$lower)
  stat.upper <- sapply(stat.objects, function(x) x$upper)
  
  if(pooling.method=="random" || pooling.method=="both") {
    pooled.stat.random <- combine.est(stat.vals, stat.se, hetero=TRUE)
    dataset.names <- c(dataset.names, random.pooled.name)
    stat.vals <- c(stat.vals, meta.random=pooled.stat.random$estimate)
    stat.se <- c(stat.se, meta.random=pooled.stat.random$se)
    stat.lower <- c(stat.lower, meta.random=pooled.stat.random$estimate - pooled.stat.random$se*qnorm(0.975))
    stat.upper <- c(stat.upper, meta.random=pooled.stat.random$estimate + pooled.stat.random$se*qnorm(0.975))
  }
  if(pooling.method=="fixed" || pooling.method=="both") {
    pooled.stat.fixed <- combine.est(stat.vals, stat.se, hetero=FALSE)
    dataset.names <- c(dataset.names, fixed.pooled.name)
    stat.vals <- c(stat.vals, meta.fixed=pooled.stat.fixed$estimate)
    stat.se <- c(stat.se, meta.fixed=pooled.stat.fixed$se)
    stat.lower <- c(stat.lower, meta.fixed=pooled.stat.fixed$estimate - pooled.stat.fixed$se*qnorm(0.975))
    stat.upper <- c(stat.upper, meta.fixed=pooled.stat.fixed$estimate + pooled.stat.fixed$se*qnorm(0.975))
  }
  
  labeltext <- cbind(dataset.names, c(rep(myspace,length(dataset.names))))
  
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
    zero <- 1
    xlab <- "D-index"
  }
  if(!just.meta) {
    survcomp::forestplot.surv(labeltext=labeltext, mean=stat.vals, lower=stat.lower, upper=stat.upper, is.summary = c(rep(FALSE, length(stat.vals)-num.summary.stats), rep(TRUE, num.summary.stats)), zero=zero, xlab=xlab)
  } else {
    survcomp::forestplot.surv(labeltext=labeltext[c(nrow(labeltext)-1, nrow(labeltext)),], mean=stat.vals[c(nrow(labeltext)-1, nrow(labeltext))], lower=stat.lower[c(nrow(labeltext)-1, nrow(labeltext))], upper=stat.upper[c(nrow(labeltext)-1, nrow(labeltext))], is.summary = c(FALSE, FALSE), zero=zero, xlab=xlab, x.ticks=x.ticks)
  }
  return(list(stat.vals=stat.vals, stat.se=stat.se, stat.lower=stat.lower, stat.upper=stat.upper))
}