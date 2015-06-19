
create.survival.plot <- function(
                                 surv.time,
                                 surv.event,
                                 groups,
                                 datasets=NULL,
                                 risk.vals=NULL, # for computing c-index and d-index, should the risk order be the same as the ordered factor 
                                 xlab="Time",
                                 ylab="Survival",
                                 main="Survival Plot",
                                 cex=0.5,
                                 time.cens=NULL,
                                 col=RColorBrewer::brewer.pal(length(surv.obj$strata), name="Dark2"),
                                 group.names=NULL, # the names to use on the legend. The order of names should correspond with levels(groups)
                                 reverse.colour.order=FALSE,
                                 reverse.legend.order=FALSE,
                                 legend.pos="topright",
                                 legend.lwd=5,
                                 stats.to.show=c("n","p","d","c", "hr"), # an ordered vector of stats to show: n for the number of samples, p for p-value of the likelihood-ratio test, d for d-index, c for c-index
                                 show.confidence.intervals=TRUE,
                                 pooling.method=c("random", "fixed"), #pooling method for concordance index, and D-index
                                 # Most relevant legend parameters are already covered, but allow further parameters
                                 legend.par=list(),
                                 ...) {
  pooling.method = match.arg(pooling.method)
  
  surv.time.to.plot <- surv.time
  surv.event.to.plot <- surv.event
  
  if(!is.null(time.cens)) {
    censored.out <- survcomp::censor.time(surv.time=surv.time, surv.event=surv.event, time.cens=time.cens)
    surv.time.to.plot <- censored.out$surv.time.cens
    surv.event.to.plot <- censored.out$surv.event.cens
  }
  
    surv.obj <- survfit(Surv(surv.time.to.plot, surv.event.to.plot) ~ groups)
  #} else {
  #  surv.obj <- survfit(Surv(surv.time, surv.event) ~ groups + strata(datasets))
  #}
 
   if(is.null(group.names)) {
    group.names <- names(surv.obj$strata)
  }
  
  #if(!is.null(stats.to.show) && !all(stats.to.show %in% c("n","p","d","c"))) {
  #  stop("Argument stats.to.show should only contain 'n', 'p', 'd', 'c'.")
  #}
  if(any(c("c", "d") %in% stats.to.show)) {
    if(is.null(risk.vals)) {
      stop("For calculation of c-index and d-index, risk.vals must be provided")
    }
  }
  
  # In all cases, make sure that the associations between colours, legend labels, and survival curves are consistent
  plot.col=col
  legend.col=col
  legend.groups=group.names
  if(reverse.colour.order == TRUE) {
    plot.col <- rev(plot.col)
    legend.col <- rev(legend.col)
  }
  if(reverse.legend.order == TRUE) {
    legend.col <- rev(legend.col)
    legend.groups <- rev(legend.groups)
  }
  
  plot(
    surv.obj,
    xlab=xlab,
    ylab=ylab,
    main=main,
    cex=cex,
    col=plot.col,
    ...)
  
  # Make sure legend.par does not contain col or legend parameters
  if("col" %in% names(legend.par)) {
    stop("Cannot specify legend colour separately from the plot - see parameters 'col', 'reverse.colour.order', 'reverse.legend.order'")
  }
  if("legend" %in% names(legend.par)) {
    stop("Cannot specify legend names separately from the plot - see parameter 'group.names', 'reverse.legend.order'")
  }
  
  legend.par$col <- legend.col
  legend.par$legend <- legend.groups
  legend.par$bty <- "n"
  
  if(!("x" %in% names(legend.par))) {
    legend.par$x <- legend.pos
  }
  if(!("lwd" %in% names(legend.par))) {
    legend.par$lwd <- legend.lwd
  }
  
  do.call(legend, legend.par)
  
  # display stats
  if(is.null(datasets)) {
    coxph.summary <- summary(survival::coxph(Surv(surv.time, surv.event) ~ groups))
  } else {
    coxph.summary <- summary(survival::coxph(Surv(surv.time, surv.event) ~ groups + strata(datasets)))
  }
  n <- coxph.summary$n
  p <- coxph.summary$logtest[["pvalue"]]
  
  hr.out <- survcomp::hazard.ratio(groups, surv.time, surv.event, strat=datasets)
  
  if(length(stats.to.show) > 0) {
    text.to.show <- ""
    for(i in 1:length(stats.to.show)) {
      if(i != 1) {
          text.to.show <- paste0(text.to.show, "\n")
      }
      if(stats.to.show[i] == "n") {
        text.to.show <- paste0(text.to.show, "n = ", n)
      } else if(stats.to.show[i] == "p") {
        text.to.show <- paste0(text.to.show, "Likelihood ratio test: p = ", round(p, digits=3))
      } else if(stats.to.show[i] == "hr") {
        if(length(hr.out$hazard.ratio) == 1) {
          text.to.show <- paste0(text.to.show, sprintf("HR: %.3f, 95%% CI: [%.3f-%.3f]", hr.out$hazard.ratio, hr.out$lower, hr.out$upper))
        } else {
          for(j in 1:length(hr.out$hazard.ratio)) {
            if(j != 1) {
                text.to.show <- paste0(text.to.show, "\n")
              }
            text.to.show <- paste0(text.to.show, sprintf("HR %d: %.3f, 95%% CI: [%.3f-%.3f]", j, hr.out$hazard.ratio[j], hr.out$lower[j], hr.out$upper[j]))
          }
        }
      } else if(stats.to.show[i] == "c") {
        if(is.null(datasets) || length(unique(datasets)) == 1) {
          ci.out <- survcomp::concordance.index(risk.vals, surv.time, surv.event, method='noether', strat=datasets)
          c.index <- ci.out$c.index
          c.lower <- ci.out$lower
          c.upper <- ci.out$upper
        } else {
          # pooled stat via meta-analysis
          stat.objects <- lapply(unique(datasets), function(current.dataset) {
              survcomp::concordance.index(x=risk.vals[datasets == current.dataset], surv.time=surv.time[datasets == current.dataset], surv.event=surv.event[datasets == current.dataset], method='noether')
              })
          names(stat.objects) <- unique(datasets)
          stat.vals <- sapply(stat.objects, function(x) x$c.index)
          stat.se <- sapply(stat.objects, function(x) x$se)
          if(pooling.method=="random") {
            pooled.stat <- combine.est(stat.vals, stat.se, hetero=TRUE)
          } else { # fixed effects
            pooled.stat <- combine.est(stat.vals, stat.se, hetero=FALSE)
          }
          c.index = pooled.stat$estimate
          c.lower <- pooled.stat$estimate - pooled.stat$se*qnorm(0.975)
          c.upper <- pooled.stat$estimate + pooled.stat$se*qnorm(0.975)
          stat.vals <- c(stat.vals, meta.random=pooled.stat$estimate)
          stat.se <- c(stat.se, meta.random=pooled.stat$se)
        }
        if(show.confidence.intervals) {
          text.to.show <- paste0(text.to.show, sprintf("Concordance index: %.3f, 95%% CI: [%.3f-%.3f]", c.index, c.lower, c.upper))
        } else {
          text.to.show <- paste0(text.to.show, sprintf("Concordance index: %.3f", c.index))
        }
      } else if(stats.to.show[i] == "d") {
        if(is.null(datasets) || length(unique(datasets)) == 1) {
          di.out <- survcomp::D.index(risk.vals, surv.time, surv.event, strat=datasets)
          d.index <- di.out$d.index
          d.lower <- di.out$lower
          d.upper <- di.out$upper
        } else {
          # pooled stat via meta-analysis
          stat.objects <- lapply(unique(datasets), function(current.dataset) {
              survcomp::D.index(x=risk.vals[datasets == current.dataset], surv.time=surv.time[datasets == current.dataset], surv.event=surv.event[datasets == current.dataset])
              })
          names(stat.objects) <- unique(datasets)
          stat.vals <- sapply(stat.objects, function(x) x$d.index)
          stat.se <- sapply(stat.objects, function(x) x$se)
          if(pooling.method=="random") {
            pooled.stat <- combine.est(stat.vals, stat.se, hetero=TRUE)
          } else { # fixed effects
            pooled.stat <- combine.est(stat.vals, stat.se, hetero=FALSE)
          }
          d.index = pooled.stat$estimate
          d.lower <- pooled.stat$estimate - pooled.stat$se*qnorm(0.975)
          d.upper <- pooled.stat$estimate + pooled.stat$se*qnorm(0.975)
          stat.vals <- c(stat.vals, meta.random=pooled.stat$estimate)
          stat.se <- c(stat.se, meta.random=pooled.stat$se)
        }
        if(show.confidence.intervals) {
          text.to.show <- paste0(text.to.show, sprintf("D-Index: %.3f, 95%% CI: [%.3f-%.3f]", d.index, d.lower, d.upper))
        } else {
          text.to.show <- paste0(text.to.show, sprintf("D-Index: %.3f", d.index))
        }
      } 
    }
    invisible(plotrix::corner.label(text.to.show, x=-1, y=-1))
  }
  return(stat.objects)
}