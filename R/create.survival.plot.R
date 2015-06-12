
create.survival.plot <- function(
                                 surv.time,
                                 surv.event,
                                 groups,
                                 datasets=NULL,
                                 risk.vals=NULL, # for computing c-index and d-index, should the risk order be the same as the ordered factor 
                                 xlab="Time",
                                 ylab="Survival",
                                 main="Survival Plot",
                                 cex=0.8,
                                 col=RColorBrewer::brewer.pal(length(surv.obj$strata), name="Dark2"),
                                 group.names=NULL,
                                 reverse.colour.order=FALSE,
                                 reverse.legend.order=FALSE,
                                 legend.pos="topright",
                                 legend.lwd=5,
                                 stats.to.show=c("n","p","d","c"), # an ordered vector of stats to show: n for the number of samples, p for p-value of the likelihood-ratio test, d for d-index, c for c-index
                                 # Most relevant legend parameters are already covered, but allow further parameters
                                 legend.par=list(),
                                 ...) {
  if(is.null(datasets)) {
    surv.obj <- survfit(Surv(surv.time, surv.event) ~ groups)
  } else {
    surv.obj <- survfit(Surv(surv.time, surv.event) ~ groups + strata(datasets))
  }
 
   if(is.null(group.names)) {
    group.names <- names(surv.obj$strata)
  }
  
  if(!is.null(stats.to.show) && !all(stats.to.show %in% c("n","p","d","c"))) {
    stop("Argument stats.to.show should only contain 'n', 'p', 'd', 'c'.")
  }
  if(any(c("c", "d") %in% stats.to.show)) {
    if(is.null(risk.vals)) {
      stop("For calculation of c-index and d-index, risk.vals must be provided")
    }
    if(!is.ordered(groups)) {
      stop("For calculation of c-index and d-index, groups must be an ordered factor")
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
  c <- survcomp::concordance.index(risk.vals, surv.time, surv.event)$c.index
  d <- survcomp::D.index(risk.vals, surv.time, surv.event)$d.index
  
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
      } else if(stats.to.show[i] == "c") {
        text.to.show <- paste0(text.to.show, "c-index: ", round(c, digits=3))
      } else if(stats.to.show[i] == "d") {
        text.to.show <- paste0(text.to.show, "d-index: ", round(d, digits=3))
      } 
    }
    invisible(plotrix::corner.label(text.to.show, x=-1, y=-1))
  }
}