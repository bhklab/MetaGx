
create.survival.plot <- function(
                                 surv.time,
                                 surv.event,
                                 groups,
                                 risk.order.same.as.groups = TRUE, # for computing c-index and d-index, should the risk order be the same as the ordered factor 
                                 xlab="Time",
                                 ylab="Survival",
                                 main="Survival Plot",
                                 cex=0.8,
                                 col=RColorBrewer::brewer.pal(length(surv.obj$strata), name="Dark2"),
                                 strata.names=NULL,
                                 reverse.colour.order=FALSE,
                                 reverse.legend.order=FALSE,
                                 legend.pos="topright",
                                 legend.lwd=5,
                                 stats.to.show=c("n","p","d","c"), # an ordered vector of stats to show: n for the number of samples, p for p-value of the likelihood-ratio test, d for d-index, c for c-index
                                 # Most relevant legend parameters are already covered, but allow further parameters
                                 legend.par=list(),
                                 ...) {
  surv.obj <- survfit(Surv(surv.time, surv.event) ~ groups)
 
   if(is.null(strata.names)) {
    strata.names <- names(surv.obj$strata)
  }
  
  if(!is.null(stats.to.show) && !all(stats.to.show %in% c("n","p","d","c"))) {
    stop("Argument stats.to.show should only contain 'n', 'p', 'd', 'c'.")
  }
  if(any(c("c", "d") %in% stats.to.show)) {
    if(!is.ordered(groups)) {
      stop("For calculation of c-index and d-index, groups must be an ordered factor")
    }
    if(risk.order.same.as.groups == TRUE) {
      risk.prediction <- as.numeric(groups)
    } else {
      risk.prediction <- -as.numeric(groups)
    }
  }
  
  # In all cases, make sure that the associations between colours, legend labels, and survival curves are consistent
  plot.col=col
  legend.col=col
  legend.strata=strata.names
  if(reverse.colour.order == TRUE) {
    plot.col <- rev(plot.col)
    legend.col <- rev(legend.col)
  }
  if(reverse.legend.order == TRUE) {
    legend.col <- rev(legend.col)
    legend.strata <- rev(legend.strata)
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
    stop("Cannot specify legend strata separately from the plot - see parameter 'strata.names', 'reverse.legend.order'")
  }
  
  legend.par$col <- legend.col
  legend.par$legend <- legend.strata
  
  if(!("x" %in% names(legend.par))) {
    legend.par$x <- legend.pos
  }
  if(!("lwd" %in% names(legend.par))) {
    legend.par$lwd <- legend.lwd
  }
  
  do.call(legend, legend.par)
  
  # display stats
  coxph.summary <- summary(survival::coxph(Surv(surv.time, surv.event) ~ groups))
  n <- coxph.summary$n
  p <- coxph.summary$logtest[["pvalue"]]
  c <- survcomp::concordance.index(risk.prediction, surv.time, surv.event)$c.index
  d <- survcomp::D.index(risk.prediction, surv.time, surv.event)$d.index
  
  if(length(stats.to.show) > 0) {
    text.to.show <- ""
    for(i in 1:length(stats.to.show)) {
      if(i != 1) {
          text.to.show <- paste0(text.to.show, "\n")
      }
      if(stats.to.show[i] == "n") {
        text.to.show <- paste0(text.to.show, "n = ", n)
      } else if(stats.to.show[i] == "p") {
        text.to.show <- paste0(text.to.show, "p = ", round(p, digits=3))
      } else if(stats.to.show[i] == "c") {
        text.to.show <- paste0(text.to.show, "c-index: ", round(c, digits=3))
      } else if(stats.to.show[i] == "d") {
        text.to.show <- paste0(text.to.show, "d-index: ", round(d, digits=3))
      } 
    }
    invisible(plotrix::corner.label(text.to.show, x=-1, y=-1))
  }
}