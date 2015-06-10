
create.survival.plot <- function(
                                 surv.obj,
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
                                 # Most relevant legend parameters are already covered, but allow further parameters
                                 legend.par=list(),
                                 ...) {
  
  if(is.null(strata.names)) {
    strata.names <- names(surv.obj$strata)
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
    ctla4.surv.obj,
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
}