# From a MetaGxOvarian eset, return a dataframe with survival

get.survival.data.for.plotting <- function(
                                  eset,
                                  gene.names=NULL,
                                  entrez.ids=NULL,
                                  num.quantiles=2,
                                  survival.type=c("overall.survival", "tumor.recurrence"),
                                  time.unit=c("days", "months", "years"),
                                  include.sample=rep(TRUE, nrow(pData(eset))),
                                  discard.cases.with.incomplete.survival=TRUE,
                                  additional.colnames.to.keep=NULL,
                                  gene.name.as.output.colname=FALSE,
                                  verbose=TRUE
                                  ) {
  days.per.month <- 30.4368
  days.per.year <- 365.242
  survival.type=match.arg(survival.type)
  time.unit=match.arg(time.unit)
  
  eset <- eset[,include.sample]
  # Find the fData index of the gene 
  if(is.null(gene.names) && is.null(entrez.ids)) {
    stop("One of gene.name and entrez.id should be provided")
  }
  if(!is.null(gene.names) && !is.null(entrez.ids)) {
    stop("Exactly one of gene.name and entrez.id should be provided")
  }
  indexes <- integer(0)
  if(!is.null(gene.names)) {
    for(gene.name in gene.names) {
      index <- which(fData(eset)$gene == gene.name)
      if(length(index) == 0) {
        stop(paste("Could not find gene name", gene.name))
      }
      if(length(index) > 1) {
        stop(paste("Found more than one gene called", gene.name))
      }
      if(verbose) {
        message(paste0("For gene ", gene.name, ", using Entrez ID: ", fData(eset)$EntrezGene.ID[index]))
      }
      indexes <- c(indexes, index)
    }
  }
  if(!is.null(entrez.ids)) {
    for(entrez.id in entrez.ids) {
      index <- which(fData(eset)$EntrezGene.ID == entrez.id)
      if(length(index) == 0) {
        stop(paste("Could not find Entrez ID", entrez.id))
      }
      if(length(index) > 1) {
        stop(paste("Found more than one Entrez ID", entrez.id))
      }
      indexes <- c(indexes, index)
    }
  }
  
  # Get the expression quantile 
  expression.values <- t(exprs(eset)[indexes,,drop=FALSE])

  # Rename columns by gene name. If gene.name.as.output.colname is FALSE, then the colnames will be geneid.[entrez id]
  if(gene.name.as.output.colname == TRUE) {
      colnames(expression.values) <- fData(eset)$gene[indexes]
  }
  
  expression.quantiles <- as.data.frame(lapply(as.data.frame(expression.values), function(x) cut(x, breaks=quantile(x, probs=seq(0,1,length=num.quantiles+1)), include.lowest=TRUE)))
  if(num.quantiles==2) {
    expression.quantiles <- as.data.frame(lapply(expression.quantiles, function(x) {levels(x) <- c("Low", "High"); return(x)}))
  } else if(num.quantiles==3) {
    expression.quantiles <- as.data.frame(lapply(expression.quantiles, function(x) {levels(x) <- c("Low", "Mid", "High"); return(x)}))
  } else {
    expression.quantiles <- as.data.frame(lapply(expression.quantiles, function(x) {levels(x) <- paste0("Expression.quantile.", 1:num.quantiles); return(x)}))
  }
  
  colnames(expression.values) <- paste0(colnames(expression.values), ".expression")
  colnames(expression.quantiles) <- paste0(colnames(expression.quantiles), ".quantile")
  pData(eset) <- cbind(pData(eset), expression.values)
  pData(eset) <- cbind(pData(eset), expression.quantiles)
  
  colnames.to.keep <- colnames(pData(eset))[c(grep(".quantile", colnames(pData(eset)), fixed=TRUE), grep(".expression", colnames(pData(eset)), fixed=TRUE))]
  if("data.source" %in% colnames(pData(eset))) {
    colnames.to.keep <- c("data.source", colnames.to.keep)
  }
  if(!is.null(additional.colnames.to.keep)) {
    colnames.to.keep <- c(additional.colnames.to.keep, colnames.to.keep)
  }
  
  
  if(survival.type=="overall.survival") {
    colnames.to.keep <- c(colnames.to.keep, "days_to_death", "vital_status")
    survival.data <- pData(eset)[colnames.to.keep]
    survival.data$vital_status <- survival.data$vital_status == "deceased"
  } else if(survival.type=="tumor.recurrence") {
    colnames.to.keep <- c(colnames.to.keep, "days_to_tumor_recurrence", "recurrence_status")
    survival.data <- pData(eset)[colnames.to.keep]
    survival.data$recurrence_status <- survival.data$recurrence_status == "recurrence"
  }
  # Rename last two columns
  colnames(survival.data)[length(colnames(survival.data))-1] <- "days_to_event"
  colnames(survival.data)[length(colnames(survival.data))] <- "event_status"
  
  # Remove cases with missing survival data
  if(discard.cases.with.incomplete.survival) {
    if(verbose) {
      paste0(length(!is.na(survival.data$days_to_event) & !is.na(survival.data$event_status)), " cases out of ", length(survival.data$days_to_event), " have complete survival information and will be kept")
    }
    survival.data <- survival.data[!is.na(survival.data$days_to_event) & !is.na(survival.data$event_status),]
  }
  
  if(time.unit=="months") {
    survival.data$months_to_event <- survival.data$days_to_event / days.per.month
    survival.data$days_to_event <- NULL
  }
  if(time.unit=="years") {
    survival.data$years_to_event <- survival.data$days_to_event / days.per.year
    survival.data$days_to_event <- NULL
  }
  
  return(survival.data)
}