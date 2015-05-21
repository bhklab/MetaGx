getVerhaakSubtypes <- function(eset) {
  ## Load gene sets from the original publication
  # Load Verhaak et al. supplementary from the package inst directory
	#supplementary.data <- read.xls(system.file(file.path("extdata", "JCI65833sd1.xls"), package="MetaGx"), sheet=7, skip=1)
  # Use this instead when running this method from source
	supplementary.data <- read.xls("../inst/extdata/JCI65833sd1.xls", sheet=7, skip=1)

	genesets <- lapply(levels(supplementary.data$CLASS), function(y) as.character(supplementary.data[supplementary.data$CLASS==y,1]))
	names(genesets) <-  levels(supplementary.data$CLASS)
  
	expression.matrix <- exprs(eset)
  rownames(expression.matrix) <- fData(eset)$gene[match(rownames(expression.matrix), rownames(fData(eset)))]
  
  ## Get ssGSEA subtype scores
	gsva.out <- gsva(expression.matrix, genesets, method="ssgsea", min.sz=10, tau=0.75, parallel.sz=1)
  gsva.out <- t(gsva.out)
  ## Classify each sample according to the max ssGSEA subtype score. Note that this differs slightly
  # from the Veerhak et al. classification which has a "first pass" for classifying
  # Immunoreactive and Mesenchymal classes, followed by classification of remaining samples into
  # the max-scoring ssGSEA value for the four subtypes.
  subclasses <- as.factor(apply(gsva.out, 1, function(x) colnames(gsva.out)[which.max(x)]))
  ## Append a new column for Verhaak subtypes
  pData(eset) <- data.frame(pData(eset), Verhaak.subtypes=subclasses)
  return(list(Annotated.eset=eset, gsva.out=gsva.out))
}