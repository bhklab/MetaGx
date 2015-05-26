getTothillSubtypes <- function(eset, probe.mapping=c("probesets","entrez.ids")) {
  probe.mapping <- match.arg(probe.mapping)
  
	## Load train data with predefined class labels
  #supplementary.data <- read.table(system.file(file.path("extdata", "../inst/extdata/tothill.supptable.2.txt"), package="MetaGx"), header=TRUE)
  supplementary.data <- read.table("../inst/extdata/tothill.supptable.2.txt", header=TRUE)
  ## Train a diagonal linear discriminant classifier using the Tothill data set and overlapping probesets / entrez gene ids
  tothill.eset <- esets$GSE9891
  
  if(probe.mapping == "probesets") {
    intersecting.probesets <- intersect(rownames(exprs(tothill.eset)), supplementary.data$probeset)
  }
  else {
    ## TODO
  }
	## Train and classify with diagonal LDA
}