source("../R/datasetMerging.R")
source("../R/getSubtype.R")

.getPooledEset <- function(only.hgs=TRUE) {
  if(only.hgs) {
    source("../inst/extdata/hgs.patientselection.config")
  } else {
    source(system.file("extdata", "patientselection.config", package="MetaGxOvarian"))
  }
  source(system.file("extdata", "createEsetList.R", package="MetaGxOvarian"))
  esets <- lapply(esets, function(x) {
    factor.indices <- sapply(pData(x), is.factor)
    pData(x)[factor.indices] <- lapply(pData(x)[factor.indices], as.character)
    return(x)
    })
  return(datasetMerging(esets))
}

pooled.eset <- .getPooledEset(only.hgs = TRUE)
