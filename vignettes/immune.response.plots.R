source("../R/datasetMerging.R")

getPooledEset(only.hgs=TRUE) {
  return(datasetMerging(esets))
}