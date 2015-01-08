########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

`checkNames` <- 
function (eset) {
  ## check feature and sample names in esets
  #
  # Args:
  #   esets: expressionSet object
  #
  # Returns:
  #   esets: updated expressionSet object with atching feature and sample names
  
  ## check feature names
  check.feature <- intersect(rownames(Biobase::exprs(eset)), rownames(Biobase::fData(eset)))
  if (length(check.feature) == 0) {
    warning("Names of features do not match between expressions and annotations")
    return (NULL)
  } else {
    if (length(check.feature) != nrow(Biobase::exprs(eset)) || length(check.feature) != nrow(Biobase::fData(eset))) {
      warning("Some features are missing between expressions and annotations")
    }
  }
  ## check sample names
  check.sample <- intersect(colnames(Biobase::exprs(eset)), rownames(Biobase::pData(eset)))
  if (length(check.sample) == 0) {
    warning("Names of samples do not match between expressions and phenotypes")
    return (NULL)
  } else {
    if (length(check.sample) != ncol(Biobase::exprs(eset)) || length(check.sample) != nrow(Biobase::pData(eset))) {
      warning("Some samples are missing between expressions and phenotypes")
    }
  }
  Biobase::exprs(eset) <- Biobase::exprs(eset)[check.feature, check.sample, drop=FALSE]
  Biobase::fData(eset) <- Biobase::fData(eset)[check.feature, , drop=FALSE]
  Biobase::pData(eset) <- Biobase::pData(eset)[check.sample, , drop=FALSE]
  Biobase::pData(eset)[ , "samplename"] <- rownames(Biobase::pData(eset))
  return (eset)
}

