########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

## FIXME 
## the package should properly extend the expressionSet class to add a slot for subtyping information
## the current workaround is using the experimentData slot
`setSubtype` <- 
function (eset, subtype.class, subtype.crisp, subtype.fuzzy) {
  if (class(eset) != "newEset") {
    stop("eset should be an expressionSet object")
  }
  ## subtype.class is a vector of subtype labels, each sample is associated to a single subtype
  if (!missing(subtype.class)) {
    if (!is.factor(subtype.class)) { stop("subtype.class must be a vector of factors") }
    if (length(subtype.class) != length(Biobase::sampleNames(eset))) { stop("Length of subtype.class must be equal to the number of samples in eset") }
    if (!all(names(subtype.class) %in% Biobase::sampleNames(eset))) { stop("Names of subtype.class must be the same than sample names in eset") }
    subtype.class <- subtype.class[Biobase::sampleNames(eset)]
  } else {
    subtype.class <- NULL
  }
  ## subtype.crisp is a matrix of subtype discrete calls
  if (!missing(subtype.crisp)) {
    if (nrow(subtype.crisp) != length(Biobase::sampleNames(eset))) { stop("Number of rows of subtype.crisp must be equal to the number of samples in eset") }
    if (!all(rownames(subtype.crisp) %in% Biobase::sampleNames(eset))) { stop("Row names of subtype.crisp must be the same than sample names in eset") }
    subtype.crisp <- subtype.crisp[Biobase::sampleNames(eset), , drop=FALSE]
  } else {
    subtype.crisp <- NULL
  }
  ## subtype.fuzzy is a matrix of subtype discrete calls
  if (!missing(subtype.fuzzy)) {
    if (nrow(subtype.fuzzy) != length(Biobase::sampleNames(eset))) { stop("Number of rows of subtype.fuzzy must be equal to the number of samples in eset") }
    if (!all(rownames(subtype.fuzzy) %in% Biobase::sampleNames(eset))) { stop("Row names of subtype.fuzzy must be the same than sample names in eset") }
    subtype.fuzzy <- subtype.fuzzy[Biobase::sampleNames(eset), , drop=FALSE]
  } else {
    subtype.fuzzy <- NULL
  }
  eset@subtype <- subtype.class
  eset@crisp <- subtype.crisp
  eset@fuzzy <- subtype.fuzzy
  
  return (eset)
}


