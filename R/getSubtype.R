########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

## FIXME 
## the package should properly extend the expressionSet class to add a slot for subtyping information
## the current workaround is using the experimentData slot
`getSubtype` <- 
function (eset, method=c("class", "crisp", "fuzzy")) {
  method <- match.arg(method)
  if (class(eset) != "newEset") {
    stop("eset should be an expressionSet object")
  }
  # if (length(Biobase::experimentData(eset)@other) == 0) { return (NULL) }
  switch(method,
    "class" = {
      res <- eset@subtype
    },
    "crisp" = {
      res <- eset@crisp
    },
    "fuzzy" = {
      res <- eset@fuzzy
    }  
  )
  return (res)
}


