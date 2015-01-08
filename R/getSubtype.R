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
  if (class(eset) != "ExpressionSet") {
    stop("eset should be an expressionSet object")
  }
  if (length(Biobase::experimentData(eset)@other) == 0) { return (NULL) }
  switch(method,
    "class" = {
      res <- Biobase::experimentData(eset)@other$class
    },
    "crisp" = {
      res <- Biobase::experimentData(eset)@other$crisp
    },
    "fuzzy" = {
      res <- Biobase::experimentData(eset)@other$fuzzy
    }  
  )
  return (res)
}


