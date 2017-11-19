#' Function to obtain cancer gene signatures from literature
#'
#' This function returns known cancer gene signatures
#' @param cancerType A string specifying what cancer to obtain signatures for. Options are currently "ovarian" and "breast"
#' @return A list containing 2 lists of vectors. List 1 corresponds to the gene signatures and list 2 corresponds to the direction of association
#' for each gene in the gene signatures
#' @export
#' @examples
#'
#' cancerSigs = obtainCancerSigs("ovarian")


obtainCancerSigs = function(cancerType)
{
  if(cancerType == "ovarian"){
    data("ovarianSigs")
    benchList = ovarianSigs
  }else if(cancerType == "breast"){
    data("breastSigs")
    benchList = breastSigs
  }
  geneSigList = list()
  geneDirecList = list()
  for(i in 1:length(benchList))
  {
    entrezIds = as.numeric(as.character(benchList[[i]]$EntrezGene.ID))
    geneDirecs = as.numeric(as.character(benchList[[i]]$coefficient))/abs(as.numeric(as.character(benchList[[i]]$coefficient)))
    #removes first gene which according to paper it appears isnt used anyways
    remInds = which(is.na(geneDirecs))
    if(length(remInds) > 0)
    {
      entrezIds = entrezIds[-remInds]
      geneDirecs = geneDirecs[-remInds]
    }
    geneSigList[[i]] = entrezIds
    geneDirecList[[i]] = geneDirecs
  }
  names(geneSigList) = names(benchList)
  names(geneDirecList) = names(benchList)

  cancerSigs = list(geneSigList, geneDirecList)
  names(cancerSigs) = c("geneSigList", "geneDirecList")
  return(cancerSigs)

}
