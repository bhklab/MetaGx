
#' Function to determine the correlations between genes of patients
#'
#' This function returns a matrix containing the pairwise correlations for gene signatures according to the correlation between the signature
#' scores in the patients that would be used for a given analysis (cancerType and survivalMetric determine the patients used)
#' @param geneEntrezIds A character vector containing the entrez IDs for the genes in the signature being tested
#' @param cancerType a string specifying what cancer type to use in the analysis. The options are "breast" and "ovarian"
#' @param survivalMetric a string specifying what type of survival data to use in the analysis. If the cancerType is set to "ovarian"
#' then enter "overall" for overall survival and "relapse" for relapse free survival. If the cancerType is set to "breast" then the options
#' are "relapse" for relapse free survival, "overallMeta" for overall survival on the metabric study data, "overallTcga" for
#' overall survival on the TCGA study data, and "overall" for overall survival on the remaining studies (No Metabric and TCGA data).
#' @param method The method to be used to compute the correlation coefficient. Options are pearson (default), spearman, and kendall.
#' @param dataList a list of esets to be used in the analysis. By default this variable is NULL and the other inputs will be used to obtain the esets for the analysis. However, it is useful to provide
#' a dataList if one is running this function multiple times so that the runtime will be reduced, or if there is a specific set of esets one would like to use in their analysis.
#' @return a matrix containing the pairwise correlations between the supplied genes with the rows and columns named after the entrez gene ids
#' @export
#' @examples
#'
#' geneEntrezIds = c("43847", "434768", "80070", "3620")
#' corMatrix = correlateGenes(geneEntrezIds, cancerType = "ovarian", survivalMetric = "overall")

correlateGenes = function(geneEntrezIds, cancerType, survivalMetric, method = "pearson", dataList = NULL)
{
  naGenes = which(is.na(geneEntrezIds))
  if(length(naGenes) > 0)
    stop("You have supplied NA entrez gene IDs. Please remove the missing genes and rerun the function")
  if(is.null(dataList))
    dataList = loadMetaData(cancerType, survivalMetric)
  numPatientsVec = sapply(1:length(dataList), function(x) ncol(dataList[[x]]@assayData$exprs))
  numPatients = sum(numPatientsVec)
  dataNameVec = unlist(sapply(1:length(dataList), function(x) rep(names(dataList)[x], numPatientsVec[x])))
  expFrame = as.data.frame(NULL)
  for(i in 1:length(geneEntrezIds))
  {
    entrez = geneEntrezIds[i]
    expVec = c()
    for(j in 1:length(dataList))
    {
      entrezIds = dataList[[j]]@featureData@data$EntrezGene.ID
      rowInd = which(entrezIds == entrez)
      if(length(rowInd) == 0){
        expVec = c(expVec, rep(NA, numPatientsVec[j]))
      }else{
        expVec = c(expVec, dataList[[j]]@assayData$exprs[rowInd, ])
      }
    }
    expFrame = rbind(expFrame, expVec)
  }
  rownames(expFrame) = geneEntrezIds
  missingExp = as.vector(colSums(apply(expFrame, 1, is.na)))
  notPresent = which(missingExp == numPatients)

  if(length(notPresent) > 0)
    expFrame = expFrame[-notPresent, ]

  corMat = matrix(rep(NA, nrow(expFrame)*nrow(expFrame)), nrow = nrow(expFrame), ncol = nrow(expFrame), byrow = TRUE)
  colnames(corMat) = rownames(expFrame)
  rownames(corMat) = rownames(expFrame)
  
  dataNames = names(dataList)
  corMatList = lapply(1:length(dataNames), function(x) corMat)
  names(corMatList) = dataNames
  standDevList = corMatList
  
  for(i in 1:nrow(expFrame))
  {
    for(j in 1:nrow(expFrame))
    {
      x = expFrame[i, ]
      names(x) = dataNameVec
      y = expFrame[j ,]
      names(y) = dataNameVec
      missingVec = c(which(is.na(expFrame[i, ])), which(is.na(expFrame[j, ])))
      missingVec = unique(missingVec)
      if(length(missingVec) > 0)
      {
        x = x[-missingVec]
        y = y[-missingVec]
      }
      #corMat[i, j] = cor(t(x), t(y), method = method)
      
      for(k in 1:length(dataNames))
      {
        dataName = as.character(dataNames[k])
        keepInds = which(grepl(dataName, names(x)))
        if(length(keepInds) > 3)
        {
          #print(length(keepInds))
          r = cor(t(x[keepInds]), t(y[keepInds]), method = method)
          seFish = 1/sqrt(length(keepInds) - 3)
          standErr = tanh(seFish)
          corMatList[[k]][i, j] = r
          standDevList[[k]][i, j] = standErr
        }
      }
    }
  }
  
  for(i in 1:nrow(expFrame))
  {
    for(j in 1:nrow(expFrame))
    {
      cors = sapply(1:length(corMatList), function(x) corMatList[[x]][i, j])
      corsSe = sapply(1:length(standDevList), function(x) standDevList[[x]][i, j])
      corMat[i, j] = combine.est(cors[!is.na(cors)], corsSe[!is.na(corsSe)])$estimate
    }
  }
  
  
  return(corMat)

}

