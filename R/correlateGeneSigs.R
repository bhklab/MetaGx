#' Function to determine the correlation between genes in a signature
#'
#' This function returns a matrix containing the pairwise correlations for the genes in a signature. The correlations returned are meta estimates as
#' the correlation is computed between the signatures within each dataset and then a meta estimate is calculated by using the correlations of the signatures in each dataset.
#' @param sigInfoList The list returned from the function getGenesProgValue. 
#' @param method The method to be used to compute the correlation coefficient. Options are pearson (default), spearman, and kendall.
#' @param subtypeName Default is NULL, in which case all the patients in the esets are used. If a string corresponding to one of the subtypes from the analysis is supplied (subtypes found in names(sigInfoList$patientSurvData[[1]])) then only patients that are the subtype "subtypeName" will be used in the analysis.
#' @param genesRequired a number less than or equal to 1 and greater than 0 specifying how many genes must be present in a dataset from the signature in order for the scores of the patients
#' from that dataset to be used when calculating the correlation between 2 signatures
#' @return a matrix containing the pairwise correlations between the supplied gene signatures with the rows and columns named after the names of the geneSigList elements
#' @export
#' @examples
#'
#' geneSigList = list(c("KLK14","RHOX8","ADAMTS20","IDO1"), c("ADAMTS20","IDO1", "FAP"))
#' names(geneSigList) = c("4 Gene Signature", "3 Gene Signature")
#' geneDirecList = list(c(1, 1, -1, -1), c(-1, -1, 1))
#' sigInfoList = getGenesProgValue(geneSigList, geneDirecList, cancerType = "ovarian", subtype = "verhaak", survivalMetric = "overall")
#' corMatrix = correlateGeneSigs(sigInfoList)

correlateGeneSigs = function(sigInfoList, method = "pearson", subtypeName = NULL, genesRequired = 0.75)
{
  #tested 2 signatures that differed by just geneDirecList all +1 for one and all -1 for the other and got the expected -1 on diagonals
  #differs by small amount from old corSigs code (ex .873 to .870 ?)
  
  noDataInds = which(sapply(1:length(sigInfoList$patientSurvData), function(x) nrow(sigInfoList$patientSurvData[[x]][[1]])) == 0)
  if(length(noDataInds) > 0)
  {
    warning(paste(names(sigInfoList$patientSurvData)[noDataInds], "was not present in the data and is being removed"))
    sigInfoList$patientSurvData[noDataInds] = NULL
    sigInfoList$genesInfo[noDataInds] = NULL
  }
  
  #remove dependency on geneSigList
  #adjust to use results form getPatientSurvData
  if(is.null(subtypeName)){
    subtypeInd = 1
  }else{
    subtypeInd = which(tolower(names(sigInfoList$patientSurvData[[1]])) == tolower(subtypeName))
    if(length(subtypeInd) == 0)
      stop(paste(subtypeName, "Was not found as a subtype, valid subtypes are those from the analysis that produced sigInfoLists which are", names(sigInfoList$patientSurvData[[1]])))
  }
  
  if((genesRequired > 1) | (genesRequired <= 0))
    stop("Variable genesRequired must be greater than 0 or less than or equal to 1")
  
  nameVec = c()
  for(i in 1:length(sigInfoList$patientSurvData))
  {
    if(is.null(names(sigInfoList$patientSurvData))){
      nameVec = c(nameVec, paste("signature", i))
    }else{
      nameVec = c(nameVec, names(sigInfoList$patientSurvData)[i])
    }
  }
  #names(geneSigList) = nameVec
  
  expFrame = as.data.frame(NULL)
  numGenesPresList = list()

  sigLengths = unlist(sapply(1:length(sigInfoList$patientSurvData), function(x) nrow(sigInfoList$genesInfo[[x]])))
  corMat = matrix(rep(NA, length(sigInfoList$patientSurvData)*length(sigInfoList$patientSurvData)), nrow = length(sigInfoList$patientSurvData), ncol = length(sigInfoList$patientSurvData), byrow = TRUE)
  colnames(corMat) = names(sigInfoList$patientSurvData)[1:length(sigInfoList$patientSurvData)]
  rownames(corMat) = names(sigInfoList$patientSurvData)[1:length(sigInfoList$patientSurvData)]
  
  dataNames = as.character(unique(unlist(lapply(1:length(sigInfoList$patientSurvData), function(x) unique(sigInfoList$patientSurvData[[x]][[subtypeInd]]$dataName)))))
  corMatList = lapply(1:length(dataNames), function(x) corMat)
  names(corMatList) = dataNames
  standDevList = corMatList
  
  numGenesPresList = list()
  for(i in 1:length(sigInfoList$patientSurvData))
  {
    survData = sigInfoList$patientSurvData[[1]][[subtypeInd]]
    dataNames = as.character(unique(survData$dataName))
    genesPresVec = c()
    for(j in 1:length(dataNames))
    {
      firstPatientRow = which(survData$dataName == dataNames[j])[1]
      genesPresVec = c(genesPresVec, survData$numGenesPresent[firstPatientRow])
    }
    names(genesPresVec) = dataNames
    numGenesPresList[[i]] = genesPresVec
  }
  
  for(i in 1:length(sigInfoList$patientSurvData))
  {
    survData = sigInfoList$patientSurvData[[i]][[subtypeInd]]
    xPatNamesOrig = survData$`Patient ID`
    xPatDataNamesOrig  = survData$dataName
    remInds = c()
    numGenesPres = numGenesPresList[[i]]
    invalidDataNames = names(numGenesPres)[which(numGenesPres < genesRequired*sigLengths[i])]
    xOrig = survData$scoreVals
    if(length(invalidDataNames) > 0)
    {
      for(k in 1:length(invalidDataNames))
        remInds = c(remInds, which(survData$dataName == invalidDataNames[k]))
      xOrig = xOrig[-remInds]
      xPatNamesOrig = xPatNamesOrig[-remInds]
      xPatDataNamesOrig  = xPatDataNamesOrig[-remInds]
    }
    
    for(j in 1:length(sigInfoList$patientSurvData))
    {
      remInds = c()
      x = xOrig
      xPatNames = xPatNamesOrig
      xPatDataNames = xPatDataNamesOrig
      survData = sigInfoList$patientSurvData[[j]][[subtypeInd]]
      yPatNames = survData$`Patient ID`
      yPatDataNames  = survData$dataName
      numGenesPres = numGenesPresList[[j]]
      invalidDataNames = names(numGenesPres)[which(numGenesPres < genesRequired*sigLengths[j])]
      y = survData$scoreVals
      if(length(invalidDataNames) > 0)
      {
        for(k in 1:length(invalidDataNames))
          remInds = c(remInds, which(survData$dataName == invalidDataNames[k]))
        y = y[-remInds]
        yPatNames = yPatNames[-remInds]
        yPatDataNames = yPatDataNames[-remInds]
      }
      
      #make sure only same patients are used
      yKeep = c()
      xKeep = c()
      for(k in 1:length(xPatNames))
      {
        yKeep = c(yKeep, which(yPatNames == xPatNames[k]))
        if(length(which(yPatNames == xPatNames[k])) > 0)
          xKeep = c(xKeep, k)
      }
      y = y[yKeep]
      yPatNames = yPatNames[yKeep]
      yPatDataNames = yPatDataNames[yKeep]
      x = x[xKeep]
      xPatNames = xPatNames[xKeep]
      xPatDataNames = xPatDataNames[yKeep]
      
      for(k in 1:length(dataNames))
      {
        dataName = as.character(dataNames[k])
        keepInds = which(as.character(xPatDataNames) == dataName)
        if(length(keepInds) > 3)
        {
          r = cor(x[keepInds], y[keepInds], method = method)
          seFish = 1/sqrt(length(keepInds) - 3)
          standErr = tanh(seFish)
          corMatList[[k]][i, j] = r
          standDevList[[k]][i, j] = standErr
        }
      }
      
      #check = c()
      #for(k in 1:length(xKeep))
      #  check = c(check, xPatNames[k] == yPatNames[k])
      
      #corMat[i, j] = cor(x, y, method = method)
    }
  }
  
  for(i in 1:length(sigInfoList$patientSurvData))
  {
    for(j in 1:length(sigInfoList$patientSurvData))
    {
      cors = sapply(1:length(corMatList), function(x) corMatList[[x]][i, j])
      corsSe = sapply(1:length(standDevList), function(x) standDevList[[x]][i, j])
      corMat[i, j] = combine.est(cors[!is.na(cors)], corsSe[!is.na(corsSe)])$estimate
    }
  }
  
  return(corMat)
  
}

