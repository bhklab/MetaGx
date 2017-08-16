
#' Function to test prognostic value of a signature relative to random signatures of the same size
#'
#' This function runs a permutation by assessing the prognostic value (via the D index) of signatures that are the same size as the one provided and determining how many of those random signatures are more prognostic than the one provided
#' @param geneIds A character vector containing the ensemble IDs, entrez IDs, or gene symbols for the genes in the signature being tested. Entrez IDs are 
#'  recommended as if gene symbols or ensemble Ids cant be mapped to entrez Ids then the gene will be omitted from the analysis
#' @param geneDirecs A numeric vector composed of +1 and - 1 indicating the direction of association for each gene supplied in the geneEntrezIds vector.
#' +1 implies that the expression of that gene in a patient will be added to the patients score and -1 implies the expression will be subtracted from their score.
#' If one is looking for high scores to be associated with good survival, than genes that high expression is believed to lead to good prognosis should be given +1 and genes that high expression is believed
#' to lead to bad prognosis should be given a minus 1.
#' @param cancerType A string representing the type of cancer one would like to test the gene signature on. Options are currently "ovarian" and "breast"
#' @param subtype a string representing the subtyping scheme the patients will be classified under. For breast cancer, the options are "scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003", "intClust", "AIMS", or "claudinLow".
#' For ovarian cancer the options are "verhaak", "bentink", "tothill", "helland", and "konecny".
#' @param numbPerms an integer specifying the number of random signatures of the same size to test the signature against
#' @param survivalMetric a string specifying the type of survival analysis that will be investigated in the reports. Either "overall" for overall survival or "relapse" for relapse free survival
#' @param dataList a list of esets to be used in the analysis. By default this variable is NULL and the other inputs will be used to obtain the esets for the analysis. However, it is useful to provide
#' a dataList if there is a specific set of esets one would like to use in their analysis.
#' @param bestProbesList a list of vectors (one for each dataList element) containing the row indices of the probes in the eset@assayData$exprs expresson matrices that should be used. By default the function obtains dataLists
#' and the probes to be used on it own. If dataLists are provided and bestprobes is NULL than getBestprobes will be used to select the probes that have the highest IQR for each probe ID
#' @param removeMid a number greater than 0 and less than .5 specifying the fraction of the patients with scores/risk predictions in the middle of the patients scores to be removed when generating survival results
#' The default is 0, all patients used in the analysis
#' @param censorTime an integer specifying the point in time (years) at which the survival data must be censored. The default is 10 years
#' @return a vector of p values for the permutation test indicating the chance of observing a signature as prognostically relevant as the one provided in all the patients and in each subtypes by chance. 
#' @export
#' @examples
#' geneIds = c("43847", "434768", "80070", "3620")
#' geneDirecs = c(1, 1, -1, -1)
#' pVec = permTestSig(geneEntrezIds, geneDirecs, cancerType = "ovarian", subtype = "verhaak", survivalMetric = "overall")
#'

permTestSig = function(geneIds, geneDirecs, cancerType, subtype, survivalMetric, numbPerms = 10000, dataList = NULL, bestProbesList = NULL, removeMid = 0, censorTime = 10)
{
  esetsAndProbes = getEsetsProbesSubtypes(cancerType, survivalMetric, subtype, dataNames = NULL)
  dataList = esetsAndProbes$esets
  dataListBestProbes = esetsAndProbes$esetsBestProbes
  infoString = names(dataList)
  
  geneInfoFrame = getGeneInfo(geneIds, geneDirecs)
  geneDataSig = getPatientSurvivalData(as.numeric(geneInfoFrame$`Entrez ID`), geneDirecs, cancerType, subtype, survivalMetric, dataList = dataList, bestProbesList = dataListBestProbes, removeMid = removeMid, censorTime = censorTime)
  #geneData
  for(i in 1:length(geneDataSig))
  {
    survInfo = geneDataSig[[i]]
    #statTabRow = names(survInfoList)[i]
    dindInfo = metaGx::getDindexOfDatasets(survInfo)
    dindInfo$dIndex = as.numeric(as.character(dindInfo$dIndex))
    dindInfo$dIndStandErr = as.numeric(as.character(dindInfo$dIndStandErr))
    dindComb = combine.est(dindInfo$dIndex[!is.na(dindInfo$dIndex)], dindInfo$dIndStandErr[!is.na(dindInfo$dIndStandErr)], hetero = TRUE)
    dIndex = dindComb$estimate
    dIndexSe = dindComb$se
    dindexP = 2*pnorm(-abs(log(dindComb$estimate)/dindComb$se))
  }
  
  sigResults = getMetaDindexAndLogP(geneDataSig, as.numeric(numGroups), numDigits = 5)
  
  dIndSigP = as.numeric(as.character(sigResults$`D Index P`))
  
  #some NAs introduced due to //
  #entList  = lapply(1:length(dataList), function(x) unique(as.numeric(as.character(dataList[[x]]@featureData@data$EntrezGene.ID))))
  #entIdsAll = Reduce(intersect, entList)
  
  entIdsAll  = unique(unlist(lapply(1:length(dataList), function(x) unique(as.numeric(as.character(dataList[[x]]@featureData@data$EntrezGene.ID))))))

  randSampList = lapply(1:numbPerms, function(x) sample(1:length(entIdsAll), length(geneIds)))
  geneEntrezList =  lapply(1:numbPerms, function(x) entIdsAll[randSampList[[x]]])
  
  geneDataList = list()
  dIndVec = c()
  for(i in 1:numbPerms)
  {
    geneData = getPatientSurvivalData(geneEntrezList[[i]], geneDirecs, cancerType, subtype, survivalMetric, dataList = dataList, bestProbesList = dataListBestProbes, removeMid = removeMid, censorTime = censorTime)
    dataResults = getMetaDindexAndLogP(geneData, as.numeric(numGroups), numDigits = 5)
    dIndDataP = as.numeric(as.character(dataResults$`D Index P`))
    dIndVec = c(dIndVec, dIndDataP)
  }
  subtypeNames = names(geneData)
  permTestPVec = c()
  for(i in 1:length(subtypeNames))
  {
    sigP = dIndSigP[i]
    dataP = dIndVec[seq(i, length(dIndVec), length(subtypeNames))]
    numbBetter = which(dataP < sigP)
    if(length(numbBetter) == 0){
      permTestP = 1/numbPerms
    }else{
      permTestP = length(numbBetter)/numbPerms
    }
    permTestPVec = c(permTestPVec, permTestP)
  }
  names(permTestPVec) = subtypeNames
  
  return(permTestPVec)
  
}


