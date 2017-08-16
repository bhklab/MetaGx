
#' Function to obtain survival info and risk predictions for patients
#'
#' This function returns a data frame with survival info and risk predictions for the patients and an additional data frame for each subtype
#' @param geneEntrezIds A numeric or character vector containing the entrez IDs for the genes in the signature being tested
#' @param geneDirecs A numeric vector composed of +1 and - 1 indicating the direction of association for each gene supplied in the geneEntrezIds vector.
#' +1 implies that the expression of that gene in a patient will be added to the patients score and -1 implies the expression will be subtracted from their score.
#' If one is looking for high scores to be associated with good survival, than genes that high expression is believed to lead to good prognosis should be given +1 and genes that high expression is believed
#' to lead to bad prognosis should be given a minus 1.
#' @param cancerType A string representing the type of cancer one would like to test the gene signature on. Options are currently "ovarian" and "breast"
#' @param subtype a string representing the subtyping scheme the patients will be classified under. For breast cancer, the options are "scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003", "intClust", "AIMS", or "claudinLow".
#' For ovarian cancer the options are "verhaak", "bentink", "tothill", "helland", and "konecny".
#' @param survivalMetric a string specifying the type of survival analysis that will be investigated in the reports. Either "overall" for overall survival or "relapse" for relapse free survival
#' @param dataList a list of esets to be used in the analysis. By default this variable is NULL and the other inputs will be used to obtain the esets for the analysis. However, it is useful to provide
#' a dataList if one is running this function multiple times so that the runtime will be reduced, or if there is a specific set of esets one would like to use in their analysis.
#' @param bestProbesList a list of vectors (one for each dataList element) containing the row indices of the probes in the eset@assayData$exprs expresson matrices that should be used. By default the function obtains dataLists
#' and the probes to be used on it own. If dataLists are provided and bestprobes is NULL than getBestprobes will be used to select the probes that have the highest IQR for each probe ID
#' @param removeMid a number greater than 0 and less than .5 specifying the fraction of the patients with scores/risk predictions in the middle of the patients scores to be removed when generating survival results
#' The default is 0, all patients used in the analysis
#' @param censorTime an integer specifying the point in time (years) at which the survival data must be censored. The default is 10 years
#' @return a list of data frames (1 for all patients and 1 for each subtype) containing the patients' risk prediction, time to the survival event, status of the survival event, the number of genes used in the analysis, the dataset the patient is from and the ID of the patient in the dataset
#' @export
#' @examples
#' geneEntrezIds = c("43847", "434768", "80070", "3620")
#' geneDirecs = c(1, 1, -1, -1)
#' survInfoList = getPatientSurvivalData(geneEntrezIds, geneDirecs, cancerType = "ovarian", subtype = "verhaak", survivalMetric = "overall")
#'

getPatientSurvivalData = function(geneEntrezIds, geneDirecs, cancerType, subtype, survivalMetric, dataList = NULL, bestProbesList = NULL, removeMid = 0, censorTime = 10)
{
  #i=3
  #geneEntrezIds = geneSigList[[i]]
  #geneEntrezIds = as.numeric(geneEntrezList[[i]])
  #geneDirecs = geneDirecList[[i]]
  #bestProbesList = dataListBestProbes
  geneEntrezIds = as.numeric(geneEntrezIds)
  
  if(removeMid > 0.5)
  {
    warning("removeMid was greater than 0.5 and has been changed to 0.5")
    removeMid = 0.5
  }
  #maybe/probably replace params in above with what is needed
  #add warning for genes that arent in the dataset
  survInfoList = list()

  naInds = which(is.na(geneEntrezIds))
  if(length(naInds) == length(geneEntrezIds)){
    survInfoList[[1]] = as.data.frame(NULL)
    names(survInfoList) = c("All Patients")
    return(survInfoList)
  }else if(length(naInds) > 0){
    geneEntrezIds = geneEntrezIds[-naInds]
    geneDirecs = geneDirecs[-naInds]
  }

  if(is.null(dataList) & is.null(bestProbesList)){
    print("bye")
    esetsAndProbes = getEsetsProbesSubtypesEvents(cancerType, survivalMetric, subtype)
    dataList = esetsAndProbes$esets
    dataListBestProbes = esetsAndProbes$esetsBestProbes
    infoString = names(dataList)
  }else if(!is.null(dataList) & is.null(bestProbesList)){
    dataListBestProbes = list()
    print("hi")
    for(i in 1:length(dataList))
      dataListBestProbes[[i]] = metaGx::getBestProbes(dataList[[i]])
  }else{
    dataListBestProbes = bestProbesList
  }

  geneSig = geneEntrezIds

  getSurvStatsGen = function(dataList, survivalMetric, subtypeName, censorTime)
  {
    survInfoDat = as.data.frame(NULL)

    for(i in 1:length(dataList))
    {
        survInfoNew = metaGx::getSurvInfo(dataList[[i]], geneSig, geneDirecs, survivalMetric, subtypeName = subtypeName, censorTime = censorTime, dataListBestProbes[[i]])
        if(!is.null(dim(survInfoNew)))
        {
          survInfoNew = cbind(survInfoNew, names(dataList)[i])
          colnames(survInfoNew)[ncol(survInfoNew)] = "dataName"
          survInfoDat = rbind(survInfoDat, survInfoNew)
        }
    }
    return(survInfoDat)
  }

  subtypeVec = c()
  for(i in 1:length(dataList))
    subtypeVec = c(subtypeVec, unique(as.character((dataList[[i]]@phenoData@data$subtypes))))
  subtypes = unique(subtypeVec)
  numSubtypes = length(subtypes)

  survStatsAll = getSurvStatsGen(dataList, survivalMetric, NULL, censorTime)
  if(is.null(survStatsAll) == FALSE)
  {
    survInfoList[[length(survInfoList) + 1]] = survStatsAll
    for(i in 1:numSubtypes)
    {
      survInfoList[[length(survInfoList) + 1]] = getSurvStatsGen(dataList, survivalMetric, subtypes[i], censorTime)
    }
    names(survInfoList) = c("All Patients", as.character(subtypes))
  }
  if(removeMid > 0)
    survInfoList = extremePatientSubset(survInfoList, removeMid)

  return(survInfoList)
}

#survInfoNew = metaGx::getSurvInfo(dataList[[1]], geneEntrezList[[1]], geneDirecList[[1]], survivalMetric, censorTime = 10, dataListBestProbes[[1]], subtypeName = NULL)
#geneData = testOvCancSurv(geneEntrezList[[i]][j], 1, params, bestProbes = geneInfo[[2]], dataList = geneInfo[[3]])
#geneEntrezIds = geneEntrezList[[1]]
#geneDirecs = geneDirecList[[1]]
#params = params
#bestProbes = dataListBestProbes
#dataList = dataList
#infoString = infoString
#geneData = testOvCancSurv(geneSigFrame[, 1], geneSigFrame[, 2], params, bestProbes = geneSigInfo[[2]], ovDataList = geneSigInfo[[3]], infoString = geneSigInfo[[4]])
