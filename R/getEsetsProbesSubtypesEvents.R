#' Function to obtain a list of expression sets (eSets), their best probes, their patients subtypes, and the event time and status of the patients
#'
#' This function returns a list of expression sets (eSets), a list of each eSets best probes (see getBestProbes for the nature of the results), 
#' a list of their event times (see getSurvEvent), a list of their event statuses (for the given survival metric specified), and adds a subtype element to the eSets (see getPatientSubtypes).
#' 
#' @param cancerType A string representing the type of cancer one would like to test the gene signature on. Options are currently "ovarian" and "breast"
#' @param survivalMetric a string specifying the type of survival analysis that will be investigated in the reports. Either "overall" for overall survival or "relapse" for relapse free survival
#' @param subtype a string representing the subtyping scheme the patients will be classified under. For breast cancer, the options are "scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003", "intClust", "AIMS", or "claudinLow".
#' For ovarian cancer the options are "verhaak", "bentink", "tothill", "helland", and "konecny".
#' @param dataNames a character vector specifying the names of the datasets to include in the analysis. The names must match the names in the column "Data Name" from the dataframe returned from the
#' obtainDataInfo function. Use loadMetaData followed by obtainDataInfo with the cancerType and survivalMetric of interest to get the appropriate data names in the obtainDataInfo table for your analysis. By default
#' all datasets are included.
#' @param remMissSubtypeData A boolean specifying whether datasets that do not contain patients with all of the subtypes in the subtyping scheme specified should be removed. TRUE by default and TRUE is required in order for the 
#' results to be used as inputs to the getPatientSurvivalData function when running a survival analysis on a given signatur eor gene
#' @param includeAll a boolean specifying whether to include the TCGA and METABRIC datasets in the breast cancer analysis. Default is FALSE
#' @return a data frame with information about the esets supplied
#' @export
#' @examples
#' 
#' esetsProbesEvents = getEsetsProbesSubtypesEvent("ovarian", "overall", "verhaak")
#' 

getEsetsProbesSubtypesEvents = function(cancerType, survivalMetric, subtype, dataNames = NULL, remMissSubtypeData = TRUE, includeAll = FALSE)
{
  
  dataList = loadMetaData(cancerType, survivalMetric, dataNames, includeAll)
  dataListBestProbes = list()
  
  for(i in 1:length(dataList))
    dataListBestProbes[[i]] = metaGx::getBestProbes(dataList[[i]])
  
  for(i in 1:length(dataList))
    dataList[[i]] = getPatientSubtypes(dataList[[i]], cancerType, subtype, intersectThresh = 0.75)
  
  missSubInds = which(sapply(1:length(dataList), function(x) sum(is.na(dataList[[x]]$subtypes))) == sapply(1:length(dataList), function(x) length(dataList[[x]]$subtypes)))
  dataList[missSubInds] = NULL
  dataListBestProbes[missSubInds] = NULL
  
  survEventList = getSurvEventData(dataList, survivalMetric)
  
  subtypeVec = c()
  for(i in 1:length(dataList))
    subtypeVec = c(subtypeVec, unique(as.character((dataList[[i]]@phenoData@data$subtypes))))
  subtypes = unique(subtypeVec)
  #in case of NAs from inability to map patient to subtype
  subtypes = subtypes[!is.na(subtypes)]
  numSubtypes = length(subtypes)
  allSubNotPresent = which(sapply(dataList, function(x) length(unique(as.vector(x$subtypes))[!is.na(unique(as.vector(x$subtypes)))])) != numSubtypes)
  
  if(length(allSubNotPresent) > 0)
  {
    dataList[allSubNotPresent] = NULL
    dataListBestProbes[allSubNotPresent] = NULL
  }
  names(dataListBestProbes) = names(dataList)
  esetsProbesEvents = list(dataList, dataListBestProbes, survEventList$eventList, survEventList$eventTimeList)
  names(esetsProbesEvents) = c("esets", "esetsBestProbes", "esetsEventStatus", "esetsEventTimes")
  return(esetsProbesEvents)
}
