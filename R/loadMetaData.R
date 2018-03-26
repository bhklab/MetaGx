#' Function to load data required for the survival reports
#'
#' This function returns the data required to generate a survival report according to the user inputs that define the cancer type and survival type being considered
#' @param cancerType a string specifying what cancer type to use in the analysis. The options are "breast" and "ovarian"
#' @param survivalMetric a string specifying what type of survival data to use in the analysis. If the cancerType is set to "ovarian"
#' then enter "overall" for overall survival and "relapse" for relapse free survival. If the cancerType is set to "breast" then the options
#' are "relapse" for relapse free survival, "overallMeta" for overall survival on the metabric study data, "overallTcga" for
#' overall survival on the TCGA study data, "overall" for overall survival on the remaining studies (No Metabric and TCGA data), and metastasis for distant metastasis free survival.
#' Note that for all cancer types the option "hierarchy" is available and this survivalMetric parameter checks for relapse free, then distant metastasis, and then overall survival data
#' in each dataset and keeps the first survival data it finds. Note that the ovarian datasets do not have distant metastasis free survival info and for the breast analyses metabric and tcga are removed.
#' @param dataNames a character vector specifying the names of the datasets to include in the analysis. The names must match the names in the column "Data Name" from the dataframe returned from the
#' obtainDataInfo function. Use loadMetaData followed by obtainDataInfo with the cancerType and survivalMetric of interest to get the appropriate data names in the obtainDataInfo table for your analysis. By default
#' all datasets are included.
#' @param includeAll a boolean specifying whether to include the TCGA and METABRIC datasets in the breast cancer analysis. Default is FALSE
#' @return a list containing the datasets that will be used in the analysis
#' @export
#' @examples
#'
#' dataList = loadMetaData(cancerType = "ovarian", survivalMetric = "hierarchy")

loadMetaData = function(cancerType, survivalMetric, dataNames = NULL, includeAll = FALSE)
{

  dataList = list()
  package.name = NULL

  if(cancerType == "ovarian")
  {
    dataListAndDups = loadOvarianEsets()
    dataList = dataListAndDups$esets
  }else if(cancerType == "breast"){
    dataListAndDups = loadBreastEsets()
    dataListTcgaMetaAndDups = loadBreastEsets(loadString = c("TCGA", "METABRIC"))
    dataList = c(dataListAndDups$esets, dataListTcgaMetaAndDups$esets)
  }
  infoString = names(dataList)

  if(!is.null(dataNames))
  {
    keepInds = c()
    for(i in 1:length(dataNames))
      keepInds = c(keepInds, which(names(dataList) == dataNames[i]))
    if(length(keepInds) == 0){
      stop(cat("None of the data names supplied were present for the given cancerType and survivalMetric. Please run 1. dataList = loadMetaData(\"",cancerType, "\",\"",survivalMetric,"\") 2. validNames = obtainDataInfo(dataList, \"", survivalMetric,"\")$`Data Name` to see valid data names", sep =""))
    }else{
      remove = c(1:length(dataList))[-keepInds]
      dataList[remove] = NULL
      infoString = infoString[-remove]
    }
  }

  survEventList = getSurvEventData(dataList, survivalMetric)
  statList = survEventList$eventList
  timeList = survEventList$eventTimeList

  remove = c()
  #remove when all vital/tumor reccurence statuses are NA or when all daysToDeath/daysToReccurrence are NA
  for(i in 1:length(dataList))
  {
    numMissing = sum(is.na(timeList[[i]]))
    numPresent = length(timeList[[i]])
    if(numMissing == numPresent)
      remove = c(remove, i)
    numMissing = sum(is.na(statList[[i]]))
    numPresent = length(statList[[i]])
    if(numMissing == numPresent)
      remove = c(remove, i)
  }
  if(length(remove) > 0)
  {
    remove = unique(remove)
    dataList[remove] = NULL
    infoString = infoString[-remove]
    names(dataList) = infoString
  }

  #Metabric and TCGA studies are very large so offer option to do them separately. Also would skew results if used them with other data
  if(cancerType == "breast")
  {
    if(survivalMetric == "overall" | survivalMetric == "hierarchy"){
      if(includeAll == FALSE)
      {
        metaTcgaInds = c(which(infoString == "METABRIC"), which(infoString == "TCGA"))
        if(length(metaTcgaInds) > 0){
          dataList[metaTcgaInds] = NULL
          infoString = infoString[-metaTcgaInds]
        }
      }
    }else if(survivalMetric == "overallTcga"){
      dataList = list(TCGA)
      infoString = c("TCGA")
    }else if(survivalMetric == "overallMeta"){
      dataList = list(METABRIC)
      infoString = c("METABRIC")
    }
  }
  names(dataList) = infoString

  return(dataList)
}
