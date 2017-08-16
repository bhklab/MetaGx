
#' Function to compute the D index for risk predictions (and their standard errors) on each dataset that was used in the survival analysis
#'
#' This function separates the patients into their respective datasets and then computes the D index for each dataset using the risk predictions, event times, and event occurence indicators of the patients
#' @param survInfo One of the data frames in the list returned from the function getPatientSurvivalData. Alternatively, a data frame with 7 columns called eventTime, eventStatus, scoreVals, scoreVals, numGenesPresent, dataName and patient ID that correspond to
#' the survival event times, the survival event's status, the score/risk prediction of the patient, the number of genes used to calculate the score, teh ID of the patient in the dataset and the name of the dataset that the patient is from. 
#' @param minPatients an integer specifying the minimum number of patients/samples required to calculate a D.index. Default is 10
#' @return a data frame with information about the datasets, including the D index and standard error of the D index for each dataset
#' @export
#' @examples
#' 
#' geneEntrezIds = c("43847", "434768", "80070", "3620")
#' geneDirecs = c(1, 1, -1, -1)
#' survInfoList = getPatientSurvivalData(geneEntrezIds, geneDirecs, cancerType = "ovarian", subtype = "verhaak", survivalMetric = "overall")
#' dIndexData = getDindexOfDatasets(survInfoList$IMR)
#' dIndexData

getDindexOfDatasets = function(survInfo, minPatients = 10)
{
  forestPlotFrame = as.data.frame(NULL)
  dataNames = as.character(unique(survInfo$dataName))
  dindVec = c()
  dindSeVec = c()
  dindLowVec = c()
  dindHighVec = c()
  sampleString = c()
  numGenesVec = c()
  for(i in 1:length(dataNames))
  {
    dataStr = dataNames[i]
    dataRows =  which(as.character(survInfo$dataName) == dataStr)
    sampleString = c(sampleString, length(dataRows))
    numGenesVec = c(numGenesVec, survInfo[dataRows, ]$numGenesPresent[1])
    #seen d indices of 300 with 8 samples, so set sample min requirement?
    if(length(dataRows) > minPatients){
      survInfoInd = survInfo[dataRows, ]
      #Note: whether you use the original scores or shifted scores, D.index returns the same results
      dIndInfo = D.index(x=survInfoInd$scoreValsOrig, surv.time=survInfoInd$timeToDeath, surv.event=survInfoInd$vitalStat)
      dInds = dIndInfo$d.index
      dIndSe = dIndInfo$se
      dIndsLow = dIndInfo$lower
      dIndsHigh = dIndInfo$upper
      #dIndsLow = dInds*exp(-qnorm(.05, lower.tail = FALSE)*dIndSe)
      #dIndsHigh = dInds*exp(qnorm(.05, lower.tail = FALSE)*dIndSe)
    }else{
      dInds = NA
      dIndsLow = NA
      dIndsHigh = NA
      dIndSe = NA
    }
    dindVec = c(dindVec, dInds)
    dindSeVec = c(dindSeVec, dIndSe)
    dindLowVec = c(dindLowVec, dIndsLow)
    dindHighVec = c(dindHighVec, dIndsHigh)
  }
  forestPlotFrame = as.data.frame(cbind(dataNames, sampleString, numGenesVec, dindVec, dindSeVec, dindLowVec, dindHighVec))
  colnames(forestPlotFrame) = c("dataset", "numbPatients", "numGenesInData", "dIndex", "dIndStandErr", "dIndLowerBound", "dIndUpperBound")
  
  return(forestPlotFrame)
}

