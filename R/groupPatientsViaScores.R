
#' Function to group patients according to their risk predictions
#'
#' This function groups patients by splitting their risk predictions/scores into a number of equally sized groups according to the number of groups specified by the user. The grouping is most commonly performed in order to conduct survival analyses on the patients
#' @param survInfo One of the data frames in the list returned from the function getPatientSurvivalData. Alternatively, a data frame with 5 columns called eventTime, eventStatus, scoreVals, scoreVals, numGenesPresent, and dataName that correspond to
#' the survival event times, the survival event's status, the score/risk prediction of the patient, the number of genes used to calculate the score and the name of the dataset that the patient is from. 
#' @param numGroups an integer specifying the number of groups the patients will be split into when generating survival curves and calculating the log rank p value. The value is 2 by default and cannot exceed 5.
#' @param normalizeEsetScores A boolean specifying whether the groups used to determine the log rank p value should be determined by 1. Splitting the scores/risk predictions of the patients from each dataset into the number of groups specified or
#' 2. normalizing the scores of the patients in each dataset such that the 2.5 percentile and 97.5 percentile for the patients scores in each dataset have the same value, then creating one vector of patient score and splitting this
#' vector of patient score sinto the number of groups specified (TRUE). Normalization/TRUE is the default.
#' @return a vector the length of the number of patients in the survInfo data frame that has the group (integer) that each patient is part of
#' @export
#' @examples
#' 
#' geneEntrezIds = c("43847", "434768", "80070", "3620")
#' geneDirecs = c(1, 1, -1, -1)
#' survInfoList = getPatientSurvivalData(geneEntrezIds, geneDirecs, cancerType = "ovarian", subtype = "verhaak", survivalMetric = "overall")
#' patientGroups = groupPatientsViaScores(survInfo = survInfoList$`All Patients`, numGroups = 2, normalizeEsetScores = TRUE)

groupPatientsViaScores = function(survInfo, numGroups, normalizeEsetScores = TRUE)
{
  getGroups = function(scoreVals, numGroups)
  {
    sortInd = sort(scoreVals, decreasing = TRUE, index.return=TRUE)$ix;
    groupSize = floor(length(scoreVals)/numGroups)
    remainder =  length(scoreVals) - groupSize*numGroups
    dataGroups = matrix(0, length(scoreVals))
    extra = 0
    for(i in 0:(numGroups - 1))
    {
      startInd = 1 + i*groupSize - extra
      if((numGroups - i) == remainder)
      {
        groupSize = groupSize + 1
        extra = i
      }
      endInd = (i + 1)*groupSize - extra
      dataGroups[sortInd[startInd:endInd]] = i;
    } 
    return(dataGroups)
  }
  
  dataGroups = c()
  if(normalizeEsetScores == TRUE)
  {
    scoreVals = survInfo$scoreVals
    dataGroups = getGroups(scoreVals, numGroups)
  }else{
    dataNames = as.character(unique(survInfo$dataName))
    dataOrd = c()
    for(i in 1:length(dataNames))
    {
      dataStr = dataNames[i]
      dataInds = which(as.character(survInfo$dataName) == dataStr)
      scoreVals = survInfo$scoreVals[dataInds]
      dataGroups = c(dataGroups, getGroups(scoreVals, numGroups))
      dataOrd = c(dataOrd, dataInds)
    }
    dataGroups = dataGroups[dataOrd]
  }
  return(dataGroups)
}

