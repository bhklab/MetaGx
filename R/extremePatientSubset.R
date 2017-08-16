
#' Function to remove risk predictions from a list of survival info frames
#'
#' This function separates the patients into their respective datasets and then computes the D index for each dataset using the risk predictions, event times, and event occurence indicators of the patients
#' @param survInfo One of the data frames in the list returned from the function getPatientSurvivalData. Alternatively, a data frame with 5 columns called eventTime, eventStatus, scoreVals, scoreVals, numGenesPresent, and dataName that correspond to
#' the survival event times, the survival event's status, the score/risk prediction of the patient, the number of genes used to calculate the score and the name of the dataset that the patient is from. 
#' @param removeMid a number greater than 0 and less than 1 specifying the fraction of the patients with scores/risk predictions in the middle of the patients scores to be removed when generating survival results
#' The default is 0, all patients used in the analysis
#' @return a data frame with information about the datasets, including the D index and standard error of the D index for each dataset
#' @export
#' @examples
#' 
#' geneEntrezIds = c("43847", "434768", "80070", "3620")
#' geneDirecs = c(1, 1, -1, -1)
#' survInfoList = getPatientSurvivalData(geneEntrezIds, geneDirecs, cancerType = "ovarian", subtype = "verhaak", survivalMetric = "overall")
#' survInfoList = extremePatientSubset(survInfoList, 0.33)
#' 

extremePatientSubset = function(survInfoList, removeMid)
{
  if(removeMid > 0 & removeMid < 1){
    keepTopBotPerc = (1 - removeMid)/2
    survInfoListNew = list()
    for(i in 1:length(survInfoList))
    {
      #testVec = c()
      survInfoNew = as.data.frame(NULL)
      survInfo = survInfoList[[i]]
      if(nrow(survInfo) > 0){
        dataNames = as.character(unique(survInfo$dataName))
        rownames(survInfo) = as.character(1:nrow(survInfo))
        for(j in 1:length(dataNames))
        {
          dataStr = dataNames[j]
          dataInds = which(as.character(survInfo$dataName) == dataStr)
          scoreVals = survInfo$scoreVals[dataInds]
          survInfoTemp = survInfo[dataInds, ]
          negToPosInds = order(scoreVals)
          numbKeepTopBot = floor(length(scoreVals)*keepTopBotPerc)
          keepInds = c(negToPosInds[1:numbKeepTopBot], negToPosInds[(length(negToPosInds) - numbKeepTopBot):length(negToPosInds)])
          survInfoNew = rbind(survInfoNew, survInfo[as.numeric(rownames(survInfoTemp)[keepInds]), ])
          #testVec = c(testVec, as.numeric(rownames(survInfoTemp)[keepInds]))
        }
        survInfoListNew[[i]] = survInfoNew
      }else{
        survInfoListNew[[i]] = survInfo
      }
    }
    names(survInfoListNew) = names(survInfoList)
    return(survInfoListNew)
    
      }else{
    stop("removeMid must be greater than 0 and less than 1")
  }
}
