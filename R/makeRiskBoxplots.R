
#' Function to create a boxplot showing the risk predictions for the patients in each subtype
#'
#' This function returns a data frame which describes the number of events (deceased or tumor recurrence), median time to an event, the number of patients, the number of genes, and the platform of an eset
#' @param survInfoList The list returned from the function getPatientSurvivalData. Alternatively, a list of data frames with 5 columns called eventTime, eventStatus, scoreVals, scoreVals, numGenesPresent, and dataName that correspond to
#' the survival event times, the survival event's status, the score/risk prediction of the patient, the number of genes used to calculate the score and the name of the dataset that the patient is from. 
#' The names of the dataframes in the list will be used for the boxplot names
#' @return a boxplot with a p value from the Kruskal-Wallis rank sum test between the scores/risk predictions in the data frames of survInfoList
#' @export
#' @examples
#' 
#' geneEntrezIds = c("43847", "434768", "80070", "3620")
#' geneDirecs = c(1, 1, -1, -1)
#' survInfoList = getPatientSurvivalData(geneEntrezIds, geneDirecs, cancerType = "ovarian", subtype = "verhaak", survivalMetric = "overall")
#' makeRiskBoxplots(survInfoList)
#' 

makeRiskBoxplots = function(survInfoList)
{
  groupNames = names(survInfoList)
  if(is.null(groupNames))
    groupNames = sapply(1:length(survInfoList), function(x) paste("Group", x))
  scoreVec = c()
  groupVec = c()
  for(i in 1:length(survInfoList))
  {
    scoreVec = c(scoreVec, survInfoList[[i]]$scoreVals)
    groupVec = c(groupVec, rep(groupNames[i], length(survInfoList[[i]]$scoreVals)))
  }
  
  kruskalP = kruskal.test(scoreVec ~ as.factor(groupVec))$p.value
  boxplot(scoreVec ~ groupVec, ylab = "Risk Predictions", cex.axis = 1.3, cex.lab = 1.3)
  text(length(survInfoList) - 1 , min(scoreVec), paste("Kruskal-Wallis Test P =", format(kruskalP, scientific =  TRUE, digits = 2)), cex = 1.4)
  
}



