#' Function to obtain information about various esets
#'
#' This function returns a data frame which describes the number of events (deceased or tumor recurrence), median time to an event, the number of patients, the number of genes, and the platform of an eset
#' @param survInfoList The list returned from the function getPatientSurvivalData. Alternatively, a list of data frames with 5 columns called eventTime, eventStatus, scoreVals, scoreVals, numGenesPresent, and dataName that correspond to
#' the survival event times, the survival event's status, the score/risk prediction of the patient, the number of genes used to calculate the score and the name of the dataset that the patient is from. 
#' an error will be thrown if the esets in dataList do not contain survival information for the requested event. Enter "overall" for
#' overall survival and "relapse" for relapse free survival.
#' @param numGroups an integer specifying the number of groups the patients will be split into when generating survival curves and calculating the log rank p value. The value is 2 by default and cannot exceed 5.
#' @param normalizeEsetScores A boolean specifying whether the groups used to determine the log rank p value should be determined by 1. Splitting the scores/risk predictions of the patients from each dataset in the number of groups specified or
#' 2. normalizing the scores of the patients in each dataset such that the 2.5 percentile and 97.5 percentile for the patients scores in each dataset have the same value, then creating one vector of patient score and splitting this
#' vector of patient scores into the number of groups specified (TRUE). Normalization/TRUE is the default.
#' @param numDigits an integer specifying the number of digits for the values in the dataframe, Default is 2
#' @param minPatients an integer specifying the minimum number of patients/samples required to calculate a D.index for one of the datasets. Default is 20. Note that for a dataset with 100 patients and 4 subtypes, something much larger than 20 would 
#' likely remove the D.index in some of the subtype analyses (samples in subtype analysis ~ 100 patients divide 4 subtypes)
#' @return a data frame with the D index, D index confidence interval, D index p value, and log rank p value for each survival frame in the survInfoList
#' @export
#' @examples
#' 
#' geneEntrezIds = c("43847", "434768", "80070", "3620")
#' geneDirecs = c(1, 1, -1, -1)
#' survInfoList = getPatientSurvivalData(geneEntrezIds, geneDirecs, cancerType = "ovarian", subtype = "verhaak", survivalMetric = "overall")
#' resultsTable = getMetaDindexAndLogP(survInfoList, 2)
#' resultsTable

getMetaDindexAndLogP = function(survInfoList, numGroups, normalizeEsetScores = FALSE, numDigits = 2, minPatients = 20)
{
  statTable = as.data.frame(NULL)
  for(i in 1:length(survInfoList))
  {
    survInfo = survInfoList[[i]]
    #statTabRow = names(survInfoList)[i]
    dindInfo = metaGx::getDindexOfDatasets(survInfo, minPatients = minPatients)
    dindInfo$dIndex = as.numeric(as.character(dindInfo$dIndex))
    dindInfo$dIndStandErr = as.numeric(as.character(dindInfo$dIndStandErr))
    dindComb = combine.est(dindInfo$dIndex[!is.na(dindInfo$dIndex)], dindInfo$dIndStandErr[!is.na(dindInfo$dIndStandErr)], hetero = TRUE)
    dIndex = dindComb$estimate
    dIndexSe = dindComb$se
    dindexP = 2*pnorm(-abs(log(dindComb$estimate)/dindComb$se))
    
    dataGroups = groupPatientsViaScores(survInfo, numGroups, normalizeEsetScores)
    survObj = survfit(Surv(survInfo$timeToDeath, survInfo$vitalStat) ~ dataGroups)
    bb = survdiff(Surv(survInfo$timeToDeath, survInfo$vitalStat) ~ dataGroups,rho=0)
    logRankP = (1 - pchisq(bb$chisq, 1))
    
    seLow = sprintf(paste0("%.", numDigits, "f"), dIndex - 2*dIndexSe)
    seHigh = sprintf(paste0("%.", numDigits, "f"), dIndex + 2*dIndexSe)
    statTabVec = data.frame(sprintf(paste0("%.", numDigits, "f"), dIndex), paste("(", seLow, ", ", seHigh, ")", sep=""), format(dindexP, scientific =  TRUE, digits = numDigits), format(logRankP, scientific =  TRUE, digits = numDigits))
    colnames(statTabVec) = c("D Index", "D Index 95% CI", "D Index P", "Log Rank Test P")
    statTable = rbind(statTable, statTabVec)
    colnames(statTable) = c("D Index", "D Index 95% CI", "D Index P", "Log Rank Test P")
    statTable[, "Log Rank Test P"] = replace(as.character(statTable[, "Log Rank Test P"]), as.character(statTable[, "Log Rank Test P"]) == "0e+00", "1e-16")
    statTable[, "Log Rank Test P"] = replace(as.character(statTable[, "Log Rank Test P"]), as.character(statTable[, "Log Rank Test P"]) == "0.0e+00", "1e-16")
  }
  rownames(statTable) = names(survInfoList)
  return(statTable)
}

