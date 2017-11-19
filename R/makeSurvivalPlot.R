
#' Function to generate a Kaplan-Meier survival curve using risk predictons and survival data
#'
#' This function generates a Kaplan-Meier survival curve using the risk predictons and survival data stored in the survInfo data frame
#' @param survInfo One of the data frames in the list returned from the function getPatientSurvivalData. Alternatively, a data frame with 6 columns called eventTime, eventStatus, scoreVals, scoreVals, numGenesPresent, dataName and patient ID that correspond to
#' the survival event times, the survival event's status, the score/risk prediction of the patient, the number of genes used to calculate the score, the name of the dataset that the patient is from and the ID of the patient int that dataset. 
#' @param numGroups an integer specifying the number of groups the patients will be split into when generating survival curves and calculating the log rank p value. The value is 2 by default and cannot exceed 5.
#' @param normalizeEsetScores A boolean specifying whether the groups used to determine the log rank p value should be determined by 1. Splitting the scores/risk predictions of the patients from each dataset in the number of groups specified or
#' 2. normalizing the scores of the patients in each dataset such that the 2.5 percentile and 97.5 percentile for the patients scores in each dataset have the same value, then creating one vector of patient score and splitting this
#' vector of patient score sinto the number of groups specified (TRUE). Normalization/TRUE is the default. Note that this generally improves results when usign 1 type of survival data (better comparison of patients across datasets), but is inappropriate/can ruin survival curve splits
#' when using multiple survival events (e.g overlall + relapse) likely because the scores do not relate to the survival event the same way so the split needs to hapen individually in each dataset where only 1 type of survival event is present.
#' @param titleStr a string to be used as the title of the plot. If not supplied a default title will be used.
#' @param xLoc the location of the D index and lo rank p value along the x-axis in the plot. By default 2
#' @param yLabel a string specifying the y axis title in the survival plot. Default is NULL resulting in "probability of survival" in the plot
#' @return a Kaplan-Meier survival curve with the corresponding D index and log rank p value in the plot
#' @export
#' @examples
#' 
#' geneEntrezIds = c("43847", "434768", "80070", "3620")
#' geneDirecs = c(1, 1, -1, -1)
#' survInfoList = getPatientSurvivalData(geneEntrezIds, geneDirecs, cancerType = "ovarian", subtype = "verhaak", survivalMetric = "overall")
#' makeSurvivalPlot(survInfoList$`All Patients`, numGroups = 2)

makeSurvivalPlot = function(survInfo, numGroups,  normalizeEsetScores = FALSE, titleStr = "Kaplan-Meier Survival Curve", xLoc = 2., yLabel = NULL)
{
  #makeSurvPlot(geneSigInfoList[[geneSigInfInd]], geneSigInfoList[[geneSigInfInd + 1]], geneSigInfoList[[geneSigInfInd + 2]], params$threeGroups, dindStr = dindStrSig, titleStr)
  
  if(is.null(yLabel))
    yLabel = "Probability of Survival"
  
  metaDindAndLogRankP = getMetaDindexAndLogP(list(survInfo), numGroups, normalizeEsetScores)
  dindStr = as.character(metaDindAndLogRankP$`D Index`)
  logRankP = as.character(metaDindAndLogRankP$`Log Rank Test P`)
  if(logRankP == 0)
    logRankP = 1.00*10^-16
  
  dataGroups = groupPatientsViaScores(survInfo, numGroups, normalizeEsetScores)
  datFrame = data.frame("surv.time" = survInfo$timeToDeath, "surv.event"= survInfo$vitalStat, "strat"= dataGroups)
  ddweights <- array(1, dim=nrow(datFrame))
  #km.coxph.plot(formula.s = Surv(survInfo$timesToDeath, survInfo$vitalStats) ~ survInfo$groups)
  if(numGroups == 2){
    colStr = c("darkblue", "darkred")
    legTextStr = paste(c("High Score", "Low Score"), "   ", sep="")
    #legTextStr = paste(c("High Score","", "Low Score"), "   ", sep="")
  }else if(numGroups == 3){
    colStr = c("darkblue", "darkgreen", "darkred")
    legTextStr = c("High Score", "Intermediate Score", "Low Score")
  }else if(numGroups == 4){
    colStr = c(1:4)
    legTextStr = c("High Score", "2nd Highest", "2nd Lowest", "Low Score")
  }else{
    colStr = c(1:numGroups)
    legTextStr = c("Highest Score")
    for(i in 2:(numGroups - 1))
      legTextStr = c(legTextStr, paste("High Score", i))
    legTextStr = c(legTextStr, "Lowest Score")
  }
  
  km.coxph.plot(formula.s=Surv(surv.time, surv.event) ~ strat, data.s=datFrame,
                weight.s=ddweights, x.label="Time (years)", y.label=yLabel,
                main.title=titleStr, leg.text=legTextStr,
                leg.pos="topright", leg.inset=0.0, .col=colStr,
                .lty=c(1,1), show.n.risk=TRUE, n.risk.step=2, n.risk.cex=0.85, verbose=FALSE, leg.bty = "n", o.text = "")
  
  #legend(y.intersp = 2)
  text(xLoc + 0.4, 0.115, paste("D index = ", dindStr))
  text(xLoc + 0.5, 0.03, paste("Log Rank P = ", format(logRankP, scientific =  TRUE, digits = 2)))

}
