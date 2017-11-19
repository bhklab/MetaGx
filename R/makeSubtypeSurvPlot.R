
#' Function to generate a Kaplan-Meier survival curve using subtype groupings and survival data
#'
#' This function generates a Kaplan-Meier survival curve using the subtypes and survival data stored in the survInfo data frames
#' @param survInfoList The data frames in the list returned from the function getPatientSurvivalData. Alternatively, a list of data frames with 6 columns called eventTime, eventStatus, scoreVals, scoreVals, numGenesPresent, dataName and patient ID that correspond to
#' the survival event times, the survival event's status, the score/risk prediction of the patient, the number of genes used to calculate the score, the name of the dataset that the patient is from and the ID of the patient int that dataset. 
#' @param titleStr a string to be used as the title of the plot. If not supplied a default title will be used.
#' @param xLoc the location of the D index and lo rank p value along the x-axis in the plot. By default 2
#' @param yLabel a string specifying the y axis title in the survival plot. Default is NULL resulting in "probability of survival" in the plot
#' @return a Kaplan-Meier survival curve with the patient groups dictated by the subtypes and a log rank p value in the plot
#' @export
#' @examples
#' 
#' geneEntrezIds = c("43847", "434768", "80070", "3620")
#' geneDirecs = c(1, 1, -1, -1)
#' survInfoList = getPatientSurvivalData(geneEntrezIds, geneDirecs, cancerType = "ovarian", subtype = "verhaak", survivalMetric = "overall")
#' makeSurvivalPlot(survInfoList)

makeSubtypeSurvPlot = function(survInfoList, titleStr = "Kaplan-Meier Survival Curve", xLoc = 2., yLabel = NULL)
{
  #makeSurvPlot(geneSigInfoList[[geneSigInfInd]], geneSigInfoList[[geneSigInfInd + 1]], geneSigInfoList[[geneSigInfInd + 2]], params$threeGroups, dindStr = dindStrSig, titleStr)
  
  if(is.null(yLabel))
    yLabel = "Probability of Survival"
  
  strat = c()
  patNameVec = rownames(survInfoList$`All Patients`)
  survInfo = as.data.frame(NULL)
  for(i in 1:length(patNameVec))
  {
    for(j in 2:length(survInfoList))
    {
      rowInd = which(rownames(survInfoList[[j]]) == patNameVec[i])
      if(length(rowInd) > 0)
      {
        strat = c(strat, names(survInfoList)[j])
        survInfo = rbind(survInfo, survInfoList[[j]][rowInd, ])
      }
    }
  }
  
  survObj = survfit(Surv(survInfo$timeToDeath, survInfo$vitalStat) ~ strat)
  bb = survdiff(Surv(survInfo$timeToDeath, survInfo$vitalStat) ~ strat, rho = 0)
  logRankP = (1 - pchisq(bb$chisq, 1))
  if(logRankP == 0)
    logRankP = 1.00*10^-16
  logRankP = as.character(format(logRankP, scientific = TRUE, digits = 2))

  datFrame = data.frame("surv.time" = survInfo$timeToDeath, "surv.event"= survInfo$vitalStat, "strat"= strat)
  ddweights <- array(1, dim=nrow(datFrame))
  
  data.s = datFrame
  weight.s = array(1, dim = nrow(data.s), dimnames = list(rownames(data.s)))
  formula.s=Surv(surv.time, surv.event) ~ strat
  data.s <- data.s[!is.na(weight.s) & weight.s > 0, , drop = FALSE]
  weight.s <- weight.s[!is.na(weight.s) & weight.s > 0]
  assign("weight.s", weight.s, envir = .GlobalEnv)
  coxInfo = survival::coxph(formula.s, data = data.s, weights = weight.s)
  survfitObj = survfit(formula.s, data = data.s, weights = weight.s)
  infoTab = summary(survfitObj)$table
  gsub("strat=", "", rownames(infoTab))
  
  legTextStr = gsub("strat=", "", rownames(infoTab))

  
  km.coxph.plot(formula.s=Surv(surv.time, surv.event) ~ strat, data.s=datFrame,
                weight.s=ddweights, x.label="Time (years)", y.label=yLabel,
                main.title=titleStr, leg.text=legTextStr,
                leg.pos="topright", leg.inset=0,
                .lty=c(1,1), show.n.risk=TRUE, n.risk.step=2, n.risk.cex=0.75, verbose=FALSE, leg.bty = "n", o.text = "")
  
  text(xLoc + 0.5, 0.03, paste("Log Rank P = ", format(logRankP, scientific =  TRUE, digits = 2)))
  #legend(y.intersp = 2)

}



