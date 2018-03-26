#' Function to obtain survival info and risk predictions for an expression set (eSet object)
#'
#' This function returns a data frame with survival info and risk predictions for the patients in a given expression set object
#' @param eset an eset object
#' @param geneEntrezIds A character vector containing the entrez IDs for the genes in the signature being tested
#' @param geneDirecs A numeric vector composed of +1 and - 1 indicating the direction of association for each gene supplied in the geneEntrezIds vector.
#' +1 implies that the expression of that gene in a patient will be added to the patients score and -1 implies the expression will be subtracted from their score.
#' If one is looking for high scores to be associated with good survival, than genes that high expression is believed to lead to good prognosis should be given +1 and genes that high expression is believed
#' to lead to bad prognosis should be given a minus 1.
#' @param subtypeName Default is NULL, in which case all the patients in the eset are used. If a string corresponding to one of the subtypes in the esets eset$subtypes vector is supplied then only patients that are the subtype "subtypeName" will be used in the analysis.
#' @param survivalMetric a string specifying the type of survival analysis that will be investigated in the reports. Either "overall" for overall survival or "relapse" for relapse free survival
#' @param censorTime an integer specifying the point in time (years) at which the survival data must be censored. The default is 10 years
#' @param bestProbes a list of vectors containing the row indices of the probes in the eset@assayData$exprs expresson matrix that should be used. If bestprobes is NULL than getBestprobes will be used to select the probes that have the highest IQR for each probe ID
#' @return a list with 3 elements. "subtype", which has Subtypes identified by the subtyping classification model. "subtype.proba", which has the probabilities that a patient belongs to each subtype estimated by the subtyping classification model. "subtype.crisp", which has Crisp classes identified by the subtyping classification model.
#' @export
#' @examples
#'
#' geneEntrezIds = c("43847", "434768", "80070", "3620")
#' geneDirecs = c(1, 1, -1, -1)
#' esetsAndProbes = getEsetsProbesSubtypesEvents("ovarian", "overall", "verhaak", dataNames = c("GSE9891))
#' eset = esetsAndProbes$esets$GSE9891
#' bestProbes = esetsAndProbes$esetsBestProbes$GSE9891
#' survInfo = getSurvInfo(eset, geneEntrezIds, geneDirecs, survivalMetric = "overall", bestProbes = bestProbes)
#'
#' firstSubtype = as.character(esetsAndProbes$esets$E.MTAB.386$subtypes[1])
#' survInfoSubtype = getSurvInfo(eset, geneEntrezIds, geneDirecs, survivalMetric = "overall", subtypeName = firstSubtype)

getSurvInfo = function(eset, geneEntrezIds, geneDirecs, survivalMetric, subtypeName = NULL, censorTime = 10, bestProbes = NULL)
{
  survData = as.data.frame(NULL)
  dataInfo = eset@featureData@data;
  
  #below probably pointless as names are already the same
#  bestProbeNames = names(bestProbes)
#  entNames = as.character(eset@featureData@data$EntrezGene.ID[as.numeric(bestProbes)])
#  entNames = gsub(" ", "", entNames)
#  entNames = gsub(",", "///", entNames)
#  if(identical(entNames, names(bestProbes))){
#    dataInfoEntrezGene.ID = names(bestProbes)
#  }else{
    dataInfoEntrezGene.ID = as.character(dataInfo$EntrezGene.ID)
    dataInfoEntrezGene.ID = dataInfoEntrezGene.ID[as.numeric(bestProbes)]
#  }
    
  
  if(is.null(bestProbes))
    bestProbes = metaGx::getBestProbes(eset)
  
  survEventList = getSurvEventData(list(eset), survivalMetric)
  dataTimeToDeath = survEventList$eventTimeList[[1]]
  dataVitalStat = survEventList$eventList[[1]]
  if(!is.null(censorTime))
  {
    newDat = censor.time(dataTimeToDeath, dataVitalStat, time.cens = censorTime)
    dataTimeToDeath = newDat[[1]]
    dataVitalStat = newDat[[2]]
  }
  dataValsOrig = eset@assayData$exprs;
  
  subtypes = as.vector(eset$subtypes);
  subInds = which(subtypes == subtypeName)
  if(!is.null(subtypeName))
  {
    if(length(subInds) > 0){
      dataValsOrig = dataValsOrig[, subInds, drop = FALSE]
      dataTimeToDeath = dataTimeToDeath[subInds]
      dataVitalStat = dataVitalStat[subInds]
    }else{
      warning(paste("No subtype", subtypeName, "found in eset$subtype. Check the names of the eset$subtype elements and ensure the supplied subtypeName variable is present"))
      survData = NA
      return(survData)
    }
  }
  
  ##new start
  #dataInfoEntrezGene.ID = as.character(dataInfo$EntrezGene.ID)
  ##some missing ids that are "///"
  #slashStarts = which(substr(dataInfoEntrezGene.ID, 1, 1) == "/")
  #slashStLeng = nchar(dataInfoEntrezGene.ID[slashStarts])
  #slashRem = substr(dataInfoEntrezGene.ID[slashStarts], 4, slashStLeng)
  #dataInfoEntrezGene.ID[slashStarts] = slashRem
  #missingEnts = which(dataInfoEntrezGene.ID == "")
  ##set the missing ents to -100 so the probe selection works
  ##dont end up using any of them anyways
  #dataInfoEntrezGene.ID[missingEnts]  = -100
  
  ##multiple entrez ids have ent1///ent2///....///entN
  #multEnts = which(grepl("/", dataInfoEntrezGene.ID))
  #multStrings = dataInfoEntrezGene.ID[multEnts]
  #firstSlashes = rapply(gregexpr(pattern = "/", multStrings), function(x) head(x, 1))
  #firstEnts = substr(multStrings, 1, firstSlashes - 1)
  #dataInfoEntrezGene.ID[multEnts] = firstEnts
  
  #multEnts = which(grepl(",", dataInfoEntrezGene.ID))
  #multStrings = dataInfoEntrezGene.ID[multEnts]
  #firstSlashes = rapply(gregexpr(pattern = ",", multStrings), function(x) head(x, 1))
  #firstEnts = substr(multStrings, 1, firstSlashes - 1)
  #dataInfoEntrezGene.ID[multEnts] = firstEnts
  
  #dataInfoEntrezGene.ID = gsub(",", "/", dataInfoEntrezGene.ID)
  #multEnts = which(grepl("/", dataInfoEntrezGene.ID))
  #multStrings = dataInfoEntrezGene.ID[multEnts]
  #firstSlashes = rapply(gregexpr(pattern = "/", multStrings), function(x) head(x, 1))
  #firstEnts = substr(multStrings, 1, firstSlashes - 1)
  #dataInfoEntrezGene.ID[multEnts] = firstEnts
  #dataInfoEntrezGene.ID = as.numeric(dataInfoEntrezGene.ID)
  
  #old line
  #dataInfoEntrezGene.ID = as.numeric(dataInfo$EntrezGene.ID)
  
  dataVals = dataValsOrig[as.numeric(bestProbes), , drop = FALSE]
  dataInfoEntrezGene.ID = gsub(" ", "", dataInfoEntrezGene.ID) 
  dataInfoEntrezGene.ID = gsub(",", "///", dataInfoEntrezGene.ID)
  dataInfoEntrezGene.ID = paste0("/", dataInfoEntrezGene.ID, "/")
  #print("hello")
  dataInfoProbeset = as.vector(dataInfo$probeset[as.numeric(bestProbes)])
  dataInfoGene = as.vector(dataInfo$gene[as.numeric(bestProbes)])
  
  geneSigInds = c()
  geneSigEntrez = c()
  geneSigNames = c()
  dataGeneSigInds = c()
  genesFoundInd = c()
  for(i in 1:length(geneEntrezIds))
  {
    geneSigEntrez = paste0("/", geneEntrezIds[i], "/")
    #geneSigNames = c(geneSigNames, xOrig$geneNames[geneSigInd])
    #print(which(dataInfoEntrezGene.ID == geneSigEntrez[i]))
    if(length(which(grepl(geneSigEntrez, dataInfoEntrezGene.ID))) > 0)
      genesFoundInd = c(genesFoundInd, i)
    dataGeneSigInds = c(dataGeneSigInds, which(grepl(geneSigEntrez, dataInfoEntrezGene.ID)));
    #print(i)
    #print(which(dataInfoEntrezGene.ID == geneSigEntrez[i]))
  }
  #geneSigInds = which(xOrig$probeIds %in% geneSig);
  #geneSigEntrez = xOrig$entrezIds[geneSigInds]
  
  #tcgaGeneSigInds = which(tcgaInfo$EntrezGene.ID %in% geneSigEntrez);
  
  scores = NULL
  topGenes = dataInfoProbeset[dataGeneSigInds]
  
  #if none of the randomly selected genes are in the data set then move on
  if(length(topGenes) < 1)
  {
    survData = NA
    return(survData)
  }
  
  rownames(dataVals) = dataInfoProbeset
  for(i in 1:dim(dataVals)[2])
  {
    testSamp = t(dataVals[, i, drop = FALSE]);
    mysig <- cbind("probe"=topGenes, "EntrezGene.ID"=NA, "coefficient"= (as.numeric(geneDirecs[genesFoundInd])/sum(abs(as.numeric(geneDirecs[genesFoundInd])))))
    rownames(mysig) = topGenes;
    scores <- rbind(scores, cbind("score"=metaGx::calcSigScore(x=mysig, testSamp, annot=NULL, doMapping=FALSE, signed=TRUE)$score, "fold"=i))
  }
  #print("hi")
  scoreVals = as.numeric(scores[,1])
  names(scoreVals) = colnames(dataVals)
  invalidPatVec = c("start")
  
  unknown = which(is.na(dataTimeToDeath) == TRUE)
  if(length(unknown) > 0)
  {
    dataTimeToDeath = dataTimeToDeath[-unknown];
    dataVitalStat = dataVitalStat[-unknown]
    scoreVals = scoreVals[-unknown]
    invalidPatVec = c(invalidPatVec, unknown)
  }
  
  #Is this dirupting something else, why didnt I add this earlier if it was so easy --> double check
  unknown = which(is.na(dataVitalStat) == TRUE)
  if(length(unknown) > 0)
  {
    dataTimeToDeath = dataTimeToDeath[-unknown];
    dataVitalStat = dataVitalStat[-unknown]
    scoreVals = scoreVals[-unknown]
    invalidPatVec = c(invalidPatVec, unknown)
  }
  
  #some patient expression values are NA, causes NA scores
  missingGenes = which(is.na(scoreVals) == TRUE)
  if(length(missingGenes > 0))
  {
    dataTimeToDeath = dataTimeToDeath[-missingGenes];
    dataVitalStat = dataVitalStat[-missingGenes]
    scoreVals = scoreVals[-missingGenes]
    invalidPatVec = c(invalidPatVec, missingGenes)
  }
  
  scoreValsOrig = scoreVals
  quantVals = quantile(scoreVals, c(.025, .975))
  oldLow = quantVals[1]
  oldHigh = quantVals[2]
  newLow = -1
  newHigh = 1
  if(length(scoreValsOrig) > 1){
    scoreVals = newLow*(1 - (scoreVals - oldLow)/(oldHigh - oldLow)) + newHigh*((scoreVals - oldLow)/(oldHigh - oldLow))
  }else{
    #cant shift data if only 1 patient is present so do nothing
    scoreVals = scoreValsOrig
  }
  topGeneVec = rep(length(topGenes), length(dataTimeToDeath))
  
  #dataMatrix = matrix(c(dataTimeToDeath, dataVitalStat, scoreVals, scoreValsOrig, topGeneVec, names(scoreVals)), nrow = length(dataVitalStat), ncol = 6)
  #survData = as.data.frame(dataMatrix, stringsAsFactors = FALSE)
  survData = data.frame(dataTimeToDeath, dataVitalStat, scoreVals, scoreValsOrig, topGeneVec, names(scoreVals), stringsAsFactors = FALSE)
  colnames(survData) = c("timeToDeath", "vitalStat", "scoreVals", "scoreValsOrig", "numGenesPresent", "Patient ID")
  
  return(survData)
}
