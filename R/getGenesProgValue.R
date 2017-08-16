
#' Function that assesses the prognostic value of a gene signature and the genes in the signature
#'
#' This function runs a survival analysis on breast or ovarian cancer patients for a gene signature and the genes in the signature using overall survival or relapse freee survival data
#' @param geneSigList A list of character vectors containing the ensemble IDs, entrez IDs, or gene symbols for the genes signatures/genes to conduct a survival analysis on. Entrez IDs are 
#'  recommended as if gene symbols or ensemble Ids cant be mapped to entrez Ids then the gene will be omitted from the analysis. Note that The names of the
#' elements in the list will correspond to the names of the signatures in the report, and default names will be provided if the list elements are not named.
#' @param geneDirecList A list of numeric vectors composed of +1 and - 1 indicating the direction of association for each vector of genes supplied in the geneSigList vectors.
#' +1 implies that the expression of that gene in a patient will be added to the patients score and -1 implies the expression will be subtracted from their score.
#' If one is looking for high scores to be associated with good survival, than genes that high expression is believed to lead to good prognosis should be given +1 and genes that high expression is believed
#' to lead to bad prognosis should be given a minus 1.
#' @param cancerType A string representing the type of cancer one would like to test the gene signature on. Options are currently "ovarian" and "breast"
#' @param subtype a string specifying what type of survival data to use in the analysis. If the cancerType is set to "ovarian"
#' then enter "overall" for overall survival and "relapse" for relapse free survival. If the cancerType is set to "breast" then the options
#' are "relapse" for relapse free survival, "overallMeta" for overall survival on the metabric study data, "overallTcga" for
#' overall survival on the TCGA study data, and "overall" for overall survival on the remaining studies (No Metabric and TCGA data).
#' @param survivalMetric a string specifying what type of survival data to use in the analysis. If the cancerType is set to "ovarian"
#' then enter "overall" for overall survival and "relapse" for relapse free survival. If the cancerType is set to "breast" then the options
#' are "relapse" for relapse free survival, "overallMeta" for overall survival on the metabric study data, "overallTcga" for
#' overall survival on the TCGA study data, "overall" for overall survival on the remaining studies (No Metabric and TCGA data), and metastasis for distant metastasis free survival.
#' Note that for all cancer types the option "hierarchy" is available and this survivalMetric parameter checks for relapse free, then distant metastasis, and then overall survival data
#' in each dataset and keeps the first survival data it finds. Note that the ovarian datasets do not have distant metastasis free survival info and for the breast analyses metabric and tcga are removed.
#' @param numGroups an integer specifying the number of groups the patients will be split into when generating survival curves and calculating the log rank p value. The value is 2 by default and cannot exceed 5.
#' @param dataNames a character vector specifying the names of the datasets to include in the analysis. The names must match the names in the column "Data Name" from the dataframe returned from the
#' obtainDataInfo function. Use loadMetaData followed by obtainDataInfo with the cancerType and survivalMetric of interest to get the appropriate data names in the obtainDataInfo table for your analysis. By default
#' all datasets are included.
#' @param soloGeneAnalysis a boolean specifying whether each unique gene in geneSigList should have an survival analysis conducted on it in order to assess its
#' prognostic value independent of its gene signature. Results will show up in a section called Individual Gene Survival Analysis and high scores in the survival curves correspond to high expression as the direction is defaulted to 1
#' to allow for easy comparison amongst all the genes. Default value is FALSE.
#' @param removeMid a number greater than 0 and less than .5 specifying the fraction of the patients with scores/risk predictions in the middle of the patients scores to be removed when generating survival results
#' The default is 0, all patients used in the analysis
#' @param censorTime a number specifying the point in time (years) at which the survival data must be censored. The default is 10 years
#' @param addBenchmarks a boolean specifying whether to add known signatures from literature to the analysis for comparison to the provided signatures. Default value is FALSE
#' @return A list containing the results used to generate the pdf report made by the createSurvivalReport function. The list has the elements patientSurvData, genesPrognosticVal, genesDindexInfo, genesPrognosticSummary, datasets, genesStatus, and genesInfo.
#' patientSurvData has the survival data frames for the signature and each gene (1 frame for each subtype and all patients), genesPrognosticValue has the D index and log p values for the signature and individual genes on the patients and subtypes.
#' genesDindexInfo has the D indices and their standard errors for the signature and the genes in each individual dataset (genesPrognosticValue D indices are the meta estimates from using the D index in each dataset), genesPrognosticSummary is the info in
#' genesPrognosticVal summarized in tables where there is one table for each subtype and the analysis on all patients, geneStatusVec is a vector of strings that indicate whether a gene was in the external datasets and used in the analysis, genesInfo is
#' a list with dataframes with info pertaining to each signature, and sigCorrelationMatrix is a matrix containing the correlations between the signatures (NULL if not multiple signatures).
#' @export
#' @examples
#'
#' geneSigList = list(c("KLK14","RHOX8","ADAMTS20","IDO1"), c("ADAMTS20","IDO1", "FAP"))
#' names(geneSigList) = c("4 Gene Signature", "3 Gene Signature")
#' geneDirecList = list(c(1, 1, -1, -1), c(-1, -1, 1))
#' sigInfoList = getGenesProgValue(geneSigList, geneDirecList, cancerType = "ovarian", subtype = "verhaak", survivalMetric = "overall")
#'

getGenesProgValue = function(geneSigList, geneDirecList, cancerType, subtype, survivalMetric, numGroups = 2, dataNames = NULL, soloGeneAnalysis = FALSE, removeMid = 0, censorTime = 10, addBenchmarks = FALSE)
{
  geneSigOrigList = geneSigList
  #geneInfoList = geneSigList
  if(addBenchmarks == TRUE)
  {
    cancerSigs = obtainCancerSigs(cancerType)
    names(cancerSigs$geneSigList) = paste(names(cancerSigs$geneSigList), "Benchmark")
    geneSigList = c(geneSigList, cancerSigs$geneSigList)
    geneDirecList = c(geneDirecList, cancerSigs$geneDirecList)
  }
  geneInfoList = geneSigList
  #bad assumption that1 gene sent to the routine will be mapped, will crash if the assumption isnt true below at geneFrame = getGeneInfo
  mapVec = rep(TRUE, length(geneSigList))
  
  geneEntrezList = list()
  geneNameList = list()

  for(i in 1:length(geneSigList))
  {
    geneFrame = metaGx::getGeneInfo(geneSigList[[i]], geneDirecList[[i]])
    notMapped = which(is.na(geneFrame), arr.ind=TRUE)[, "row"]
    #below wont work/is useless as getGeneInfo wouldve had an error if none of the entrez ids could be mapped
    if(length(notMapped) == nrow(geneFrame))
      mapVec[i] = FALSE
    geneEntrezList[[i]] = geneFrame$`Entrez ID`
    if(is.null(names(geneSigList))){
      geneNameList[[i]] = paste("signature", i)
    }else if(is.na(names(geneSigList)[i])){
      geneNameList[[i]] = paste("signature", i)
    }else{
      geneNameList[[i]] = names(geneSigList)[i]
    }
  }

  if(soloGeneAnalysis == TRUE)
  {
    uniqueIds = unique(unlist(geneSigOrigList))
    geneSoloFrame = getGeneInfo(uniqueIds, rep(1, length(uniqueIds)))
    notMapped = which(is.na(geneSoloFrame), arr.ind=TRUE)[, "row"]
    for(i in 1:nrow(geneSoloFrame))
    {
      geneEntrezList[[length(geneEntrezList) + 1]] = geneSoloFrame$`Entrez ID`[i]
      geneDirecList[[length(geneEntrezList)]] = 1
      if(!is.na(geneSoloFrame$Symbol[i]))
        geneNameList[[length(geneEntrezList)]] = geneSoloFrame$Symbol[i]
      else
        geneNameList[[length(geneEntrezList)]] = uniqueIds[i]
      geneInfoList[[length(geneEntrezList)]] = notMapped[i]
      mapVec = c(mapVec, !is.na(geneSoloFrame$Description[i]))
    }
  }

  geneStatusVec = c(rep("In datasets", length(geneEntrezList)))
  notMapped = which(is.na(geneEntrezList))
  geneStatusVec[notMapped] = "Could not map to entrez gene ID"
  if(length(notMapped) == 0)
    notMapped = 0

  esetsAndProbes = getEsetsProbesSubtypesEvents(cancerType, survivalMetric, subtype, dataNames)
  dataList = esetsAndProbes$esets
  dataListBestProbes = esetsAndProbes$esetsBestProbes
  infoString = names(dataList)

  geneDataList = list()
  forestPlotDataList = list()
  geneTabList = list()
  geneFrameList = list()
  for(i in 1:length(geneEntrezList))
  {
    print(i)
    if(mapVec[i] == FALSE){
      geneFrameList[[i]] = NA
    }else{
      if(i <= length(geneSigList))
        geneFrameList[[i]] = getGeneInfo(geneSigList[[i]], geneDirecList[[i]])
      else
        geneFrameList[[i]] = getGeneInfo(geneEntrezList[[i]], geneDirecList[[i]])
    }
    geneDataList[[i]] = getPatientSurvivalData(as.numeric(geneEntrezList[[i]]), geneDirecList[[i]], cancerType, subtype, survivalMetric, dataList = dataList, bestProbesList = dataListBestProbes, removeMid = removeMid, censorTime = censorTime)
    names(geneDataList)[i] = geneNameList[[i]]
    names(geneFrameList)[i] = geneNameList[[i]]

    if(nrow(geneDataList[[i]]$`All Patients`) > 0){
      geneTabList[[i]] = getMetaDindexAndLogP(geneDataList[[i]], as.numeric(numGroups))
      forestPlotDataList[[i]] = lapply(geneDataList[[i]], function(x) getDindexOfDatasets(x))
    }else{
      geneTabList[[i]] = NA
      forestPlotDataList[[i]] = NA
      if(nrow(geneDataList[[i]]$`All Patients`) == 0)
        geneStatusVec[i] = "Gene not present in the datasets used"
    }
    names(geneTabList)[i] = geneNameList[[i]]
    names(forestPlotDataList)[i] = geneNameList[[i]]
    names(geneStatusVec)[i] = geneNameList[[i]]
  }

  genesSubtypeTabList = list()
  subtypeNames = unique(unlist(lapply(geneTabList, function(x) rownames(x))))
  for(i in 1:length(subtypeNames))
  {
    genesSubtypeTab = as.data.frame(NULL)
    subName = subtypeNames[i]
    for(j in 1:length(geneTabList))
    {
      rowInd = which(rownames(geneTabList[[j]]) == subName)
      geneTab = geneTabList[[j]]
      if(length(rowInd) > 0)
      {
        genesSubtypeTab = rbind(genesSubtypeTab, geneTabList[[j]][rowInd, ])
        rownames(genesSubtypeTab)[nrow(genesSubtypeTab)] = names(geneTabList)[j]
      }
    }
    genesSubtypeTabList[[i]] = genesSubtypeTab
  }
  names(genesSubtypeTabList) = subtypeNames

  setupTable = function(table)
  {
    geneRowNames = rownames(table)
    table = cbind(table, geneRowNames)
    table = unique(table)
    table = table[,1:4]

    table = cbind(table, format(p.adjust(as.vector(table[,3]), method = "fdr"), scientific =  TRUE, digits = 2))
    table = cbind(table, format(p.adjust(as.vector(table[,4]), method = "fdr"), scientific =  TRUE, digits = 2))

    colnames(table) = c("D Index", "D Index 95% CI", "D Index P", "Log Rank Test P", "D Index FDR", "Log Rank FDR")
    #colnames(table) = c("D Index P Val", " D Index", "D Index Standard Error", "Log Rank Test P Val")
    rownames(table) = unique(geneRowNames)
    sortInd = order(as.numeric(as.character(table[,"Log Rank Test P"])), decreasing = FALSE)
    table = table[sortInd, ]
    return(table)
  }
  genesSubtypeTabList = lapply(genesSubtypeTabList, function(x) setupTable(x))
  
  corMatrix = list()
  
  reportList = list(geneDataList, geneTabList, forestPlotDataList, genesSubtypeTabList, dataList, geneStatusVec, geneFrameList, corMatrix)
  names(reportList) = c("patientSurvData", "genesPrognosticVal", "genesDindexInfo", "genesPrognosticSummary", "datasets", "genesStatus", "genesInfo", "sigCorrelationList")
  
#  geneEntrezList = lapply(1:length(geneSigList), function(x) as.numeric(geneEntrezList[[x]]))
#  names(geneEntrezList) = names(geneSigList)[1:length(geneEntrezList)]
  #should add patient Ids to survInfo frames and use result of getPatientSurviuvalData (scoreVals) for corMatrix
  #as its slower than entire rest of function

  subtypeNames = tolower(names(reportList$patientSurvData[[1]]))
  if(length(geneSigList) > 1)
  {
    numSigs = length(geneSigOrigList)
    if(addBenchmarks == TRUE)
      numSigs = numSigs + length(obtainCancerSigs(cancerType)$geneSigList)
    reportListSig = reportList
    if(length(reportListSig$patientSurvData) > numSigs)
      reportListSig$patientSurvData[(numSigs + 1):length(reportList$patientSurvData)] = NULL
    for(i in 1:length(subtypeNames))
    {
      corMatrix[[i]] = correlateGeneSigs(reportListSig, subtypeName = subtypeNames[i])  
    }
    names(corMatrix) = names(reportList$patientSurvData[[1]])
  }
  reportList$sigCorrelationList = corMatrix
  #want to return survivalFrames, geneTables, geneSubtypeSummaryTables, ForestPlots (sig and genes), dataList, and geneStatus (found, not mapped, or not in datasets)
  #generate remaining needed things in Rnw file
  #have geneSig stuff just be another element in the lists, aka get rid of analyze signature routine
  #feed setup report into this function
  #geneInfoLists = list(geneInfoList, list(geneSymbVec, pathwayNameVec, geneNameList, geneDescripList, geneNameList), geneFileInfo, statTableList, summaryTabList, indexInfo, geneSigInfoList, geneSigFrame, forestPlotFrame, dataTabInfo, statTableSigList)

  return(reportList)
}


#geneSigList = list(c("KLK14","RHOX8","ADAMTS20","IDO1"), c("ADAMTS20","IDO1", "FAP"))
#names(geneSigList) = c("4 Gene Signature", "3 Gene Signature")
#geneDirecList = list(c(1, 1, -1, -1), c(-1, -1, 1))
#survInfoList = getGenesProgValue(geneSigList, geneDirecList, cancerType = "breast", subtype = "intclust", survivalMetric = "overall", addBenchmarks = TRUE, soloGeneAnalysis = TRUE)
#cancerType = "breast"
#subtype = "intClust"
#survivalMetric = "overall"
#numGroups = 2
#dataNames = NULL
#soloGeneAnalysis = TRUE
#removeMid = 0
#censorTime = 10
#addBenchmarks = TRUE

