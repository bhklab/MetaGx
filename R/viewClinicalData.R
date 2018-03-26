
#' Function to view the available clinical data available for each dataset
#'
#' This function returns the data required to generate a survival report according to the user inputs that define the cancer type and survival type being considered
#' @param cancerType a string specifying what cancer type to use in the analysis. The options are "breast" and "ovarian"
#' @param plots a boolean specifying whether to graphically illustrate what data is available in the datasets. Default is TRUE
#' @return a matrix containing the perentage of patients with the given clinical info in each dataset. A graphical illustration
#' of the info will also be generated. This can be turned off by setting plots to FALSE
#' @export
#' @examples
#'
#' dataClinicalInfo = viewClinicalData = function("ovarian")
#' 

viewClinicalData = function(cancerType, plots = TRUE)
{
  #pData Variables
  if(cancerType == "breast")
  {
    pDataID <- c("er","pgr", "her2", "age_at_initial_pathologic_diagnosis", 
                 "grade", "dmfs_days", "dmfs_status", "days_to_tumor_recurrence", 
                 "recurrence_status", "days_to_death", "vital_status", "sample_type", "treatment")
    
    dataListAndDups = loadBreastEsets()
    dataListTcgaMetaAndDups = loadBreastEsets(loadString = c("TCGA", "METABRIC"))
    esets = c(dataListAndDups$esets, dataListTcgaMetaAndDups$esets)
  }else if(cancerType == "ovarian"){
    pDataID <- c("sample_type", "histological_type", "primarysite", "summarygrade",
                 "summarystage", "tumorstage", "grade",
                 "age_at_initial_pathologic_diagnosis", "pltx", "tax",
                 "neo", "days_to_tumor_recurrence", "recurrence_status",
                 "days_to_death", "vital_status")
    dataListAndDups = loadOvarianEsets()
    esets = dataListAndDups$esets
  }else{
    stop("invalid cancerType string, only ovarian and breast are currenly available")
  }

  numSamples <- vapply(seq_along(esets), FUN=function(i, esets){
    length(sampleNames(esets[[i]]))
  }, numeric(1), esets=esets)
  
  
  SampleNumberSummaryAll <- data.frame(NumberOfSamples = numSamples,
                                       row.names = names(esets))
  total <- sum(SampleNumberSummaryAll[,"NumberOfSamples"])
  SampleNumberSummaryAll <- rbind(SampleNumberSummaryAll, total)
  rownames(SampleNumberSummaryAll)[nrow(SampleNumberSummaryAll)] <- "Total"
  
  pDataPercentSummaryTable <- NULL
  pDataSummaryNumbersTable <- NULL
  
  pDataSummaryNumbersList = lapply(esets, function(x)
    vapply(pDataID, function(y) sum(!is.na(pData(x)[,y])), numeric(1)))
  
  pDataPercentSummaryList = lapply(esets, function(x)
    vapply(pDataID, function(y)
      sum(!is.na(pData(x)[,y]))/nrow(pData(x)), numeric(1))*100)
  
  pDataSummaryNumbersTable = sapply(pDataSummaryNumbersList, function(x) x)
  pDataPercentSummaryTable = sapply(pDataPercentSummaryList, function(x) x)
  
  rownames(pDataSummaryNumbersTable) <- pDataID
  rownames(pDataPercentSummaryTable) <- pDataID
  colnames(pDataSummaryNumbersTable) <- names(esets)
  colnames(pDataPercentSummaryTable) <- names(esets)
  
  pDataSummaryNumbersTable <- rbind(pDataSummaryNumbersTable, total)
  rownames(pDataSummaryNumbersTable)[nrow(pDataSummaryNumbersTable)] <- "Total"
  
  
  # Generate a heatmap representation of the pData
  pDataPercentSummaryTable<-t(pDataPercentSummaryTable)
  pDataPercentSummaryTable<-cbind(Name=(rownames(pDataPercentSummaryTable))
                                  ,pDataPercentSummaryTable)
  
  nba<-pDataPercentSummaryTable
  gradient_colors = c("#ffffff","#ffffd9","#edf8b1","#c7e9b4","#7fcdbb",
                      "#41b6c4","#1d91c0","#225ea8","#253494","#081d58")
  
  library(lattice)
  nbamat<-as.matrix(nba)
  rownames(nbamat)<-nbamat[,1]
  nbamat<-nbamat[,-1]
  Interval<-as.numeric(c(10,20,30,40,50,60,70,80,90,100))
  
  print(levelplot(nbamat,col.regions=gradient_colors,
            main="Available Clinical Annotation",
            scales=list(x=list(rot=90, cex=0.5),
                        y= list(cex=0.5),key=list(cex=0.2)),
            at=seq(from=0,to=100,length=10),
            cex=0.2, ylab="", xlab="", lattice.options=list(),
            colorkey=list(at=as.numeric(factor(c(seq(from=0, to=100, by=10)))),
                          labels=as.character(c( "0","10%","20%","30%", "40%","50%",
                                                 "60%", "70%", "80%","90%", "100%"),
                                              cex=0.2,font=1,col="brown",height=1,
                                              width=1.4), col=(gradient_colors))))
  
  return(nbamat)
}



