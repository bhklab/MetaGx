#' Function to load data required for the survival reports
#'
#' This function returns the data required to generate a survival report according to the user inputs that define the cancer type and survival type being considered
#' @param cancerType a string specifying what cancer type to use in the analysis. The options are "breast" and "ovarian"
#' @param survivalMetric a string specifying what type of survival data to use in the analysis. If the cancerType is set to "ovarian"
#' then enter "overall" for overall survival and "relapse" for relapse free survival. If the cancerType is set to "breast" then the options
#' are "relapse" for relapse free survival, "overallMeta" for overall survival on the metabric study data, "overallTcga" for
#' overall survival on the TCGA study data, "overall" for overall survival on the remaining studies (No Metabric and TCGA data), and metastasis for distant metastasis free survival.
#' Note that for all cancer types the option "hierarchy" is available and this survivalMetric parameter checks for relapse free, then distant metastasis, and then overall survival data
#' in each dataset and keeps the first survival data it finds. Note that the ovarian datasets do not have distant metastasis free survival info and for the breast analyses metabric and tcga are removed.
#' @param dataNames a character vector specifying the names of the datasets to include in the analysis. The names must match the names in the column "Data Name" from the dataframe returned from the
#' obtainDataInfo function. Use loadMetaData followed by obtainDataInfo with the cancerType and survivalMetric of interest to get the appropriate data names in the obtainDataInfo table for your analysis. By default
#' all datasets are included.
#' @param includeAll a boolean specifying whether to include the TCGA and METABRIC datasets in the breast cancer analysis. Default is FALSE
#' @return a list containing the datasets that will be used in the analysis
#' @export
#' @examples
#'
#' dataList = loadMetaData(cancerType = "breast", survivalMetric = "hierarchy")

loadMetaData = function(cancerType, survivalMetric, dataNames = NULL, includeAll = FALSE)
{

  #will need to
  #1. change so breast data is loaded via metaGxBreast
  #2. consider adapting so that appropriate data (overall/relapse free) is found automatically
  #likely uneccessary/not good as some data does not cooperate with subtyping methods
  #would have to adapt subtyping methods to remove data and return new data?
  dataList = list()
  package.name = NULL

  #setwd("D:\\PMH Data\\Breast Cancer Data\\Breast Cancer Data Haibe-Kains\\wetransfer-474706")
  #install.packages("MetaGxBreast_0.99.0.tar.gz", type ="source", repos = NULL)
  if(cancerType == "ovarian")
  {
    #library(survcomp)
    #library("MetaGxOvarian")
    #is below line out okay?

    source(system.file("extdata", "patientselection.config", package="MetaGxOvarian"))
    source(system.file("extdata", "createEsetList.R", package="MetaGxOvarian"))
    
    dataList[[1]] = E.MTAB.386
    dataList[[2]] = GSE2109
    dataList[[3]] = GSE6008
    dataList[[4]] = GSE6822
    dataList[[5]] = GSE8842
    dataList[[6]] = GSE9891
    dataList[[7]] = GSE12418
    dataList[[8]] = GSE12470
    dataList[[9]] = GSE13876
    dataList[[10]] = GSE14764
    dataList[[11]] = GSE17260
    dataList[[12]] = GSE18520
    #below dataset every row has NA values
    #dataList[[13]] = GSE19829
    dataList[[13]] = GSE20565
    dataList[[14]] = GSE26193
    dataList[[15]] = GSE26712
    dataList[[16]] = GSE30009
    dataList[[17]] = GSE30161
    dataList[[18]] = GSE32062
    dataList[[19]] = GSE32063
    dataList[[20]] = GSE44104
    dataList[[21]] = GSE49997
    dataList[[22]] = GSE51088
    dataList[[23]] = PMID15897565
    dataList[[24]] = PMID17290060
    dataList[[25]] = PMID19318476
    dataList[[26]] = TCGAOVARIAN
    #tcga.rnaseq is a subset of the tcgaovarian data
    #dataList[[27]] = TCGA.RNASeqV2

    #PMID158 and PMID 193 are not mentioned in Deena's metaGx paper?

    infoString = c("E.MTAB.386", "GSE2109" ,"GSE6008", "GSE6822", "GSE8842", "GSE9891", "GSE12418",
                   "GSE12470", "GSE13876", "GSE14764", "GSE17260", "GSE18520",
                    "GSE20565", "GSE26193", "GSE26712", "GSE30009", "GSE30161", "GSE32062"
                   , "GSE32063", "GSE44104", "GSE49997", "GSE51088", "PMID15897565", "PMID17290060", "PMID19318476", "TCGAOVARIAN")

  }else if(cancerType == "breast"){
    #if(mord == FALSE)
    #setwd("C:\\Users\\micha\\Documents\\MetaGxBreast_0.99.0\\MetaGxBreast\\data")
    #if(mord == TRUE)
    #  setwd("//mnt//work1//users//home2//mzon//sweaveBreastAll//data")

    #if(mord = TRUE)
    
    source(system.file("extdata", "patientselection.config", package="MetaGxBreast"))
    source(system.file("extdata", "createEsetList.R", package="MetaGxBreast"))

    infoString = c("CAL", "DFHCC", "DFHCC2", "DFHCC3", "DUKE", "DUKE2", "EMC2", "EORTC10994", "EXPO",
                   "FNCLCC", "GSE25066", "GSE32646","GSE48091","GSE58644","HLP","IRB", "KOO","LUND",
                   "LUND2","MAINZ","MAQC2","MCCC","MDA4","METABRIC","MSK","MUG","NCCS","NCI","NKI",
                   "PNC","STK","STNO2","TCGA","TRANSBIG","UCSF","UNC4","UNT","UPP","VDX")

    #for(i in 1:length(infoString))
    #  data(list = infoString[i])

    dataList[[1]] = CAL
    dataList[[2]] = DFHCC
    dataList[[3]] = DFHCC2
    dataList[[4]] = DFHCC3
    dataList[[5]] = DUKE
    dataList[[6]] = DUKE2
    dataList[[7]] = EMC2
    dataList[[8]] = EORTC10994
    dataList[[9]] = EXPO
    dataList[[10]] = FNCLCC
    dataList[[11]] = GSE25066
    dataList[[12]] = GSE32646
    dataList[[13]] = GSE48091
    dataList[[14]] = GSE58644
    dataList[[15]] = HLP
    dataList[[16]] = IRB
    dataList[[17]] = KOO
    dataList[[18]] = LUND
    dataList[[19]] = LUND2
    dataList[[20]] = MAINZ
    dataList[[21]] = MAQC2
    dataList[[22]] = MCCC
    dataList[[23]] = MDA4
    dataList[[24]] = METABRIC
    dataList[[25]] = MSK
    dataList[[26]] = MUG
    dataList[[27]] = NCCS
    dataList[[28]] = NCI
    dataList[[29]] = NKI
    dataList[[30]] = PNC
    dataList[[31]] = STK
    dataList[[32]] = STNO2
    dataList[[33]] = TCGA
    dataList[[34]] = TRANSBIG
    dataList[[35]] = UCSF
    dataList[[36]] = UNC4
    dataList[[37]] = UNT
    dataList[[38]] = UPP
    dataList[[39]] = VDX

  }
  names(dataList) = infoString

  if(!is.null(dataNames))
  {
    keepInds = c()
    for(i in 1:length(dataNames))
      keepInds = c(keepInds, which(names(dataList) == dataNames[i]))
    if(length(keepInds) == 0){
      stop(cat("None of the data names supplied were present for the given cancerType and survivalMetric. Please run 1. dataList = loadMetaData(\"",cancerType, "\",\"",survivalMetric,"\") 2. validNames = obtainDataInfo(dataList, \"", survivalMetric,"\")$`Data Name` to see valid data names", sep =""))
    }else{
      remove = c(1:length(dataList))[-keepInds]
      dataList[remove] = NULL
      infoString = infoString[-remove]
    }
  }
  
  survEventList = getSurvEventData(dataList, survivalMetric)
  statList = survEventList$eventList
  timeList = survEventList$eventTimeList

  remove = c()
  #remove when all vital/tumor reccurence statuses are NA or when all daysToDeath/daysToReccurrence are NA
  for(i in 1:length(dataList))
  {
    numMissing = sum(is.na(timeList[[i]]))
    numPresent = length(timeList[[i]])
    if(numMissing == numPresent)
      remove = c(remove, i)
    numMissing = sum(is.na(statList[[i]]))
    numPresent = length(statList[[i]])
    if(numMissing == numPresent)
      remove = c(remove, i)
  }
  if(length(remove) > 0)
  {
    remove = unique(remove)
    dataList[remove] = NULL
    infoString = infoString[-remove]
    names(dataList) = infoString 
  }

  #Metabric and TCGA studies are very large so offer option to do them separately. Also would skew results if used them with other data
  if(cancerType == "breast")
  {
    if(survivalMetric == "overall" | survivalMetric == "hierarchy"){
      if(includeAll == FALSE)
      {
        metaTcgaInds = c(which(infoString == "METABRIC"), which(infoString == "TCGA"))
        dataList[metaTcgaInds] = NULL
        infoString = infoString[-metaTcgaInds] 
      }
    }else if(survivalMetric == "overallTcga"){
      dataList = list(TCGA)
      infoString = c("TCGA")
    }else if(survivalMetric == "overallMeta"){
      dataList = list(METABRIC)
      infoString = c("METABRIC")
    }
  }
  names(dataList) = infoString

  return(dataList)
}
