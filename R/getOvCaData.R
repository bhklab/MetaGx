########################
## Benjamin Haibe-Kains
## All rights Reserved
## January 14, 2015
## Initial modification of getBrCaData for Ovarian (OvarianCuratedDataset) by  Deena and Natchar
########################

#################
## loading and changing curatedOvarianData
## Natchar January 20, 2015
#################
`getOvCaData` <- 
  function (resdir="cache", probegene.method, remove.duplicates=TRUE, topvar.genes=1000, duplicates.cor=0.98, datasets, sbt.model=c("scmgene", "scmod2", "scmod1", "pam50", "ssp2006", "ssp2003"), merging.method=c("union", "intersection"), merging.std=c("quantile", "robust.scaling", "scaling", "none"), nthread=1, verbose=TRUE) {  
    
# Load the curatedOvarianData package
library(curatedOvarianData)

# Load the filtering rules from patientselection.config
source(system.file("extdata", "patientselection.config", package = "curatedOvarianData"))

# Load logging package
library(logging)

# Get rid of duplicates and load esets into the environment
source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))


OvarianEsets <- list()
for(i in 1:length(esets)){
  currenteset <- esets[[i]]
  ################# columns of pdata we want
  PHENOdata <- pData(esets[[i]]) [,c("alt_sample_name", "sample_type", "histological_type", "summarygrade", "summarystage" , "grade", "age_at_initial_pathologic_diagnosis", "days_to_tumor_recurrence", "days_to_death", "os_binary", "relapse_binary", "batch")]
  
  
  ################# change sample_type: all "healthy" to "normal"
  # initialize lists
  
  oldList <- PHENOdata[,"sample_type"]
  listIndex <- which(oldList == "healthy")
  # listIndex <- NULL
  # for(i in 1:length(oldList)) {
  #   if (oldList[i] == "healthy") {listIndex <- c(listIndex, i)}
  # }
  
  newList <- replace(oldList, listIndex, "normal")
  
  # replacing the old sample_type column with the new column
  PHENOdata[,"sample_type"] <- newList
  
  
  ################# add dataset column 
  dataset <- matrix(names(esets)[i], nrow= nrow(pData(esets[[i]])), ncol=1)
  PHENOdata <- cbind(PHENOdata, dataset)
  
  
  
  ################# add treatment column made from pltx/tax/neo
  # # treatment column will be a list of tables 
  # 
  # pltx <- pData(currenteset)[,"pltx"]
  # tax <- pData(currenteset)[,"tax"]
  # neo <- pData(currenteset)[,"neo"]
  # 
  # treatment_values <- list(pltx, tax, neo)
  # names(treatment_values) <- c("pltx", "tax", "neo")
  
  # treatments <- list()
  # for (i in 1:length(pltx)){
  #   treatments[i]<- data.frame(treatment_values[[1]][i], treatment_values[[2]][i], treatment_values [[3]][i])
  # }
  
  # 
  # treatment_table <- cbind(pltx, tax, neo)
  # treatment_values <- NULL
  # for(i in 1:length(pltx)){
  #   treatment_values <- list(treatment_values, treatment_table[i,])
  # }
  
  
  ################# add treatment column made from pltx/tax/neo
  # save columns as a table lists of "y", "n" or NA, if mix of n and NA, produce n
  
  pltx <- pData(esets[[i]])[,"pltx"]
  tax <- pData(esets[[i]])[,"tax"]
  neo <- pData(esets[[i]])[,"neo"]
  
  treatment_table <- cbind(pltx, tax, neo)
  row <- NULL
  treatment_values <- NULL
  single_value <- NULL
  naindex <-NULL
  yindex <-NULL
  IndexedRow <- NULL
  
  for(k in 1:length(pltx)){
    row <- treatment_table[k,]
    # map NA values to "na" to avoid missing value errors
    naindex <- which(is.na(row))
    IndexedRow <- replace(row, naindex, "na")
    
    # if they are all NA, put NA, if all n, put n, if all y, then list all of them
    if (all(IndexedRow == IndexedRow[1])){
      if ( IndexedRow[1] == "n"){ #all "n"
        single_value <- IndexedRow[1]
      } 
      else if (IndexedRow[1] == "na"){ #all NA
        single_value <- NA
      }
      else {# all are y!
        single_value <- paste(as.character(colnames(treatment_table)), collapse = "/")
      }
      
      ### if entries are not all the same, take y's  
    } else {
      #refresh single_value
      single_value <- NULL 
      # if none of the elements are y, then mix of n and NA
      if (all(which(row == "y")) ==0){ 
        single_value <- "n"
        #       single_value <- NULL
      } else {# y's are present
        yindex <- which(row == "y")
        row <- NULL # reset row to clean
        for(j in 1:length(yindex)){ # go through yindex
          row <- c(row, colnames(treatment_table)[j])
          single_value <- paste(as.character(row), collapse= "/") # example: produces pltx/neo
        }
      }
    }
    treatment_values <- c(treatment_values, single_value)
  }
  
  #### adding treatment column, adding platform
  PHENOdata <- cbind(PHENOdata, treatment_values)
  
  platform_values <- matrix(annotation(esets[[i]]), nrow= nrow(pData(esets[[i]])), ncol=1)
  PHENOdata <- cbind(PHENOdata, platform_values)
  
  # phenodata and featuredata in AnnotatedDataFrame
  PData <- AnnotatedDataFrame(data = PHENOdata)
  Exprs <- exprs(currenteset)
  
  ################# Rename pData to standard names
  colnames(PData) <- c("samplename", "tissue", "histological_type","summarygrade", "summarystage","grade", "age", "t.rfs", "t.os", "e.os", "e.rfs", "series", "dataset", "treatment", "platform")
  
  
  ################# make the expression set with exprs and PData
  neweset <- list()
  neweset[1]<- ExpressionSet(Exprs, phenoData =PData, experimentData=experimentData(currenteset), featureData = featureData(currenteset),  protocolData= protocolData(currenteset), annotation=annotation(currenteset))
  
  OvarianEsets[i] <- neweset
}
 return (OvarianEsets)
}
