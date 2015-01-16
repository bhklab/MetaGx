########################
## Benjamin Haibe-Kains
## All rights Reserved
## January 14, 2015
## Initial modification of getBrCaData for Ovarian (OvarianCuratedDataset) by  Deena and Natchar
########################

`getOvCaData` <- 
function (resdir="cache", probegene.method, remove.duplicates=TRUE, topvar.genes=1000, duplicates.cor=0.98, datasets, sbt.model=c("scmgene", "scmod2", "scmod1", "pam50", "ssp2006", "ssp2003"), merging.method=c("union", "intersection"), merging.std=c("quantile", "robust.scaling", "scaling", "none"), nthread=1, verbose=TRUE) {  

#Need to remove InSilicoDB reference
# Keep Normalization, merging?
  stop("getOvCaData is not implemented yet")



#################
## loading and changing curatedOvarianData
## Natchar January 16, 2015
#################

# Load the curatedOvarianData package
library(curatedOvarianData)

# Load the filtering rules from patientselection.config
source(system.file("extdata", "patientselection.config", package = "curatedOvarianData"))

# Load logging package
library(logging)

# Get rid of duplicates and load esets into the environment
source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))


################# columns of pdata we want
PHENOdata <- pData(esets[[1]]) [,c("sample_type", "grade", "age_at_initial_pathologic_diagnosis", "days_to_tumor_recurrence", "days_to_death", "os_binary", "relapse_binary", "batch")]


################# change sample_type: all "healthy" to "normal"
# initialize lists

oldList <- PHENOdata[,"sample_type"]
listIndex <- NULL
for(i in 1:length(oldList)) {
  if (oldList[i] == "healthy") {listIndex <- c(listIndex, i)}
}

newList <- replace(oldList, listIndex, "normal")

# replacing the old sample_type column with the new column
PHENOdata[,"sample_type"] <- newList
  
  
################# add dataset column 
dataset <- matrix(names(esets)[1], nrow= nrow(pData(esets[[1]])), ncol=1)
PHENOdata <- cbind(PHENOdata, dataset)


################# add treatment column made from pltx/tax/neo
# save columns as a table lists of "y", "n" or NA

pltx <- pData(esets[[1]])[,"pltx"]
tax <- pData(esets[[1]])[,"tax"]
neo <- pData(esets[[1]])[,"neo"]

treatment_table <- cbind(pltx, tax, neo)
row <- NULL
treatment_values <- NULL
single_value <- NULL
yindex <-NULL
for(i in 1:length(pltx)){
  row <- treatment_table[i,]
  # if they are all NA, put NA, if all n, put n, if all y, then list all of them
  if ((all(row == row[1])) == TRUE || all(lapply(row, is.na) == TRUE)){
    if (is.na(row[1]) || row[1] == "n"){ #all NA, all "n"
      single_value <- row[1]
    } else {# all are y!
      single_value <- paste(as.character(colnames(treatment_table)), collapse = "/")
    }
    
  ### if entries are not all the same, take y's  
  } else {
    #refresh single_value
    single_value <- NULL 
    # if none of the elements are y, then mix of n and NA
    if (which(row == "y") ==0){ 
      single_value <- "n"
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

#### adding treatment column

PHENOdata <- cbind(PHENOdata, treatment_values)

# assayData matrix
exprs <- matrix(0, nrow(fData(esets[[1]])), ncol = nrow(PHENOdata))


# phenodata and featuredata in AnnotatedDataFrame
PData <- AnnotatedDataFrame(data = PHENOdata)
FData <- AnnotatedDataFrame(data = fData(esets[[1]]))

################# Rename pData to standard names
colnames(PData) <- c("tissue", "grade", "age", "t.trs", "t.os", "e.os", "e.rfs", "series", "dataset", "treatment")


################# make the expression set with exprs and PData
neweset<- ExpressionSet(exprs, phenoData =PData, experimentData=experimentData(esets[[1]]), featureData = featureData(esets[[1]]),  protocolData= protocolData(esets[[1]]), annotation=annotation(esets[[1]]))
