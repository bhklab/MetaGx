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
## Natchar January 15, 2015
#################

# Load the curatedOvarianData package
library(curatedOvarianData)

# Load the filtering rules from patientselection.config
source(system.file("extdata", "patientselection.config", package = "curatedOvarianData"))

# Load logging package
library(logging)

# Get rid of duplicates and load esets into the environment
source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))

#columns of pdata we want
PHENOdata <- pData(esets[[1]]) [,c("sample_type", "grade", "age_at_initial_pathologic_diagnosis", "days_to_tumor_recurrence", "days_to_death", "os_binary", "relapse_binary", "batch")]

# assayData matrix
exprs <- matrix(0, nrow(fData(esets[[1]])), ncol = nrow(pData(esets[[1]])))

# phenodata in AnnotatedDataFrameß
PData <- AnnotatedDataFrame(data = PHENOdata)
FData <- AnnotatedDataFrame(data = fData(esets[[1]]))

# Rename pData to standard names
colnames(PData) <- c("tissue", "grade", "age", "t.trs", "t.os", "e.os", "e.rfs", "series")

# make the expression set with exprs and PData
neweset<- ExpressionSet(exprs, phenoData =PData, experimentData=experimentData(esets[[1]]), featureData = featureData(esets[[1]]),  protocolData= protocolData(esets[[1]]), annotation=annotation(esets[[1]]))






}

#####



