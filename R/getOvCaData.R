########################
## Benjamin Haibe-Kains
## All rights Reserved
## January 14, 2015
## Initial modification of getBrCaData for Ovarian (OvarianCuratedDataset) by  Deena and Natchar
########################
#################
## loading and changing curatedOvarianData
## Natchar January 30, 2015
#################
`getOvCaData` <- 
  function (resdir="cache", probegene.method, remove.duplicates=TRUE, topvar.genes=1000, duplicates.cor=0.98, datasets, sbt.model=c("scmgene", "scmod2", "scmod1", "pam50", "ssp2006", "ssp2003"), merging.method=c("union", "intersection"), merging.std=c("quantile", "robust.scaling", "scaling", "none"), nthread=1, verbose=TRUE) {  
    
# Load the curatedOvarianData package
library(curatedOvarianData)

# Load logging package
library(logging)

# Load org.HS.eg.db package
library(org.Hs.eg.db)

# Load the filtering rules from patientselection.config
source(system.file("extdata", "patientselection.config", package = "curatedOvarianData"))

# Get rid of duplicates and load esets into the environment, list of expression sets esets
source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))

########################################## Manipulate and change esets ####################

OvarianEsets <- list()
for(i in 1:length(esets)){
  currenteset <- esets[[i]]
  
  #######################################
  ######### pData to include columns in common with BreastData
  #######################################
  
  
  ################# columns of pdata we want
  PHENOdata <- pData(currenteset) [,c("alt_sample_name", "sample_type", "histological_type", "summarygrade", "summarystage" , "grade", "age_at_initial_pathologic_diagnosis", "days_to_tumor_recurrence", "days_to_death", "vital_status","os_binary", "relapse_binary", "batch")]
  
  
  ################# change sample_type: all "healthy" to "normal"
  # initialize lists
  
  oldList <- PHENOdata[,"sample_type"]
  listIndex <- which(oldList == "healthy")
  
  newList <- replace(oldList, listIndex, "normal")
  
  # replacing the old sample_type column with the new column
  PHENOdata[,"sample_type"] <- newList
  
  ################# change vital_status from living/deceased to 1 0
  oldList <- PHENOdata[,"vital_status"]
  listIndex <- which(oldList =="living")
  newList <- replace(oldList, listIndex, 1)
  listIndex <- which(oldList == "deceased")
  newList <- replace(newList, listIndex, 0)
  
  PHENOdata[,"vital_status"] <- newList
  
  
  ################# add dataset column 
  dataset <- matrix(names(esets)[i], nrow= nrow(pData(currenteset)), ncol=1)
  PHENOdata <- cbind(PHENOdata, dataset)
  
  ################# add treatment column made from pltx/tax/neo
  # save columns as a table lists of "y", "n" or NA, if mix of n and NA, produce n
  
  pltx <- pData(currenteset)[,"pltx"]
  tax <- pData(currenteset)[,"tax"]
  neo <- pData(currenteset)[,"neo"]
  
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
  
  platform_values <- matrix(annotation(currenteset), nrow= nrow(pData(currenteset)), ncol=1)
  PHENOdata <- cbind(PHENOdata, platform_values)
  
  #### adding platform2 (GPL)
  platform2 <-matrix(experimentData(currenteset)@other$platform_accession, nrow = nrow(pData(currenteset)), ncol=1)
  PHENOdata <- cbind(PHENOdata, platform2)
  
  
  # phenodata and featuredata in AnnotatedDataFrame
  PData <- AnnotatedDataFrame(data = PHENOdata)
  Exprs <- exprs(currenteset)
  FData <- AnnotatedDataFrame(fData(currenteset))
  
  ################# Rename pData to standard names
  colnames(PData) <- c("samplename", "tissue", "histological_type","summarygrade", "summarystage","grade", "age", "t.rfs", "t.os", "e.os","os_binary", "e.rfs", "series", "dataset", "treatment", "platform1", "platform2")
  
  
  ################# make the expression set with exprs and PData
  neweset <- list()
  neweset[1]<- ExpressionSet(Exprs, phenoData =PData, experimentData=experimentData(currenteset), featureData = FData,  protocolData= protocolData(currenteset), annotation=annotation(currenteset))

  OvarianEsets[i] <- neweset

  #######################################
  ######## fData to include entrezgeneID
  #######################################
  
  ################# add etrez gene id to fData
  entrezgene <- list()
  gs <- toTable(org.Hs.egSYMBOL)
  gs <- gs[!is.na(gs[ , "symbol"]) & !duplicated(gs[ , "symbol"]), , drop=FALSE]
  gs <- gs[fData(OvarianEsets[[i]])[,"gene"], "gene_id"]

  fData(OvarianEsets[[i]])$entrezgene <- gs
  rownames(fData(OvarianEsets[[i]])) <- fData(OvarianEsets[[i]])[,"probeset"]

  
  #######################################
  ######## Experiment Data to include cancer_type, Angiogenic class, fuzzy, crisp
  #######################################
  
  
  ################# Specify cancer_type in experimentData
  
  experimentData(OvarianEsets[[i]])@other$cancer_type <- "ovarian"
  
  ################# ovcAngiogenic (Angiogenic vs non-Angiogenic Subtyping)
  data <- t(exprs(OvarianEsets[[i]]))
  annot <- fData(OvarianEsets[[i]])
  colnames(data) <- rownames(annot)
  hgs <- vector()
  for (l in 1:length(pData(OvarianEsets[[i]])$summarygrade)){
   hgs[l] <- (pData(OvarianEsets[[i]])$summarygrade[l] == "high" || pData(OvarianEsets[[1]])$summarystage[l] == "late" || pData(OvarianEsets[[1]])$histological_type[l] == "ser")
  }
  angio <- ovcAngiogenic(data = data, annot=annot, hgs=hgs, gmap="entrezgene", do.mapping = TRUE)
  experimentData(OvarianEsets[[i]])@other$class <- angio$subtype$subtype
  experimentData(OvarianEsets[[i]])@other$fuzzy <- angio$subtype$Angiogenic.proba
  experimentData(OvarianEsets[[i]])@other$crisp <- list()
  
  entrezgene <- NULL
  for(entrez in 1:nrow(fData(OvarianEsets[[i]]))){
    entrezgene <- c(entrezgene, paste("geneid", as.character(fData(OvarianEsets[[i]])[entrez,"entrezgene"] ), sep="."))
    
  }
  rownames(fData(OvarianEsets[[i]])) <- entrezgene
  
}
 names(OvarianEsets) <- names(esets)
 return (OvarianEsets)


## End of getOvCaData
}
