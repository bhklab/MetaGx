#' Function to determine the molecular subtypes of patients in an eSet according to their gene expression
#'
#' This function determines the molecular subtypes of patients in an eSet according to their gene expression and the subtyping scheme specified by the user.
#' 
#' @param eset an eset object
#' @param cancerType a string representing the cancer type of the patients in the eset object. Must be one of "ovarian" or "breast".
#' @param subtype a string representing the subtyping scheme the patients will be classified under. For breast cancer, the options are "scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003", "intclust", "AIMS", or "claudinLow".
#' For ovarian cancer the options are "verhaak", "bentink", "tothill", "helland", "konecny", and consensusov.
#' @param intersectThresh the fraction of genes from the eset that must be present in each verhaak subtype gene set in order for the patients in the eset to be classified according to the verhaak subtypes. Default value is 0.75 and the variable is not relevant for the other subtyping schemes.
#' @return the original eset with a new subtype factor (eset$subtype) whose j'th element corresponds to the subtype for the patient in the j'th column of the eSets expression data (eset@assayData$exprs[, j])
#' @export
#' @examples
#' 
#' eset = loadMetaData("ovarian", "overall")[[1]]
#' eset$subtype
#' eset = getPatientSubtypes(eset, cancerType = "ovarian", subtype = "verhaak")
#' eset$subtype[1:10]

getPatientSubtypes = function(eset, cancerType, subtype, intersectThresh = 0.75)
{
  subtype = tolower(subtype)
  #need to adapt the reports to handle all the schemes
  if(cancerType == "breast"){
    validSchemes = tolower(c("scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003", "intclust", "AIMS", "claudinLow"))
    if(sum(grepl(subtype, validSchemes)) == 0){
      stop(paste(subtype, "is not an available subtyping scheme for breast cancer data...... Exiting the function."))
    }else{
      subtypePredictions = metaGx::getBreastSubtypes(data = eset@assayData$exprs, subtypeModel = subtype,
                                               annot = eset@featureData@data, doMapping = TRUE)
      #confirmed that order of data columns is the same as the order of the subtype results, i.e subtye result 1 is for patient in column 1
      eset$subtypes = subtypePredictions$subtype
    }
  }else if(cancerType == "ovarian"){
    validSchemes = c("verhaak", "bentink","tothill","helland","konecny", "consensusov")
    if(sum(grepl(subtype, validSchemes)) == 0){
      stop(paste(subtype, "is not an available subtyping scheme for ovarian cancer data...... Exiting the function."))
    }else{
      subtypePredictions = metaGx::getOvarianSubtypes(eset = eset, subtypeModel = subtype,
                                                      intersectThresh = 0.75)
      #confirmed that order of data columns is the same as the order of the subtype results, i.e subtye result 1 is for patient in column 1
      eset$subtypes = subtypePredictions$subtypes
    }
  }
  
  
  return(eset)
}
