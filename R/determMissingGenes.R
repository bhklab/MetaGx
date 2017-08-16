
#' Function to determine which genes are present and absent from the datasets
#'
#' This function checks each dataset to determine whether the genes supplied are present for the patients
#' @param geneIds A character vector containing the ensemble IDs, entrez IDs, or gene symbols for the genes in the signature being tested
#' @param dataList a list containing the datasets that the genes will be searched for in
#' @return a data frame with information pertaining to whether the genes were found in the supplied datasets, note that if the ID could not be mapped
#' to an entrez ID the supplied ID will be found in both the "Symbols of Missing Genes" and "Entrez IDs of Missing Genes" section as its ID is unknown   
#' @export
#' @importFrom matrixStats rowIQRs
#' @importFrom GSVA gsva
#' @importFrom survcomp censor.time combine.est D.index km.coxph.plot
#' @importFrom survival Surv
#' @importFrom forestplot forestplot fpTxtGp fpColors
#' @importFrom grid gpar unit
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @examples
#' dataList = loadMetaData("ovarian", "overall")
#' geneIds = c("KLK14","RHOX8","ADAMTS20","IDO1")
#' genesMissFrame = determMissingGenes(geneIds, dataList)
#'

#want a frame that has columns dataset name, number present, number missing genes, names of missing genes, entrez of missing genes, maybe prognostic value of those genes also?
determMissingGenes = function(geneIds, dataList)
{
  #geneEntrezIds = geneEntrezList[[1]]
  geneFrame = getGeneInfo(geneIds, rep(1, length(geneIds)))
  geneEntrezIds = geneFrame$`Entrez ID`
  geneType = "ENTREZID"
  mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneEntrezIds), columns=c("SYMBOL"))
  mapping = mapping[!duplicated(mapping[,1]),]
  geneSymbols = mapping[, 2]

  numbMissVec = c()
  missSymbVec = c()
  missEntVec = c()
  numbPresVec = c()
  for(i in 1:length(dataList))
  {
    numbMiss = 0
    missSymbStr = ""
    missEntStr = ""
    dataInfo = dataList[[i]]@featureData@data;
    dataInfoEntrez = as.character(dataInfo$EntrezGene.ID)

    dataGeneInds = c()
    for(j in 1:length(geneEntrezIds))
    {
      genePresent = sum(dataInfoEntrez == geneEntrezIds[j]) > 0
      if(is.na(geneEntrezIds[j]))
      {
        genePresent = FALSE 
        geneSymbols[j] = paste("ID", geneIds[j])
        geneEntrezIds[j] = paste("ID", geneIds[j])
      }

      if(!genePresent)
      {
        numbMiss = numbMiss + 1
        missSymbStr = paste0(missSymbStr, geneSymbols[j], ",")
        missEntStr = paste0(missEntStr, geneEntrezIds[j], ",")
      }

      dataGeneInds = c(dataGeneInds, which(dataInfoEntrez == geneEntrezIds[j]));
    }
    if(grepl(",", substr(missEntStr, nchar(missEntStr), nchar(missEntStr))))
    {
      missEntStr = substr(missEntStr, 1, nchar(missEntStr) - 1)
      missSymbStr = substr(missSymbStr, 1, nchar(missSymbStr) - 1)
    }
    if(missSymbStr == "")
      missSymbStr = "NA"
    if(missEntStr == "")
      missEntStr = "NA"
    numbMissVec = c(numbMissVec, numbMiss)
    numbPresVec = c(numbPresVec, length(geneEntrezIds) - numbMiss)
    missSymbVec = c(missSymbVec, missSymbStr)
    missEntVec = c(missEntVec, missEntStr)
  }

  geneMissFrame = as.data.frame(cbind(names(dataList), numbPresVec, numbMissVec, missSymbVec, missEntVec))
  colnames(geneMissFrame) = c("Dataset", "Genes Present", "Genes Missing", "Symbols of Missing Genes", "Entrez IDs of Missing Genes")
  return(geneMissFrame)
}



