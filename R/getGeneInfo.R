
#' Function to create a data frame with information pertaining to the signature
#'
#' This function maps the geneIds supplied to ensemble IDs, entrez IDs, and gene symbols and returns a data frame with this information as well as information about the genes, when available
#' @param geneIds A character vector containing the ensemble IDs, entrez IDs, or gene symbols for the genes in the signature being tested
#' @param geneDirecs A numeric vector composed of +1 and - 1 indicating the direction of association for each gene supplied in the geneIds vector. 
#' This is an optional argument for those that would like to append the direction of association of the genes to the dataframe (default is NULL).
#' @return a data frame with information about the genes supplied
#' @importFrom matrixStats rowIQRs
#' @importFrom GSVA gsva
#' @importFrom survcomp censor.time combine.est D.index km.coxph.plot
#' @importFrom survival Surv
#' @importFrom forestplot forestplot fpTxtGp fpColors
#' @importFrom grid gpar unit
#' @importFrom mclust estep estepE estepEEI
#' @export
#' @examples
#'
#' geneIds = c("KLK14","RHOX8","ADAMTS20","IDO1")
#' geneDirecs = c(1, 1, -1, -1)
#' geneInfoFrame = getGeneInfo(geneIds, geneDirecs)
#' geneInfoFrame
#'

getGeneInfo = function(geneIds, geneDirec = NULL)
{
  #geneIds = genesOfInt$Gene.Symbol
  #geneDirec = genesOfInt$t
  #inGeneSig = genesOfInt$geneInSig
  library(org.Hs.eg.db)
  geneDirecOrig = geneDirec
  if(is.null(geneDirecOrig))
    geneDirec = rep(1, length(geneIds))
  
  geneDirec = sign(geneDirec)
  geneIds = as.character(geneIds)
  geneDataDf = as.data.frame(cbind(geneIds, geneDirec), stringsAsFactors=FALSE)
  colnames(geneDataDf) = c("geneIds", "geneDirecs")

  if(suppressWarnings(sum(is.na(as.numeric(geneIds))) < length(geneIds)/2))
  {
    geneType = "ENTREZID"
    names(geneDataDf)[1] = "entrez"
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$entrez), columns=c("ENSEMBL"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$ensemble = mapping[, 2]
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$entrez), columns=c("SYMBOL"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$symbol = mapping[, 2]
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$entrez), columns=c("GENENAME"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$geneDescrip = mapping[, 2]
  } else if(sum(sum(grepl("ENSG00", geneIds))) > length(geneIds)/2){
    geneType  = "ENSEMBL"
    names(geneDataDf)[1] = "ensemble"
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$ensemble), columns=c("ENTREZID"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$entrez = mapping[, 2]
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$ensemble), columns=c("SYMBOL"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$symbol = mapping[, 2]
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$ensemble), columns=c("GENENAME"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$geneDescrip = mapping[, 2]
  }else{
    geneType = "SYMBOL"
    names(geneDataDf)[1] = "symbol"
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$symbol), columns=c("ENTREZID"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$entrez = mapping[, 2]
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$symbol), columns=c("ENSEMBL"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$ensemble = mapping[, 2]
    mapping = AnnotationDbi::select(org.Hs.eg.db, keytype=geneType, keys=as.character(geneDataDf$symbol), columns=c("GENENAME"))
    mapping = mapping[!duplicated(mapping[,1]),]
    geneDataDf$geneDescrip = mapping[, 2]
  }
  geneDataDf = as.data.frame(cbind(geneDataDf$symbol, geneDataDf$geneDirecs, geneDataDf$entrez, geneDataDf$ensemble, geneDataDf$geneDescrip), stringsAsFactors=FALSE)
  colnames(geneDataDf) = c("Symbol", "Direction", "Entrez ID", "Ensemble ID", "Description")
  
  if(is.null(geneDirecOrig))
    geneDataDf = geneDataDf[, -which(colnames(geneDataDf) == "Direction")]
  return(geneDataDf)

}
