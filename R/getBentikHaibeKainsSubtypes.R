getBentinkHaibeKainsSubtypes <- function(eset) {
	## Classify new samples
  expression.matrix <- t(exprs(eset))
  annot <- fData(eset)
  colnames(annot)[which(colnames(annot) == "EntrezGene.ID")] <- "entrezgene"
  angio <- ovcAngiogenic(data = expression.matrix, annot=annot, gmap="entrezgene", do.mapping = TRUE)
  pData(eset) <- data.frame(pData(eset), Bentink.Haibe.Kains.subtypes=angio$subtype$subtype)
  return(list(Annotated.eset=eset, angio=angio))
}