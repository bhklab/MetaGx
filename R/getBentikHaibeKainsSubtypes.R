getBentinkHaibeKainsSubtypes <- function(eset) {
	## Use predefined class labels to train diagonal LDA
	## Classify new samples
  expression.matrix <- t(exprs(eset))
  annot <- fData(eset)
  colnames(annot)[which(colnames(annot) == "EntrezGene.ID")] <- "entrezgene"
  angio <- ovcAngiogenic(data = data, annot=annot, gmap="entrezgene", do.mapping = TRUE)
  pData(eset) <- data.frame(pData(eset), Bentink.Haibe.Kains.subtypes=angio$subtype)
  return(list(Annotated.eset=eset, angio=angio))
}