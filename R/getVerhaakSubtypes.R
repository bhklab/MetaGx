getVerhaakSubtypes <- function(eset) {
  ## Load gene sets from the original publication
  # Load Verhaak et al. supplementary from the package inst directory
	#supplementary.data <- read.xls(system.file(file.path("extdata", "JCI65833sd1.xls"), package="MetaGx"), sheet=7, skip=1)
  # Use this instead when running this method from source
	supplementary.data.sheet7 <- read.xls(system.file(file.path("extdata", "JCI65833sd1.xls", package="MetaGx"), sheet=7, skip=1))
	supplementary.data.sheet1 <- read.xls(system.file(file.path("extdata", "JCI65833sd1.xls", package="MetaGx"), skip=1))
  
	genesets <- lapply(levels(supplementary.data.sheet7$CLASS), function(y) as.character(supplementary.data.sheet7[supplementary.data.sheet7$CLASS==y,1]))
	names(genesets) <-  levels(supplementary.data.sheet7$CLASS)
	
	# For ssGSEA scores for the new samples, use the intersecting genes
	genesets <- lapply(genesets, function(x) intersect(x, fData(eset)$gene) )
  
  ## Determine the ssGSEA cutoffs for the IMR and MES subtypes
	supplementary.tcga.discovery <- supplementary.data.sheet1[ supplementary.data.sheet1$DATASET=="TCGA-discovery", 
                                                             c("ID", "SUBTYPE") ]
	supplementary.tcga.discovery <- supplementary.tcga.discovery[ supplementary.tcga.discovery$SUBTYPE %in% c("Mesenchymal", "Immunoreactive"), ]
  
  #tcga.eset <- esets$TCGA
  #tcga.eset <- tcga.eset[,tcga.eset$unique_patient_ID %in% supplementary.tcga.discovery$ID]
  #tcga.expression.matrix <- exprs(tcga.eset)
  #rownames(tcga.expression.matrix) <- fData(tcga.eset)$gene
  #tcga.gsva.out <- gsva(tcga.expression.matrix, genesets, method="ssgsea", min.sz=10, tau=0.75, parallel.sz=4)
  #tcga.gsva.out <- as.data.frame(t(tcga.gsva.out))
  #tcga.gsva.out$ID = tcga.eset$unique_patient_ID
	
  #tcga.gsva.out.with.published.subtype <- merge(tcga.gsva.out, supplementary.tcga.discovery, by="ID")
  IMR.threshold <- 0.63 #min(tcga.gsva.out.with.published.subtype$IMR[ tcga.gsva.out.with.published.subtype$SUBTYPE=="Immunoreactive" ])
  MES.threshold <- 0.56 #min(tcga.gsva.out.with.published.subtype$MES[ tcga.gsva.out.with.published.subtype$SUBTYPE=="Mesenchymal" ])
 
	expression.matrix <- exprs(eset)
  rownames(expression.matrix) <- fData(eset)$gene
  
  ## Get ssGSEA subtype scores
	gsva.out <- gsva(expression.matrix, genesets, method="ssgsea", tau=0.75, parallel.sz=4, mx.diff=FALSE, ssgsea.norm=FALSE)
  gsva.out <- t(gsva.out)
  
  gsva.out <- apply(gsva.out, 2, function(x) ( x - min(x) ) / ( max(x) - min(x) ))
  
  ## Classify each sample according to the max ssGSEA subtype score, using the scheme provided in the methods.
  
  subclasses <- apply(gsva.out, 1, function(x) {
    if(x[which(colnames(gsva.out)=="IMR")] > IMR.threshold && x[which(colnames(gsva.out)=="MES")] > MES.threshold) {
      return(c("IMR", "MES")[which.max(x[c("IMR", "MES")])])
      }
    colnames(gsva.out)[which.max(x)]
    })
  
  subclasses <- factor(subclasses, levels=levels(supplementary.data.sheet7$CLASS))
  ## Append a new column for Verhaak subtypes
  eset$Verhaak.subtypes <- subclasses
  return(list(Annotated.eset=eset, gsva.out=gsva.out))
}
