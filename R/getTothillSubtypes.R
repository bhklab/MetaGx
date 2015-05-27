getTothillSubtypes <- function(eset) {
  ## Load train data with predefined class labels
  #supplementary.data <- scan(system.file(file.path("extdata", "tothill.supptable.2.probesets.txt"), package="MetaGx"))
  #supplementary.probesets <- scan("../inst/extdata/tothill.supptable.2.probesets.txt", what=character(0), quiet=TRUE)
  
  # Get the probeset - Entrez ID mapping for the platform used in the Tothill et al. study
  #probe.entrez.mappings <- as.list(hgu133plus2ENTREZID[mappedkeys(hgu133plus2ENTREZID)])
  #supplementary.entrez.ids <- unlist(probe.entrez.mappings[supplementary.probesets])
  
  supp.table.2 <- read.table("../inst/extdata/tothill.supptable.probes.entrez.txt", header=TRUE)
  supplementary.entrez.ids <- unique(supp.table.2[supp.table.2$Entrez.ID != "---",]$Entrez.ID)
  
  ## Train a diagonal linear discriminant classifier using the Tothill data set and overlapping probesets / entrez gene ids
  # Get the expression matrix of this eset and the Tothill eset, with columns named by Entrez gene ids
  expression.matrix <- exprs(eset)
  rownames(expression.matrix) <- fData(eset)$EntrezGene.ID[match(rownames(expression.matrix), rownames(fData(eset)))]
 
  tothill.eset <- getGeneMapping(esets$GSE9891)
  tothill.expression.matrix <- exprs(tothill.eset)
  rownames(tothill.expression.matrix) <- fData(tothill.eset)$EntrezGene.ID[match(rownames(tothill.expression.matrix), rownames(fData(tothill.eset)))]
  
  
  colnames(tothill.expression.matrix) <- sub("X", "", pData(tothill.eset)$alt_sample_name)
  
  intersecting.entrez.ids <- intersect(supplementary.entrez.ids, intersect(rownames(expression.matrix), rownames(tothill.expression.matrix)))
  expression.matrix <- expression.matrix[rownames(expression.matrix) %in%  intersecting.entrez.ids,]
  tothill.expression.matrix <- tothill.expression.matrix[rownames(tothill.expression.matrix) %in%  intersecting.entrez.ids,]
  
  #Transpose matrices, so each row is a sample and columns are genes  
  expression.matrix <- t(expression.matrix)
  tothill.expression.matrix <- t(tothill.expression.matrix)
  
	## Train and classify with diagonal LDA based on the cla
  #supplementary.classes <- read.table(system.file(file.path("extdata", "tothill.supptable.1.classes.txt"), package="MetaGx"), header=TRUE)
  supplementary.classes <- read.table("../inst/extdata/tothill.supptable.1.classes.txt", header=TRUE)
  supplementary.classes$group <- as.character(supplementary.classes$group)
  supplementary.classes <- supplementary.classes[supplementary.classes$group != "NC",]
  supplementary.classes$group <- as.factor(supplementary.classes$group)
  levels(supplementary.classes$group) <- paste0("C", levels(supplementary.classes$group))
  rownames(supplementary.classes) <- supplementary.classes$ID
  supplementary.classes <- supplementary.classes[,-1,drop=FALSE]
  
  tothill.train.data <- merge(supplementary.classes, tothill.expression.matrix, by="row.names")
  rownames(tothill.train.data) <- tothill.train.data$Row.names
  tothill.train.data <- tothill.train.data[,-1]
  print("Training a DLDA classifier based on common genes...")
  trained.dlda <- HiDimDA::Dlda(data=tothill.train.data[,-1], grouping=tothill.train.data$group)
  print("Finished training.")
  subclasses <- predict(trained.dlda, expression.matrix, grpcodes=levels(tothill.train.data$group))$class
  pData(eset) <- data.frame(pData(eset), Tothill.subtypes=subclasses)
  return(list(Annotated.eset=eset, trained.dlda=trained.dlda))
}