library(Biobase)
library(NMF)
library(clue)

task.id <- as.integer(Sys.getenv("SGE_TASK_ID"))
set.seed(660 + task.id * 100)

.getFilteredEsetByMAD <- function(eset, num.genes) {
  expression.matrix <- exprs(eset)
  mad.vals <- apply(expression.matrix, 1, mad)
  # expression.matrix has genes as rows, patients as columns
  eset <- eset[mad.vals >= tail(sort(mad.vals),num.genes)[1],]
  return(eset)
}

.getFilteredEsetByGeneList <- function(eset, gene.list, list.type=c("gene.symbol", "entrez.id")) {
  list.type <- match.arg(list.type)
  if(list.type == "gene.symbol") {
    eset <- eset[fData(eset)$gene %in% gene.list,]
  } else if(list.type == "entrez.id") {
    eset <- eset[fData(eset)$EntrezGene.ID %in% gene.list,]
  }
  return(eset)
}


.getNMFClasses <- function(eset, filter.genes=TRUE, num.genes=2000, rank=4, nrun=100) {
  # rescale eset by z-score per gene
  expression.matrix <- exprs(eset)
  if(filter.genes) {
    mad.vals <- apply(exprs(eset), 1, mad)
    expression.matrix <- exprs(eset)[mad.vals >= tail(sort(mad.vals),num.genes)[1],]
  }
  if(any(expression.matrix < 0)) {
    expression.matrix <- expression.matrix + abs(min(expression.matrix))
  }
  
  nmf.out <- nmf(expression.matrix, rank=rank, nrun=nrun)
  # Clustering using consensus
  classes <- cutree(consensushc(nmf.out, dendrogram=FALSE), k=rank)
  
  classes <- as.factor(paste0("ConsensusNMF_", classes))
  # Clustering using matrix factorization
  #h.mat <- coef(nmf.out)
  #classes <- apply(h.mat, 2, which.max)
  #classes <- as.factor(paste0("NMF_", classes))
  return(classes)
}

.getConsensusKMeansClasses <- function(eset, filter.genes=TRUE, num.genes=2000, k=4, num.iterations=100) {
  # rescale eset by z-score per gene
  expression.matrix <- exprs(eset)
  if(filter.genes) {
    mad.vals <- apply(exprs(eset), 1, mad)
    expression.matrix <- exprs(eset)[mad.vals >= tail(sort(mad.vals),num.genes)[1],]
  }
  expression.matrix <- t(expression.matrix)
  kmeans.out <- lapply(1:100, function(x) kmeans(expression.matrix, centers = k))
  cl.ensemble <- cl_ensemble(list=kmeans.out)
  consensus.out <- cl_consensus(cl.ensemble)
  membership.matrix <- matrix(as.vector(consensus.out$.Data), nrow=ncol(eset))
  
  classes <- apply(membership.matrix, 1, which.max)
  classes <- as.factor(paste0("kmeans_", classes))
  return(classes)
}

load("esets.not.rescaled.RData")

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

#config.grid <- expand.grid(
#  gene.set=c("tcga", "tothill", "konecny", "1000", "1500", "2000", "2500", "3000"),
#  algorithm=c("nmf", "kmeans"),
#  dataset.index=1:16,
#  k=4)

config.grid <- expand.grid(
  gene.set=c("tcga", "tothill", "konecny"),
  algorithm=c("nmf", "kmeans"),
  dataset.index=1:16,
  k=4)

config.grid <- data.frame(
  gene.set=rep(c("tcga", "tothill", "konecny"), each=16),
  algorithm=rep(c("nmf", "kmeans", "nmf"), each=16),
  dataset.index=rep(1:16,3),
  k=4
  )

config.id <- task.id %% nrow(config.grid)
if(config.id == 0) {
  config.id <- nrow(config.grid)
}

config.grid$gene.set <- as.character(config.grid$gene.set)
config.grid$algorithm <- as.character(config.grid$algorithm)
config.grid$dataset.index <- as.integer(as.character(config.grid$dataset.index))
config.grid$k <- as.integer(as.character(config.grid$k))

gene.set <- config.grid$gene.set[config.id]
algorithm <- config.grid$algorithm[config.id]
dataset.index <- config.grid$dataset.index[config.id]
k <- config.grid$k[config.id]

out.dir <- paste0("aug8clusters/", gene.set, "_", algorithm, "_", k)
dir.create(out.dir)

current.eset <- esets.not.rescaled[[dataset.index]]
current.eset.name <- names(esets.not.rescaled)[dataset.index]

if(gene.set == "tcga") {
  tcga.matrix <- read.delim("../../inst/extdata/TCGA_489_UE.top1500.txt", sep="\t", header=TRUE)
  tcga.gene.symbols <- rownames(tcga.matrix)
  current.eset <- .getFilteredEsetByGeneList(current.eset, tcga.gene.symbols, list.type="gene.symbol")
} else if(gene.set == "tothill") {
  tothill.gene.table <- read.table("../../inst/extdata/tothill.clustering.genes.txt")
  current.eset <- .getFilteredEsetByGeneList(current.eset, as.character(tothill.gene.table$entrez.id), list.type="entrez.id")
} else if (gene.set == "konecny"){
  konecny.unique.entrez.ids <- scan("../../inst/extdata/konecny.unique.entrez.ids.txt", what=character(0))
  current.eset <- .getFilteredEsetByGeneList(current.eset, konecny.unique.entrez.ids, list.type="entrez.id")
}else {
  num.genes <- as.integer(gene.set)
  current.eset <- .getFilteredEsetByMAD(current.eset, num.genes)
}

if(algorithm == "kmeans") {
  out.classes <- .getConsensusKMeansClasses(current.eset, filter.genes=FALSE, k=k, num.iterations=1000)
} else if(algorithm == "nmf") {
  out.classes <- .getNMFClasses(current.eset, filter.genes=FALSE, rank=k, nrun=100)
}

write(as.character(out.classes), file=paste0(out.dir, "/", current.eset.name, "_", task.id, "_classes.txt"), ncolumns = 1)
