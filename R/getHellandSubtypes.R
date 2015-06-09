getHellandSubtypes <- function(eset) {
  supplementary.type.1 <- read.xls("../inst/extdata/journal.pone.0018064.s015.XLS", sheet=1)
  supplementary.type.2 <- read.xls("../inst/extdata/journal.pone.0018064.s015.XLS", sheet=2)
  supplementary.type.4 <- read.xls("../inst/extdata/journal.pone.0018064.s015.XLS", sheet=3)
  supplementary.type.5 <- read.xls("../inst/extdata/journal.pone.0018064.s015.XLS", sheet=4)
  supplementary.tables <- list(C1=supplementary.type.1, C2=supplementary.type.2, C4=supplementary.type.4, C5=supplementary.type.5)
  
  entrez.id.logFC.list <- lapply(supplementary.tables, function(x) {
    ## Use the supplementary table's listed probe id and gene name to determine the Entrez ID
    # If there is only one EntrezID that maps to a probe in hgu133plus2.db, use that Entrez ID.
    # If there are multiple EntrezIDs that map to a probe, then use the EntrezID (if any) that corresponds to the provided gene symbol.
    current.mapping <- suppressWarnings(select(hgu133plus2.db, as.character(x$ID), c("ENTREZID", "SYMBOL")))
    current.mapping <- current.mapping[ !is.na(current.mapping$ENTREZID), ]
    colnames(x)[1:2] <- c("PROBEID", "SYMBOL")
    mappings.with.unique.probeid <- current.mapping[ !(current.mapping$PROBEID %in% current.mapping$PROBEID[duplicated(current.mapping$PROBEID)]),]
    mappings.with.duplicate.probeid <- current.mapping[ current.mapping$PROBEID %in% current.mapping$PROBEID[duplicated(current.mapping$PROBEID)],]
    mappings.with.duplicate.probeid <- merge(x, mappings.with.duplicate.probeid, by=c("PROBEID", "SYMBOL"))[, c("PROBEID", "ENTREZID", "SYMBOL")]
    mappings.with.duplicate.probeid <- unique(mappings.with.duplicate.probeid)
    current.mapping <- rbind(mappings.with.unique.probeid, mappings.with.duplicate.probeid)
    to.return <- merge(x, current.mapping, by="PROBEID")[, c("ENTREZID", "PROBEID", "logFC")]
    return(to.return)
    })
  
  gene.union.in.supplementary <- Reduce(function(x,y) union(x, y), lapply(entrez.id.logFC.list, function (x) x$ENTREZID))
  
  intersecting.entrez.ids <- intersect(gene.union.in.supplementary, fData(eset)$EntrezGene.ID)
  
  # Only keep genes present in both the supplementary and this eset
  entrez.id.logFC.list <- lapply(entrez.id.logFC.list, function(x) x[x$ENTREZID %in% intersecting.entrez.ids, ])
  
  subtype.scores <- sapply(entrez.id.logFC.list, function(x) {
    ordered.expression.subset <- exprs(eset)[match(x$ENTREZID, fData(eset)$EntrezGene.ID),]
    
    return(apply(ordered.expression.subset, 2, function(y) sum((y * x$logFC))))
    })
  old.rownames <- rownames(subtype.scores)
  # Scale to mean=0, variance=1
  subtype.scores <- apply(subtype.scores, 2, scale)
  rownames(subtype.scores) <- old.rownames
  subclasses <- factor(colnames(subtype.scores)[apply(subtype.scores, 1, which.max)], levels=names(supplementary.tables))
  eset$Helland.subtypes <- subclasses
  return(list(Annotated.eset=eset, subtype.scores=subtype.scores))
}