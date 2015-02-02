########################
## Benjamin Haibe-Kains
## All rights Reserved
## January 14, 2015
## Initial modification of getBrCaData for Ovarian (OvarianCuratedDataset) by  Deena and Natchar
########################

## using sigOvcAngiogenic

sigOvcAngiogenic <- sigOvcAngiogenic[complete.cases(sigOvcAngiogenic),]
sigs <- data.frame(cbind(sigOvcAngiogenic[,"entrezgene"], sigOvcAngiogenic[,"weight"]))
sigs <- data.frame(paste("geneid", sigs[,1], sep="."), sigs)
colnames(sigs) <- c("probe", "EntrezGene.ID", "coefficient")


######### removing duplicate entrezgenes, only keeping the highest coefficient (otherwise duplicate row.names error occurs due to duplicates from multiple functions)
sigs <- do.call(rbind, 
        by(sigs, INDICES=list(sigs$probe), 
           FUN=function(x) head(x[order(x$coefficient), ], 1)))

rownames(sigs) <- sigs[,"probe"]

sigs <- c(list("DONE_SIG"=sigs), sigs <- lapply(apply(sigs, 1, list), function (x) { return (data.frame("probe"=x[[1]]["probe"], "EntrezGene.ID"=as.character(as.numeric(x[[1]]["EntrezGene.ID"])), "coefficient"=as.numeric(x[[1]]["coefficient"]))) }))


ss <- NULL
for (i in 1:length(sigs)) {
  ss <- c(ss, list(data.frame("feature"=paste("geneid", as.character(sigs[[i]][complete.cases(sigs[[i]]), "EntrezGene.ID"]), sep="."), "coefficient"=sigs[[i]][complete.cases(sigs[[i]]), "coefficient"])))
}

names(ss) <- c(names(sigs))


##### number of warnings equal to number of NA's - Feature Identifiers from the signature do not map with those of the expressionSet object - origin: sigScore.R


##### calling subtypeAssociation

for (i in 1:length(OvarianEsets)){
  MetaGx2::subtypeAssociation(eset=OvarianEsets[[i]], sig= ss, plot=TRUE, weighted=FALSE, condensed=TRUE, resdir=sprintf("saveres_%s", names(OvarianEsets)[i] ),  method="weighted.average", scaling="robust")
}

MetaGx2::subtypeAssociation(eset=OvarianEsets[[16]], sig= ss, plot=TRUE, weighted=FALSE, condensed=TRUE, resdir="saveres",  method="weighted.average", scaling="robust")

#### calling Correlation works
subtypeCorrelation(eset=OvarianEsets[[16]], sig=ss, weighted=TRUE, condensed=TRUE, resdir="saveres", method="spearman", sig.method="weighted.average", sig.scaling="robust")

#### calling Survival works, but gives empty xls (subtype_dindex , subtype_cindex)

subtypeSurvival(eset=OvarianEsets[[14]], sig=ss, weighted=FALSE, surv.type="os", time.cens=10 * 365, condensed=TRUE, resdir="saveres", nthread=nbcore, sig.method="weighted.average", sig.scaling="robust")

subtypeSurvival(eset=datall$merged, sig=ss, weighted=FALSE, surv.type="dmfs", time.cens=10 * 365, condensed=TRUE, resdir=saveres2, nthread=nbcore, sig.method="weighted.average", sig.scaling="robust")
