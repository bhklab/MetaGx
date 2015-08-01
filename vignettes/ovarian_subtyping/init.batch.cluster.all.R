library(Biobase)
library(MetaGxOvarian)

source("../../inst/extdata/reproduce.results.patientselection.config")
rule.2 <- c("histological_type","^ser$")
rule.3 <- c("summarystage","^late$")
rule.4 <- c("summarygrade","^high$")
### use this line if you do not want to get rid of duplicates
rm(remove.duplicates)
rescale <- FALSE

source(system.file("extdata", "createEsetList.R", package="MetaGxOvarian"))

esets.not.rescaled <- esets
# remove any genes with NA values
esets.not.rescaled <- lapply(esets.not.rescaled, function(eset) eset[apply(exprs(eset), 1, function(x) all(!is.na(x))),])
# only keep esets with at least 10000 genes
esets.not.rescaled <- esets.not.rescaled[sapply(esets.not.rescaled, function(x) nrow(x) > 10000)]

save(esets.not.rescaled, file="esets.not.rescaled.RData")
