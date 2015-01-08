MetaGx
======

Scripts to download and curate breast cancer datasets from InSilicoDB

Dependencies:

Install the R/Bioconductior dependencies:


pp <- c("InSilicoDb", "Biobase", "BiocGenerics", "org.Hs.eg.db", "survival", "survcomp", "genefu", "mRMRe", "WriteXLS")

source("http://bioconductor.org/biocLite.R")

myrepos <- biocinstallRepos()

rr <- biocLite(pkgs=pp, dependencies=TRUE, type="source", destdir=".")


TODO
- Weighted survival does not work properly (number of patients per time points should be "unweighted")
- METABRIC and STAT1?
