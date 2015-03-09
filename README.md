MetaGx2
=======

Scripts to perform meta-analysis of cancer gene expression datasets

Dependencies:

Install the R/Bioconductior dependencies:


pp <- c("Biobase", "BiocGenerics", "org.Hs.eg.db", "survival", "survcomp", "genefu", "mRMRe", "WriteXLS")

source("http://bioconductor.org/biocLite.R")

myrepos <- biocinstallRepos()

rr <- biocLite(pkgs=pp, dependencies=TRUE, type="source", destdir=".")


If using Windows:
- WriteXLS needs Perl to be installed. To check if perl has been installed:
    library(WriteXLS)
    testPerl(perl="perl", verbose=TRUE)

TODO
  - Adapt the package to handle curatedOvarianData
  - Weighted survival does not work properly (number of patients per time points should be "unweighted")
