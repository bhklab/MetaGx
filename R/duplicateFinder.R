########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

         
`duplicateFinder` <- 
function (eset, topvar.genes=1000, dupl.cor=0.99, method=c("pearson", "spearman", "kendall"), cor.matrix=FALSE, nthread=1) {
  ## find duplicates based on correlation of gene expresison profiles
  #
  # Arga:
  #   eset: an expressionSet object
  #   topvar.genes: number of most variant genes used to define the expression profiles
  #
  # Returns
  #   list of duplicates sample names
  
  method <- match.arg(method)
  if (topvar.genes < 3) { topvar.genes <- length(featureNames(eset)) }
  ## select the most variant genes
  ## at least in 80% of the datasets
  iix <- apply(exprs(eset), 1, function (x, y) {
    return ((sum(is.na(x)) / length(x)) < ( 1 - y))
  }, y=0.8)
  varg <- Biobase::featureNames(eset)[iix][order(apply(exprs(eset)[iix, , drop=FALSE], 1, var, na.rm=TRUE), decreasing=TRUE)[1:topvar.genes]]
  
  ## alternative, inefficient approach
  # pairs <- t(combn(1:length(Biobase::sampleNames(eset)), 2, simplify=TRUE))
  # splitix <- parallel::splitIndices(nx=nrow(pairs), ncl=nthread)
  # splitix <- splitix[sapply(splitix, length) > 0]
  # mcres <- parallel::mclapply(splitix, function(x, pairs, expr, method) {
  #   res <- apply(pairs[x, , drop=FALSE], 1, function (x, expr, method) {
  #     return (cor(x=expr[ , x[1], drop=FALSE], y=expr[ , x[2], drop=FALSE], method=method, use="complete.obs"))
  #   }, expr=expr, method=method)
  #   return (res)
  # }, pairs=pairs, expr=Biobase::exprs(eset)[varg, , drop=FALSE], method=method)
  # res <- t(do.call(cbind, mcres))
  # res <- unlist(apply(res, 2, function (x, y, gid) {
  #   rr <- matrix(NA, nrow=length(gid), ncol=length(gid), dimnames=list(gid, gid))
  #   rr[y] <- x
  #   rr[y[ , 2:1]] <- x
  #   diag(rr) <- 1
  #   return (list(rr))
  # }, y=pairs, gid), recursive=FALSE)
  
  ## more efficient but still slow approach
  # splitix <- parallel::splitIndices(nx=length(Biobase::sampleNames(eset)), ncl=nthread)
  # splitix <- splitix[sapply(splitix, length) > 0]
  # mcres <- parallel::mclapply(splitix, function(splitix, expr, method) {
  #     cores <- cor(x=expr[ , splitix, drop=FALSE], y=expr, method=method, use="pairwise.complete.obs")
  #   }, expr=Biobase::exprs(eset)[varg, , drop=FALSE], method=method)
  # cor.samples <- do.call(rbind, mcres)
  
  ## using mRMRe
  nn <- mRMRe::get.thread.count()
  mRMRe::set.thread.count(nthread)
  expr <- mRMRe::mRMR.data(data=data.frame(Biobase::exprs(eset)[varg, , drop=FALSE]))
  cor.samples <- mRMRe::mim(object=expr, continuous_estimator=method, method="cor")
  mRMRe::set.thread.count(nn)
  
  if (cor.matrix) { return (cor.samples) }
  diag(cor.samples) <- NA
  ## create list of duplicates for each sample
  duplix <- apply(cor.samples, 1, function (x, y) {
    res <- names(x)[!is.na(x) & x > y]
    return (res)
  }, y=dupl.cor)
  duplix <- duplix[sapply(duplix, length) > 0]
  return (duplix)
}

