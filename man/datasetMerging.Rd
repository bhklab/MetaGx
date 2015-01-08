\name{datasetMerging}
\alias{datasetMerging}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
x
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
datasetMerging(esets, method = c("union", "intersect"), standardization = c("quantile", "robust.scaling", "scaling", "none"), nthread = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{esets}{
%%     ~~Describe \code{esets} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{standardization}{
%%     ~~Describe \code{standardization} here~~
}
  \item{nthread}{
%%     ~~Describe \code{nthread} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (esets, method = c("union", "intersect"), standardization = c("quantile", 
    "robust.scaling", "scaling", "none"), nthread = 1) 
{
    require(Biobase)
    require(genefu)
    if (nthread > 1) {
        require(parallel)
    }
    method <- match.arg(method)
    standardization <- match.arg(standardization)
    ugid <- lapply(esets, function(x) {
        return(Biobase::featureData(x)@data)
    })
    ugid <- do.call(rbind, ugid)
    ugid <- ugid[!is.na(ugid[, "ENTREZID"]) & !duplicated(as.character(ugid[, 
        "ENTREZID"])), , drop = FALSE]
    rownames(ugid) <- gsub(sprintf("(\%s).", paste(names(esets), 
        collapse = "|")), "", rownames(ugid))
    switch(method, union = {
        feature.merged <- ugid
    }, intersect = {
        feature.merged <- lapply(esets, function(x) {
            return(as.character(Biobase::featureData(x)@data[, 
                "ENTREZID"]))
        })
        feature.merged <- table(unlist(feature.merged))
        feature.merged <- names(feature.merged)[feature.merged == 
            length(esets)]
        feature.merged <- ugid[match(feature.merged, as.character(ugid[, 
            "ENTREZID"])), , drop = FALSE]
    }, {
        stop("Unknown method")
    })
    exprs.merged <- lapply(esets, function(x, y) {
        ee <- Biobase::exprs(x)
        eem <- matrix(NA, nrow = length(y), ncol = ncol(ee), 
            dimnames = list(y, colnames(ee)))
        eem[rownames(ee), colnames(ee)] <- ee
        return(eem)
    }, y = rownames(feature.merged))
    exprs.merged <- do.call(cbind, exprs.merged)
    ucid <- lapply(esets, function(x) {
        return(colnames(phenoData(x)@data))
    })
    ucid <- table(unlist(ucid))
    ucid <- names(ucid)[ucid == length(esets)]
    clinicinfo.merged <- lapply(esets, function(x, y) {
        ee <- Biobase::pData(x)[, y, drop = FALSE]
    }, y = ucid)
    clinicinfo.merged <- do.call(rbind, clinicinfo.merged)
    rownames(clinicinfo.merged) <- gsub(sprintf("(\%s).", paste(names(esets), 
        collapse = "|")), "", rownames(clinicinfo.merged))
    eset.merged <- ExpressionSet(assayData = exprs.merged, phenoData = AnnotatedDataFrame(data = clinicinfo.merged), 
        featureData = AnnotatedDataFrame(data = feature.merged))
    experimentData(eset.merged)@preprocessing <- list(normalization = "mixed", 
        package = "unspecified", version = "0")
    annotation(eset.merged) <- "mixed"
    switch(standardization, none = {
    }, quantile = {
        require(limma)
        require(genefu)
        ee <- exprs(eset.merged)
        splitix <- parallel::splitIndices(nx = ncol(ee), ncl = nthread)
        mcres <- parallel::mclapply(splitix, function(x, data) {
            res <- apply(data[, x, drop = FALSE], 2, function(dx) {
                return((genefu::rescale(dx, q = 0.05, na.rm = TRUE) - 
                  0.5) * 2)
            })
            return(res)
        }, data = ee, mc.cores = nthread)
        ee <- do.call(cbind, mcres)
        ee <- limma::normalizeBetweenArrays(object = ee, method = "quantile")
        exprs(eset.merged) <- ee
    }, robust.scling = {
        require(genefu)
        ee <- exprs(eset.merged)
        splitix <- parallel::splitIndices(nx = ncol(ee), ncl = nthread)
        mcres <- parallel::mclapply(splitix, function(x, data) {
            res <- apply(data[, x, drop = FALSE], 2, function(dx) {
                return((genefu::rescale(dx, q = 0.05, na.rm = TRUE) - 
                  0.5) * 2)
            })
            return(res)
        }, data = ee, mc.cores = nthread)
        ee <- do.call(cbind, mcres)
        exprs(eset.merged) <- ee
    }, scaling = {
        ee <- exprs(eset.merged)
        splitix <- parallel::splitIndices(nx = ncol(ee), ncl = nthread)
        mcres <- parallel::mclapply(splitix, function(x, data) {
            return(apply(data[, x, drop = FALSE], 2, scale))
        }, data = ee, mc.cores = nthread)
        ee <- do.call(cbind, mcres)
        exprs(eset.merged) <- ee
    }, {
        stop("Unknown data standardization method")
    })
    return(eset.merged)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
