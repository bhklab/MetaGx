\name{duplicateFinder}
\alias{duplicateFinder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
x
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
duplicateFinder(eset, var.genes = 1000, dupl.cor = 0.95, nthread = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
%%     ~~Describe \code{eset} here~~
}
  \item{var.genes}{
%%     ~~Describe \code{var.genes} here~~
}
  \item{dupl.cor}{
%%     ~~Describe \code{dupl.cor} here~~
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
function (eset, var.genes = 1000, dupl.cor = 0.95, nthread = 1) 
{
    require(Biobase)
    if (nthread > 1) {
        require(parallel)
    }
    iix <- apply(exprs(eset), 1, function(x, y) {
        return((sum(is.na(x))/length(x)) < (1 - y))
    }, y = 0.8)
    varg <- Biobase::featureNames(eset)[iix][order(apply(exprs(eset)[iix, 
        , drop = FALSE], 1, var, na.rm = TRUE), decreasing = TRUE)[1:var.genes]]
    splitix <- parallel::splitIndices(nx = length(sampleNames(eset)), 
        ncl = nthread)
    mcres <- parallel::mclapply(splitix, function(splitix, data) {
        cores <- cor(x = data[, splitix, drop = FALSE], y = data, 
            use = "complete.obs")
    }, data = exprs(eset)[varg, , drop = FALSE])
    cor.samples <- do.call(rbind, mcres)
    diag(cor.samples) <- NA
    duplix <- apply(cor.samples, 1, function(x, y) {
        res <- names(x)[!is.na(x) & x > y]
        return(res)
    }, y = dupl.cor)
    duplix <- duplix[sapply(duplix, length) > 0]
    return(duplix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
