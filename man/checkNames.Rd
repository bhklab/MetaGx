\name{checkNames}
\alias{checkNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
x
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
checkNames(eset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
%%     ~~Describe \code{eset} here~~
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
function (eset) 
{
    require(Biobase)
    check.feature <- intersect(rownames(Biobase::exprs(eset)), 
        rownames(Biobase::fData(eset)))
    if (length(check.feature) == 0) {
        warning("Names of features do not match between expressions and annotations")
        return(NULL)
    }
    else {
        if (length(check.feature) != nrow(Biobase::exprs(eset)) || 
            length(check.feature) != nrow(Biobase::fData(eset))) {
            warning("Some features are missing between expressions and annotations")
        }
    }
    check.sample <- intersect(colnames(Biobase::exprs(eset)), 
        rownames(Biobase::pData(eset)))
    if (length(check.sample) == 0) {
        warning("Names of samples do not match between expressions and phenotypes")
        return(NULL)
    }
    else {
        if (length(check.sample) != ncol(Biobase::exprs(eset)) || 
            length(check.sample) != nrow(Biobase::pData(eset))) {
            warning("Some samples are missing between expressions and phenotypes")
        }
    }
    Biobase::exprs(eset) <- Biobase::exprs(eset)[check.feature, 
        check.sample, drop = FALSE]
    Biobase::fData(eset) <- Biobase::fData(eset)[check.feature, 
        , drop = FALSE]
    Biobase::pData(eset) <- Biobase::pData(eset)[check.sample, 
        , drop = FALSE]
    Biobase::pData(eset)[, "samplename"] <- rownames(Biobase::pData(eset))
    return(eset)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
