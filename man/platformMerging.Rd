\name{platformMerging}
\alias{platformMerging}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
x
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
platformMerging(esets)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{esets}{
%%     ~~Describe \code{esets} here~~
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
function (esets) 
{
    if (class(esets) == "ExpressionSet") {
        return(esets)
    }
    if (is.list(esets)) {
        esets.check <- sapply(esets, function(x) {
            return(class(x) == "ExpressionSet")
        })
        if (any(!esets.check)) {
            stop("Some esets in the list are not ExpressionSet")
        }
        if (length(esets) == 1) {
            eset <- esets[[1]]
        }
        else {
            for (j in 2:length(esets)) {
                Biobase::exprs(eset) <- rbind(Biobase::exprs(eset), 
                  Biobase::exprs(esets[[j]]))
                Biobase::featureData(eset)@data <- rbind(Biobase::featureData(eset)@data, 
                  Biobase::featureData(esets[[j]])@data)
            }
            duplix <- duplicated(as.character(Biobase::featureData(eset)@data[, 
                "ENTREZID"]))
            Biobase::exprs(eset) <- Biobase::exprs(eset)[!duplix, 
                , drop = FALSE]
        }
        return(eset)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
