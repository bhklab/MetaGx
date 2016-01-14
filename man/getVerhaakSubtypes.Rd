\name{getVerhaakSubtypes}
\alias{getVerhaakSubtypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getVerhaakSubtypes(eset)
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
    supplementary.data <- read.xls("../inst/extdata/JCI65833sd1.xls", 
        sheet = 7, skip = 1)
    genesets <- lapply(levels(supplementary.data$CLASS), function(y) as.character(supplementary.data[supplementary.data$CLASS == 
        y, 1]))
    names(genesets) <- levels(supplementary.data$CLASS)
    gsva.out <- gsva(exprs(eset), genesets, method = "ssgsea", 
        min.sz = 10, tau = 0.75, parallel.sz = 1)
    gsva.out <- t(gsva.out)
    subclasses <- as.factor(apply(gsva.out, 1, function(x) colnames(gsva.out)[which.max(x)]))
    pData(eset) <- data.frame(pData(eset), Verhaak.subtypes = subclasses)
    return(list(Annotated.eset = eset, gsva.out = gsva.out))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
