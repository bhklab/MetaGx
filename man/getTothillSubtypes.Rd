\name{getTothillSubtypes}
\alias{getTothillSubtypes}

\title{
Get Ovarian Cancer Subtypes (as defined by Tothill et al. (2008))
}
\description{
Classifies gene expression samples into subtypes using the classification scheme described by Tothill et al. (2008).
}
\usage{
getTothillSubtypes(eset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
An ExpressionSet object to be classified.
}
  \item{gene.mapping}{
The gene mapping of the given ExpressionSet, which will be used as features for the DLDA classifier: "Entrez.ID" (default) or "Probe.ID"
}
}
\details{
Trains an diagonal linear discriminant analysis (DLDA) classifer on the Tothill dataset. Class labels for the training data are from Supplementary Table 1, and Entrez IDs and probe sets reported in their Supplementary Table 2 of the original study. The DLDA is trained on the common Entrez IDs or probe sets between the list from the supplementary text, the train data, and the provided ExpressionSet.
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
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
