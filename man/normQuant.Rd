\name{normQuant}
\alias{normQuant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
x
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
normQuant(A, ties = TRUE, normvector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
numeric matrix. Missing values are allowed.
}
  \item{ties}{
If ‘TRUE’, ties in each row of ‘A’ are treated in careful way. tied values will be normalized to the mean of the corresponding pooled quantiles.
}
  \item{normvector}{
numeric vector of values corresponding to the quantiles distribution to fit. Note that 'normvector' must have the same length than the number of columns of 'A'
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A matrix of the same dimensions as ‘A’ containing the normalized values.
}
\references{
Bolstad, B. M., Irizarry R. A., Astrand, M., and Speed, T. P. (2003), A comparison of normalization methods for high density oligonucleotide array data based on bias and variance. \emph{Bioinformatics}, \bold{19}:185-193.
}
\author{
Benjamin Haibe-Kains
}
\note{
Function inpired from normalizaeQuantiles from the limma package (Gordon Smyth)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
dd <- matrix(runif(1500), ncol=15, nrow=10)
dd.nq <- normQuant(A=dd)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
