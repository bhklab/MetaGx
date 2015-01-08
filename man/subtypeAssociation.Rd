\name{subtypeAssociation}
\alias{subtypeAssociation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
x
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
subtypeAssociation(eset, geneid, boxp = TRUE, subtype.col, resdir, nthread = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
%%     ~~Describe \code{eset} here~~
}
  \item{geneid}{
%%     ~~Describe \code{geneid} here~~
}
  \item{boxp}{
%%     ~~Describe \code{boxp} here~~
}
  \item{subtype.col}{
%%     ~~Describe \code{subtype.col} here~~
}
  \item{resdir}{
%%     ~~Describe \code{resdir} here~~
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
function (eset, geneid, boxp = TRUE, subtype.col, resdir, nthread = 1) 
{
    if (class(eset) != "ExpressionSet") {
        stop("Handling list of expressionSet objects is not implemented yet")
    }
    if (missing(geneid)) {
        gened <- as.character(Biobase::fData(eset)[, "ENTREZID"])
    }
    sbts <- Biobase::pData(eset)[, "subtype"]
    if (sum(table(sbts) > 3) < 2) {
        warning("Not enough tumors in each subtype")
        return(NULL)
    }
    sbtu <- levels(sbts)
    if (missing(subtype.col)) {
        subtype.col <- rainbow(length(sbtu), alpha = 0.6)
    }
    else {
        if (length(subtype.col) < length(sbtu)) {
            stop(sprintf("Not enough color for \%i subtypes", 
                length(sbtu)))
        }
    }
    gid <- intersect(geneid, as.character(Biobase::fData(eset)[, 
        "ENTREZID"]))
    if (length(gid) == 0) {
        stop("Genes not in the expressionSet object")
    }
    if (length(gid) < length(geneid)) {
        warning(sprintf("\%i/\%i genes were present in the expressionSet object", 
            length(gid), length(geneid)))
    }
    splitix <- parallel::splitIndices(nx = length(gid), ncl = nthread)
    splitix <- splitix[sapply(splitix, length) > 0]
    mcres <- parallel::mclapply(splitix, function(x, ...) {
        pp <- lapply(gid[x], function(gid, eset, sbts, boxp, 
            resdir) {
            gsymb <- Biobase::fData(eset)[match(gid, Biobase::fData(eset)[, 
                "ENTREZID"]), "SYMBOL"]
            xx <- Biobase::exprs(eset)[paste("geneid", gid, sep = "."), 
                ]
            kt <- kruskal.test(x = xx, g = sbts)$p.value
            wt <- matrix(NA, nrow = length(sbtu), ncol = length(sbtu), 
                dimnames = list(sbtu, sbtu))
            wt1 <- pairwise.wilcox.test(x = xx, g = sbts, p.adjust.method = "none", 
                paired = FALSE, alternative = "greater")$p.value
            wt2 <- pairwise.wilcox.test(x = xx, g = sbts, p.adjust.method = "none", 
                paired = FALSE, alternative = "less")$p.value
            nix <- !is.na(wt1)
            wt[rownames(wt1), colnames(wt1)][nix] <- wt1[nix]
            nix <- !is.na(t(wt2))
            wt[colnames(wt2), rownames(wt2)][nix] <- t(wt2)[nix]
            diag(wt) <- 1
            if (boxp) {
                pdf(file.path(resdir, sprintf("subtype_association_boxplot_\%s.pdf", 
                  gsymb)))
                par(las = 2, mar = c(5, 4, 4, 2) + 0.1, xaxt = "n")
                mp <- boxplot(xx ~ sbts, las = 3, outline = FALSE, 
                  ylim = c(-2, 2), main = sprintf("\%s", gsymb), 
                  col = subtype.col)
                axis(1, at = 1:length(mp$names), tick = TRUE, 
                  labels = T)
                text(x = 1:length(mp$names), y = par("usr")[3] - 
                  (par("usr")[4] * 0.05), pos = 2, labels = mp$names, 
                  srt = 45, xpd = NA, font = 2, col = c("black"))
                text(x = 1:length(mp$names), y = par("usr")[3], 
                  pos = 3, labels = sprintf("n=\%i", table(sbts)[mp$names]), 
                  col = c("black"))
                legend("topleft", legend = sprintf("Kruskal-Wallis p-value = \%.1E", 
                  kt), bty = "n")
                dev.off()
            }
            return(list(kruskal.pvalue = kt, wilcoxon.pvalue = wt))
        }, eset = eset, sbts = sbts, boxp = boxp, resdir = resdir)
    }, gid = gid, eset = eset, sbts = sbts, boxp = boxp, resdir = resdir)
    pp <- do.call(c, mcres)
    names(pp) <- Biobase::fData(eset)[match(gid, Biobase::fData(eset)[, 
        "ENTREZID"]), "SYMBOL"]
    dd <- sapply(pp, function(x) {
        return(x[[1]])
    })
    dd <- data.frame(Kruskal.Wallis.pvalue = dd, Kruskal.Wallis.fdr = p.adjust(dd, 
        method = "fdr"), Biobase::fData(eset)[match(gid, Biobase::fData(eset)[, 
        "ENTREZID"]), ])
    write.csv(dd, file = file.path(resdir, "subtype_association_kruskal.csv"))
    mapply(function(x, y, resdir) {
        write.csv(x, file = file.path(resdir, sprintf("subtype_association_wilcoxon_\%s.csv", 
            y)))
    }, x = lapply(pp, function(x) {
        return(x[[2]])
    }), y = names(pp), resdir = resdir)
    return(pp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
