\name{probeGeneMapping}
\alias{probeGeneMapping}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
x
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
probeGeneMapping(eset, platform = c("MISC", "GPL8300", "GPL96", "GPL97", "GPL570", "GPL1352"), method = c("variance", "jetset"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
%%     ~~Describe \code{eset} here~~
}
  \item{platform}{
%%     ~~Describe \code{platform} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (eset, platform = c("MISC", "GPL8300", "GPL96", "GPL97", 
    "GPL570", "GPL1352"), method = c("variance", "jetset")) 
{
    require(org.Hs.eg.db)
    platform <- match.arg(platform)
    method <- match.arg(method)
    platf.map <- rbind(c("MISC", "variance", ""), c("GPL8300", 
        "jetset", "hgu95av2"), c("GPL96", "jetset", "hgu133a"), 
        c("GPL97", "jetset", "hgu133b"), c("GPL570", "jetset", 
            "hgu133plus2"), c("GPL1352", "jetset", "u133x3p"))
    dimnames(platf.map) <- list(platf.map[, 1], c("platform", 
        "method", "parameters"))
    if (!is.element(method, platf.map[platf.map[, "platform"], 
        "method"])) {
        stop(sprintf("Method \%s cannot be applied on platform \%s\nUse the following method(s) instead: \%s", 
            method, platform, paste(x = platf.map[platf.map[, 
                "platform"] == platform, "method"], collapse = ", ")))
    }
    params <- platf.map[which(platform == platf.map[, "platform"]), 
        "parameters"]
    Biobase::fData(eset) <- Biobase::fData(eset)[, c("ENTREZID", 
        "SYMBOL"), drop = FALSE]
    switch(method, jetset = {
        require(jetset.bhk)
        js <- jetset.bhk::jscores(chip = params, probeset = rownames(Biobase::exprs(eset)))
        js <- js[rownames(Biobase::exprs(eset)), , drop = FALSE]
        geneid1 <- as.character(js[, "EntrezID"])
        names(geneid1) <- rownames(js)
        geneid2 <- sort(unique(geneid1))
        names(geneid2) <- paste("geneid", geneid2, sep = ".")
        gix1 <- !is.na(geneid1)
        gix2 <- !is.na(geneid2)
        geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
        gg <- names(geneid1)[is.element(geneid1, geneid.common)]
        gid <- geneid1[is.element(geneid1, geneid.common)]
        gid.dupl <- unique(gid[duplicated(gid)])
        gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
        gid.uniq <- gid[!is.element(gid, gid.dupl)]
        gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
        js <- data.frame(js, best = FALSE)
        js[gg.uniq, "best"] <- TRUE
        if (length(gid.dupl) > 0) {
            myscore <- js[gg.dupl, "overall"]
            myscore <- cbind(probe = gg.dupl, gid = geneid1[gg.dupl], 
                score = myscore)
            myscore <- myscore[order(as.numeric(myscore[, "score"]), 
                decreasing = TRUE, na.last = TRUE), , drop = FALSE]
            myscore <- myscore[!duplicated(myscore[, "gid"]), 
                , drop = FALSE]
            js[myscore[, "probe"], "best"] <- TRUE
        }
        probes <- rownames(Biobase::exprs(eset))[js[, "best"]]
        names(probes) <- paste("geneid", js[js[, "best"], "EntrezID"], 
            sep = ".")
        gid <- js[js[, "best"], "EntrezID"]
        gsymb <- js[js[, "best"], "symbol"]
        Biobase::exprs(eset) <- Biobase::exprs(eset)[probes, 
            , drop = FALSE]
        rownames(Biobase::exprs(eset)) <- names(probes)
        Biobase::featureData(eset)@data <- Biobase::featureData(eset)@data[probes, 
            , drop = FALSE]
        rownames(Biobase::featureData(eset)@data) <- names(probes)
        Biobase::featureData(eset)@data[, "ENTREZID"] <- gid
        Biobase::featureData(eset)@data[, "SYMBOL"] <- gsymb
    }, variance = {
        gid <- as.character(Biobase::featureData(eset)@data[, 
            "ENTREZID"])
        names(gid) <- rownames(Biobase::exprs(eset))
        ugid <- sort(unique(gid))
        rr <- genefu::geneid.map(geneid1 = gid, data1 = t(Biobase::exprs(eset)), 
            geneid2 = ugid)
        probes <- colnames(rr$data1)
        names(probes) <- paste("geneid", rr$geneid1, sep = ".")
        Biobase::exprs(eset) <- Biobase::exprs(eset)[probes, 
            , drop = FALSE]
        rownames(Biobase::exprs(eset)) <- names(probes)
        Biobase::featureData(eset)@data <- Biobase::featureData(eset)@data[probes, 
            , drop = FALSE]
        rownames(Biobase::featureData(eset)@data) <- names(probes)
        gs <- toTable(org.Hs.egSYMBOL)
        gs <- gs[!duplicated(gs[, "gene_id"]), , drop = FALSE]
        rownames(gs) <- gs[, "gene_id"]
        gs <- gs[as.character(Biobase::featureData(eset)@data[, 
            "ENTREZID"]), "symbol"]
        Biobase::featureData(eset)@data[, "SYMBOL"] <- as.character(gs)
    }, {
        stop(sprintf("Unknow method for probe-gene mapping for platform \%s", 
            platform))
    })
    return(eset)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
