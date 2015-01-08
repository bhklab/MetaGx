\name{subtypeClassification}
\alias{subtypeClassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
x
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
subtypeClassification(eset, model = c("scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
%%     ~~Describe \code{eset} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
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
function (eset, model = c("scmgene", "scmod1", "scmod2", "pam50", 
    "ssp2006", "ssp2003")) 
{
    require(genefu)
    require(Biobase)
    model <- match.arg(model)
    sbt.conv <- rbind(c("ER+/HER2- Low Prolif", "LumA"), c("ER+/HER2- High Prolif", 
        "LumB"), c("HER2+", "Her2"), c("ER-/HER2-", "Basal"))
    colnames(sbt.conv) <- c("SCM.nomenclature", "SSP.nomenclature")
    sbtn2 <- c("LumA", "LumB", "Her2", "Basal", "Normal")
    datage <- t(Biobase::exprs(eset))
    annotge <- cbind(probe = rownames(Biobase::featureData(eset)@data), 
        EntrezGene.ID = as.character(Biobase::featureData(eset)@data[, 
            "ENTREZID"]))
    rownames(annotge) <- as.character(annotge[, "probe"])
    switch(model, scmgene = {
        sbts <- genefu::subtype.cluster.predict(sbt.model = scmgene.robust, 
            data = datage, annot = annotge, do.mapping = TRUE)[c("subtype2", 
            "subtype.proba2")]
        names(sbts) <- c("subtype", "subtype.proba")
        ss <- factor(x = sbts$subtype)
        levels(ss)[match(sbt.conv[!is.na(sbt.conv[, "SCM.nomenclature"]), 
            "SCM.nomenclature"], levels(ss))] <- sbt.conv[!is.na(sbt.conv[, 
            "SCM.nomenclature"]), "SSP.nomenclature"]
        sbts$subtype <- as.character(ss)
        iix <- match(sbt.conv[!is.na(sbt.conv[, "SCM.nomenclature"]), 
            "SCM.nomenclature"], colnames(sbts$subtype.proba))
        colnames(sbts$subtype.proba)[iix] <- sbt.conv[!is.na(sbt.conv[, 
            "SCM.nomenclature"]), "SSP.nomenclature"]
    }, scmod1 = {
        sbts <- genefu::subtype.cluster.predict(sbt.model = scmod1.robust, 
            data = datage, annot = annotge, do.mapping = TRUE)[c("subtype2", 
            "subtype.proba2")]
        names(sbts) <- c("subtype", "subtype.proba")
        ss <- factor(x = sbts$subtype)
        levels(ss)[match(sbt.conv[!is.na(sbt.conv[, "SCM.nomenclature"]), 
            "SCM.nomenclature"], levels(ss))] <- sbt.conv[!is.na(sbt.conv[, 
            "SCM.nomenclature"]), "SSP.nomenclature"]
        sbts$subtype <- as.character(ss)
        iix <- match(sbt.conv[!is.na(sbt.conv[, "SCM.nomenclature"]), 
            "SCM.nomenclature"], colnames(sbts$subtype.proba))
        colnames(sbts$subtype.proba)[iix] <- sbt.conv[!is.na(sbt.conv[, 
            "SCM.nomenclature"]), "SSP.nomenclature"]
    }, scmod2 = {
        sbts <- genefu::subtype.cluster.predict(sbt.model = scmod2.robust, 
            data = datage, annot = annotge, do.mapping = TRUE)[c("subtype2", 
            "subtype.proba2")]
        names(sbts) <- c("subtype", "subtype.proba")
        ss <- factor(x = sbts$subtype)
        levels(ss)[match(sbt.conv[!is.na(sbt.conv[, "SCM.nomenclature"]), 
            "SCM.nomenclature"], levels(ss))] <- sbt.conv[!is.na(sbt.conv[, 
            "SCM.nomenclature"]), "SSP.nomenclature"]
        sbts$subtype <- as.character(ss)
        iix <- match(sbt.conv[!is.na(sbt.conv[, "SCM.nomenclature"]), 
            "SCM.nomenclature"], colnames(sbts$subtype.proba))
        colnames(sbts$subtype.proba)[iix] <- sbt.conv[!is.na(sbt.conv[, 
            "SCM.nomenclature"]), "SSP.nomenclature"]
    }, pam50 = {
        sbts <- genefu::intrinsic.cluster.predict(sbt.model = pam50.robust, 
            data = datage, annot = annotge, do.mapping = TRUE)[c("subtype", 
            "subtype.proba")]
    }, ssp2006 = {
        sbts <- genefu::intrinsic.cluster.predict(sbt.model = ssp2006.robust, 
            data = datage, annot = annotge, do.mapping = TRUE)[c("subtype", 
            "subtype.proba")]
    }, ssp2003 = {
        sbts <- genefu::intrinsic.cluster.predict(sbt.model = ssp2003.robust, 
            data = datage, annot = annotge, do.mapping = TRUE)[c("subtype", 
            "subtype.proba")]
    }, {
        stop("Unknown subtype classification model")
    })
    colnames(sbts$subtype.proba) <- paste("subtyproba", colnames(sbts$subtype.proba), 
        sep = ".")
    Biobase::phenoData(eset)@data <- cbind(Biobase::phenoData(eset)@data, 
        subtype = sbts$subtype, sbts$subtype.proba)
    return(eset)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
