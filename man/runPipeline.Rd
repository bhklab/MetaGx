\name{runPipeline}
\alias{runPipeline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
x
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runPipeline(sbt.model = c("scmgene", "scmod2", "scmod1", "pam50", "ssp2006", "ssp2003"), resdir = "cache", probegene.method, remove.duplicates = TRUE, topvar.genes = 1000, duplicates.cor = 0.975, datasets, nthread = 1, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sbt.model}{
%%     ~~Describe \code{sbt.model} here~~
}
  \item{resdir}{
%%     ~~Describe \code{resdir} here~~
}
  \item{probegene.method}{
%%     ~~Describe \code{probegene.method} here~~
}
  \item{remove.duplicates}{
%%     ~~Describe \code{remove.duplicates} here~~
}
  \item{topvar.genes}{
%%     ~~Describe \code{topvar.genes} here~~
}
  \item{duplicates.cor}{
%%     ~~Describe \code{duplicates.cor} here~~
}
  \item{datasets}{
%%     ~~Describe \code{datasets} here~~
}
  \item{nthread}{
%%     ~~Describe \code{nthread} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (sbt.model = c("scmgene", "scmod2", "scmod1", "pam50", 
    "ssp2006", "ssp2003"), resdir = "cache", probegene.method, 
    remove.duplicates = TRUE, topvar.genes = 1000, duplicates.cor = 0.975, 
    datasets, nthread = 1, verbose = TRUE) 
{
    badchars <- "[\265]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[\%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
    if (!file.exists(resdir)) {
        dir.create(resdir, showWarnings = FALSE, recursive = TRUE)
    }
    if (!file.exists(file.path(resdir, "processed"))) {
        dir.create(file.path(resdir, "processed"), showWarnings = FALSE, 
            recursive = TRUE)
    }
    availcore <- parallel::detectCores()
    if (nthread > availcore) {
        nthread <- availcore
    }
    options(mc.cores = nthread)
    if (missing(probegene.method)) {
        probegene.method <- rbind(c("MISC", "variance"), c("GPL8300", 
            "jetset"), c("GPL96", "jetset"), c("GPL97", "jetset"), 
            c("GPL570", "jetset"), c("GPL1352", "jetset"))
        colnames(probegene.method) <- c("platform", "method")
    }
    sbt.model <- match.arg(sbt.model)
    if (missing(datasets)) {
        datasets <- read.csv(system.file(file.path("extdata", 
            "datasets.csv"), package = "bcmeta"), stringsAsFactors = FALSE)
    }
    datasets <- datasets[datasets[, "Include"], , drop = FALSE]
    clin.info <- c("samplename", "id", "series", "dataset", "age", 
        "node", "er", "e.rfs", "e.os", "e.dmfs", "grade", "size", 
        "pgr", "her2", "t.rfs", "t.os", "t.dmfs", "treatment", 
        "tissue")
    InSilicoLogin(login = "bhaibeka@gmail.com", password = "747779bec8a754b91076d6cc1f700831")
    cinfo <- clin.info
    eset.all <- NULL
    for (i in 1:nrow(datasets)) {
        ddn <- as.character(datasets[i, "Dataset.ID"])
        if (verbose) {
            message(sprintf("Get dataset \%s", ddn))
        }
        dataset.fn <- file.path(resdir, "processed", sprintf("\%s_processed.RData", 
            ddn))
        if (!file.exists(dataset.fn)) {
            platf <- InSilicoDb::getPlatforms(dataset = ddn)
            esets <- InSilicoDb::getDatasets(dataset = ddn, 
                norm = as.character(datasets[i, "Normalization"]), 
                curation = datasets[i, "Curation.ID"], features = "PROBE")
            if (is.null(unlist(esets))) {
                stop("Rerun the script when data are ready to download from InSilicoDB")
            }
            platf <- unlist(platf[!sapply(esets, is.null)])
            esets <- esets[!sapply(esets, is.null)]
            if (verbose) {
                message("\tCheck feature and sample names")
            }
            esets <- lapply(esets, checkNames)
            if (verbose) {
                message("\tProbe-gene mapping")
            }
            platf2 <- platf
            platf2[!is.element(platf, probegene.method[, "platform"])] <- "MISC"
            esets <- mapply(probeGeneMapping, eset = esets, platform = platf2, 
                method = probegene.method[match(platf2, probegene.method[, 
                  "platform"]), "method"])
            eset <- platformMerging(esets)
            colnames(Biobase::phenoData(eset)@data) <- gsub(" ", 
                ".", colnames(Biobase::phenoData(eset)@data))
            cinfo <- intersect(cinfo, colnames(Biobase::phenoData(eset)@data))
            Biobase::phenoData(eset)@data <- data.frame(apply(Biobase::phenoData(eset)@data, 
                2, function(x) {
                  if (is.factor(x)) {
                    x <- as.character(x)
                  }
                  return(x)
                }), stringsAsFactors = FALSE)
            Biobase::phenoData(eset)@data[Biobase::phenoData(eset)@data == 
                "NA"] <- NA
            save(list = c("eset"), compress = TRUE, file = dataset.fn)
            if (verbose) {
                message("")
            }
        }
        else {
            load(file = dataset.fn)
        }
        switch(datasets[i, "Dataset.ID"], GSE2034 = {
            if ("GSE5327" \%in\% names(eset.all)) {
                if (verbose) {
                  message("\tMerge GSE2034 and GSE5327")
                }
                cg <- union(rownames(Biobase::exprs(eset.all[["GSE5327"]])), 
                  rownames(Biobase::exprs(eset)))
                cs <- c(colnames(Biobase::exprs(eset.all[["GSE5327"]])), 
                  colnames(Biobase::exprs(eset)))
                ee <- matrix(NA, nrow = length(cg), ncol = length(cs), 
                  dimnames = list(cg, cs))
                ee[rownames(Biobase::exprs(eset.all[["GSE5327"]])), 
                  colnames(Biobase::exprs(eset.all[["GSE5327"]]))] <- Biobase::exprs(eset.all[["GSE5327"]])
                ee[rownames(Biobase::exprs(eset)), colnames(Biobase::exprs(eset))] <- Biobase::exprs(eset)
                Biobase::exprs(eset.all[["GSE5327"]]) <- ee
                cf <- c(colnames(Biobase::featureData(eset.all[["GSE5327"]])@data), 
                  colnames(Biobase::featureData(eset)@data))
                ee <- data.frame(matrix(NA, nrow = length(cg), 
                  ncol = length(cf), dimnames = list(cg, cf)))
                ee[rownames(Biobase::featureData(eset.all[["GSE5327"]])@data), 
                  colnames(Biobase::featureData(eset.all[["GSE5327"]])@data)] <- Biobase::featureData(eset.all[["GSE5327"]])@data
                ee[rownames(Biobase::featureData(eset)@data), 
                  colnames(Biobase::featureData(eset)@data)] <- Biobase::featureData(eset)@data
                Biobase::featureData(eset.all[["GSE5327"]]) <- ee
                cp <- intersect(colnames(Biobase::phenoData(eset.all[["GSE5327"]])@data), 
                  colnames(Biobase::phenoData(eset)@data))
                ee <- data.frame(matrix(NA, nrow = length(cs), 
                  ncol = length(cp), dimnames = list(cs, cp)))
                ee[rownames(Biobase::phenoData(eset.all[["GSE5327"]])@data), 
                  cp] <- Biobase::phenoData(eset.all[["GSE5327"]])@data[, 
                  cp]
                ee[rownames(Biobase::phenoData(eset)@data), cp] <- Biobase::phenoData(eset)@data[, 
                  cp]
                Biobase::phenoData(eset.all[["GSE5327"]]) <- ee
                names(eset.all)[names(eset.all) == "GSE2034"] <- "GSE2034.GSE5327"
            } else {
                eset.all <- c(eset.all, eset)
                names(eset.all)[length(eset.all)] <- datasets[i, 
                  "Dataset.ID"]
            }
        }, GSE5327 = {
            if ("GSE2034" \%in\% names(eset.all)) {
                if (verbose) {
                  message("\tMerge GSE2034 and GSE5327")
                }
                cg <- union(rownames(Biobase::exprs(eset.all[["GSE2034"]])), 
                  rownames(Biobase::exprs(eset)))
                cs <- c(colnames(Biobase::exprs(eset.all[["GSE2034"]])), 
                  colnames(Biobase::exprs(eset)))
                ee <- matrix(NA, nrow = length(cg), ncol = length(cs), 
                  dimnames = list(cg, cs))
                ee[rownames(Biobase::exprs(eset)), colnames(Biobase::exprs(eset))] <- Biobase::exprs(eset)
                ee[rownames(Biobase::exprs(eset.all[["GSE2034"]])), 
                  colnames(Biobase::exprs(eset.all[["GSE2034"]]))] <- Biobase::exprs(eset.all[["GSE2034"]])
                Biobase::exprs(eset.all[["GSE2034"]]) <- ee
                cf <- intersect(colnames(Biobase::featureData(eset.all[["GSE2034"]])@data), 
                  colnames(Biobase::featureData(eset)@data))
                ee <- data.frame(matrix(NA, nrow = length(cg), 
                  ncol = length(cf), dimnames = list(cg, cf)))
                ee[rownames(Biobase::featureData(eset)@data), 
                  colnames(Biobase::featureData(eset)@data)] <- Biobase::featureData(eset)@data
                ee[rownames(Biobase::featureData(eset.all[["GSE2034"]])@data), 
                  colnames(Biobase::featureData(eset.all[["GSE2034"]])@data)] <- Biobase::featureData(eset.all[["GSE2034"]])@data
                Biobase::featureData(eset.all[["GSE2034"]])@data <- ee
                cp <- intersect(colnames(Biobase::phenoData(eset.all[["GSE2034"]])@data), 
                  colnames(Biobase::phenoData(eset)@data))
                ee <- data.frame(matrix(NA, nrow = length(cs), 
                  ncol = length(cp), dimnames = list(cs, cp)))
                ee[rownames(Biobase::phenoData(eset.all[["GSE2034"]])@data), 
                  cp] <- Biobase::phenoData(eset.all[["GSE2034"]])@data[, 
                  cp]
                ee[rownames(Biobase::phenoData(eset)@data), cp] <- Biobase::phenoData(eset)@data[, 
                  cp]
                Biobase::phenoData(eset.all[["GSE2034"]])@data <- ee
                names(eset.all)[names(eset.all) == "GSE2034"] <- "GSE2034.GSE5327"
            } else {
                eset.all <- c(eset.all, eset)
                names(eset.all)[length(eset.all)] <- datasets[i, 
                  "Dataset.ID"]
            }
        }, GSE2109 = {
            iix <- Biobase::phenoData(eset)@data[, "Anatomical.site"] == 
                "breast"
            Biobase::exprs(eset) <- Biobase::exprs(eset)[, iix, 
                drop = FALSE]
            Biobase::phenoData(eset)@data <- Biobase::phenoData(eset)@data[iix, 
                , drop = FALSE]
            eset.all <- c(eset.all, eset)
            names(eset.all)[length(eset.all)] <- datasets[i, 
                "Dataset.ID"]
        }, {
            eset.all <- c(eset.all, eset)
            names(eset.all)[length(eset.all)] <- datasets[i, 
                "Dataset.ID"]
        })
    }
    if (verbose) {
        message("Update clinical information")
    }
    for (j in 1:length(eset.all)) {
        Biobase::phenoData(eset.all[[j]])@data <- Biobase::phenoData(eset.all[[j]])@data[, 
            cinfo, drop = FALSE]
        tt <- intersect(colnames(Biobase::phenoData(eset.all[[j]])@data), 
            c("age", "size", "er", "her2", "pgr", "grade", c(paste(c("t", 
                "e"), rep(c("rfs", "dmfs", "os"), each = 2), 
                sep = "."))))
        Biobase::phenoData(eset.all[[j]])@data[, tt] <- data.frame(apply(Biobase::phenoData(eset.all[[j]])@data[, 
            tt, drop = FALSE], 2, as.numeric), stringsAsFactors = FALSE)
        surv.time <- Biobase::phenoData(eset.all[[j]])@data[, 
            "t.rfs"]
        surv.time[is.na(Biobase::phenoData(eset.all[[j]])@data[, 
            "t.rfs"])] <- Biobase::phenoData(eset.all[[j]])@data[is.na(Biobase::phenoData(eset.all[[j]])@data[, 
            "t.rfs"]), "t.dmfs"]
        surv.event <- Biobase::phenoData(eset.all[[j]])@data[, 
            "e.rfs"]
        surv.event[is.na(Biobase::phenoData(eset.all[[j]])@data[, 
            "e.rfs"])] <- Biobase::phenoData(eset.all[[j]])@data[is.na(Biobase::phenoData(eset.all[[j]])@data[, 
            "e.rfs"]), "e.dmfs"]
        Biobase::phenoData(eset.all[[j]])@data <- cbind(Biobase::phenoData(eset.all[[j]])@data, 
            t.dfs = surv.time, e.dfs = surv.event)
    }
    if (verbose) {
        message("Subtype classification")
    }
    for (j in 1:length(eset.all)) {
        eset.all[[j]] <- subtypeClassification(eset = eset.all[[j]], 
            model = sbt.model)
    }
    if (verbose) {
        message("Merge datasets")
    }
    eset.merged <- datasetMerging(esets = eset.all, nthread = nthread)
    duplicates <- duplicateFinder(eset = eset.merged, var.genes = topvar.genes, 
        dupl.cor = duplicates.cor)
    tt <- sapply(duplicates, paste, collapse = "///")
    Biobase::pData(eset.merged) <- cbind(Biobase::pData(eset.merged), 
        duplicates = NA)
    Biobase::pData(eset.merged)[names(tt), "duplicates"] <- tt
    eset.all <- lapply(eset.all, function(x, y) {
        Biobase::pData(x) <- cbind(Biobase::pData(x), duplicates = NA)
        nn <- intersect(rownames(pData(x)), y)
        Biobase::pData(x)[nn, "duplicates"] <- y[nn]
        return(x)
    }, y = tt)
    if (remove.duplicates) {
        rmix <- duplicates
        ii <- 1
        while (length(rmix) > ii) {
            rmix <- rmix[!is.element(names(rmix), rmix[[ii]])]
            ii <- ii + 1
        }
        rmix <- unique(unlist(rmix))
        keepix <- setdiff(sampleNames(eset.merged), rmix)
        Biobase::exprs(eset.merged) <- Biobase::exprs(eset.merged)[, 
            keepix, drop = FALSE]
        Biobase::pData(eset.merged) <- Biobase::pData(eset.merged)[keepix, 
            , drop = FALSE]
        eset.all <- lapply(eset.all, function(x, y) {
            keepix <- setdiff(sampleNames(x), y)
            Biobase::exprs(x) <- Biobase::exprs(x)[, keepix, 
                drop = FALSE]
            Biobase::pData(x) <- Biobase::pData(x)[keepix, , 
                drop = FALSE]
            return(x)
        }, y = rmix)
    }
    InSilicoLogout()
    return(list(merged = eset.merged, each = eset.all))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
