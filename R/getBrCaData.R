########################
## Benjamin Haibe-Kains
## All rights Reserved
## December 23, 2013
########################

`getBrCaData` <- 
function (resdir="cache", probegene.method, remove.duplicates=TRUE, topvar.genes=1000, duplicates.cor=0.98, datasets, sbt.model=c("scmgene", "scmod2", "scmod1", "pam50", "ssp2006", "ssp2003"), merging.method=c("union", "intersection"), merging.std=c("quantile", "robust.scaling", "scaling", "none"), nthread=1, verbose=TRUE) {  

  merging.method <- match.arg(merging.method)
  merging.std <- match.arg(merging.std)
  
  badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

  ## directory where all the analysis results will be stored
  if(!file.exists(resdir)) { dir.create(resdir, showWarnings=FALSE, recursive=TRUE) }
  if(!file.exists(file.path(resdir, "processed"))) { dir.create(file.path(resdir, "processed"), showWarnings=FALSE, recursive=TRUE) }

  ## number of cpu cores available for the analysis pipeline
  ## set to 'NULL' if all the available cores should be used
  availcore <- parallel::detectCores()
  if (nthread > availcore) { nthread <- availcore }
  options("mc.cores"=nthread)

  ## probe gene mapping method for each platform
  if (missing(probegene.method)) {
    probegene.method <- rbind(
      c("MISC", "variance"),
      c("GPL8300", "jetset"),
      c("GPL96", "jetset"),
      c("GPL3921", "jetset"),
      c("GPL97", "jetset"),
      c("GPL570", "jetset"),
      c("GPL1352", "jetset"))
    colnames(probegene.method) <- c("platform", "method")
  }

  ## which subtype classification model
  ## scmgene, scmod1, scmod2, pam50, ssp2006, ssp2003
  sbt.model <- match.arg(sbt.model)

  ## number of top variant genes to use when comparing gene expression profiles

  ## maximum correlation between gene expression profiles of different samples

  ## read info about datasets
  if (missing(datasets)) {
    datasets <- read.csv(system.file(file.path("extdata", "datasets.csv"), package="MetaGx"), stringsAsFactors=FALSE)
  }
  datasets <- datasets[datasets[ , "Include"], , drop=FALSE]

  # clinical info
  clin.info <- c("samplename", "id", "series", "dataset", "age", "node", "er", "e.rfs", "e.os", "e.dmfs", "grade", "size", "pgr", "her2", "t.rfs", "t.os", "t.dmfs", "treatment", "tissue")

  ## log in inSilicoDb
  inSilicoDb::InSilicoLogin(login="bhaibeka@gmail.com", password="747779bec8a754b91076d6cc1f700831")

  ## clinical information
  cinfo <- clin.info

  ## download all datasets
  eset.all <- NULL
  for (i in 1:nrow(datasets)) {
    ddn <- stripWhiteSpace(as.character(datasets[i, "Dataset.ID"]))
    if (verbose) {
      message(sprintf("Get dataset %s", ddn))
    }
    dataset.fn <- file.path(resdir, "processed", sprintf("%s_processed.RData", ddn))
    if (!file.exists(dataset.fn)) {
      ## get dataset
      # inSilicoDb::getCurationInfo(dataset=as.character(datasets[i, "Dataset.ID"]))
      platf <- inSilicoDb::getPlatforms(dataset=ddn)
      esets <- inSilicoDb::getDatasets(dataset=ddn, norm=stripWhiteSpace(as.character(datasets[i, "Normalization"])), curation=datasets[i, "Curation.ID"], features="PROBE")
      if (is.null(unlist(esets))) {
        stop("Rerun the script when data are ready to download from inSilicoDb")
      }
      platf <- unlist(platf[!sapply(esets, is.null)])
      esets <- esets[!sapply(esets, is.null)]
      ## check sample and feature names
      if (verbose) {
        message("\tCheck feature and sample names")
      }
      esets <- lapply(esets, checkNames)

      ## probe gene mapping
      if (verbose) {
        message("\tProbe-gene mapping")
      }
      platf2 <- platf
      platf2[!is.element(platf, probegene.method[ , "platform"])] <- "MISC"
      esets <- mapply(probeGeneMapping, eset=esets, platform=platf2, method=probegene.method[match(platf2, probegene.method[ , "platform"]), "method"])
  
      ## merge datasets if they are on different platforms
      eset <- platformMerging(esets)
      ## eset is a gene-centric expressionSet object
  
      ## format clinical info
      colnames(Biobase::pData(eset)) <- gsub(" ", ".", colnames(Biobase::pData(eset)))
      cinfo <- intersect(cinfo, colnames(Biobase::pData(eset)))
      ## transform factors into characters
      Biobase::pData(eset) <- data.frame(apply(Biobase::pData(eset), 2, function(x) {
        if (is.factor(x)) { x <- stripWhiteSpace(as.character(x)) }
        return(x)
        }), stringsAsFactors=FALSE)
      Biobase::pData(eset)[!is.na(Biobase::pData(eset)) & (Biobase::pData(eset) == "NA" | Biobase::pData(eset) == "N/A" | Biobase::pData(eset) == "NILL" | Biobase::pData(eset) == "NULL" | Biobase::pData(eset) == "UNKNOWN" | Biobase::pData(eset) == "UNK")] <- NA
      save(list=c("eset"), compress=TRUE, file=dataset.fn)
      if (verbose) {
        message("")
      }
    } else {
      load(file=dataset.fn)
    }
  
    ## special action may be required for some datasets
    switch (datasets[i, "Dataset.ID"],
      "GSE2034" = {
        if ("GSE5327" %in% names(eset.all)) {
          ## merge datasets GSE2034 and GSE5327
          if (verbose) {
            message("\tMerge GSE2034 and GSE5327")
          }
          ## update gene expression data
          cg <- union(rownames(Biobase::exprs(eset.all[["GSE5327"]])), rownames(Biobase::exprs(eset)))
          cs <- c(colnames(Biobase::exprs(eset.all[["GSE5327"]])), colnames(Biobase::exprs(eset)))
          ee <- matrix(NA, nrow=length(cg), ncol=length(cs), dimnames=list(cg, cs))
          ee[rownames(Biobase::exprs(eset.all[["GSE5327"]])), colnames(Biobase::exprs(eset.all[["GSE5327"]]))] <- Biobase::exprs(eset.all[["GSE5327"]])
          ee[rownames(Biobase::exprs(eset)), colnames(Biobase::exprs(eset))] <- Biobase::exprs(eset)
          Biobase::exprs(eset.all[["GSE5327"]]) <- ee
          ## update feature data
          cf <- c(colnames(Biobase::fData(eset.all[["GSE5327"]])), colnames(Biobase::fData(eset)))
          ee <- data.frame(matrix(NA, nrow=length(cg), ncol=length(cf), dimnames=list(cg, cf)))
          ee[rownames(Biobase::fData(eset.all[["GSE5327"]])), colnames(Biobase::fData(eset.all[["GSE5327"]]))] <- Biobase::fData(eset.all[["GSE5327"]])
          ee[rownames(Biobase::fData(eset)), colnames(Biobase::fData(eset))] <- Biobase::fData(eset)
          Biobase::fData(eset.all[["GSE5327"]]) <- ee
          ## update pheno data
          cp <- intersect(colnames(Biobase::pData(eset.all[["GSE5327"]])), colnames(Biobase::pData(eset)))
          ee <- data.frame(matrix(NA, nrow=length(cs), ncol=length(cp), dimnames=list(cs, cp)))
          ee[rownames(Biobase::pData(eset.all[["GSE5327"]])), cp] <- Biobase::pData(eset.all[["GSE5327"]])[ , cp]
          ee[rownames(Biobase::pData(eset)), cp] <- Biobase::pData(eset)[ , cp]
          Biobase::pData(eset.all[["GSE5327"]]) <- ee
          ## rename the merged dataset
          names(eset.all)[names(eset.all) == "GSE2034"] <- "GSE2034.GSE5327"
        } else {
          eset.all <- c(eset.all, eset)
          names(eset.all)[length(eset.all)] <- datasets[i, "Dataset.ID"]
        }
      },
      "GSE5327" = {
        if ("GSE2034" %in% names(eset.all)) {
          ## merge datasets GSE2034 and GSE5327
          if (verbose) {
            message("\tMerge GSE2034 and GSE5327")
          }
          ## update gene expression data
          cg <- union(rownames(Biobase::exprs(eset.all[["GSE2034"]])), rownames(Biobase::exprs(eset)))
          cs <- c(colnames(Biobase::exprs(eset.all[["GSE2034"]])), colnames(Biobase::exprs(eset)))
          ee <- matrix(NA, nrow=length(cg), ncol=length(cs), dimnames=list(cg, cs))
          ee[rownames(Biobase::exprs(eset)), colnames(Biobase::exprs(eset))] <- Biobase::exprs(eset)
          ee[rownames(Biobase::exprs(eset.all[["GSE2034"]])), colnames(Biobase::exprs(eset.all[["GSE2034"]]))] <- Biobase::exprs(eset.all[["GSE2034"]])
          Biobase::exprs(eset.all[["GSE2034"]]) <- ee
          ## update feature data
          cf <- intersect(colnames(Biobase::fData(eset.all[["GSE2034"]])), colnames(Biobase::fData(eset)))
          ee <- data.frame(matrix(NA, nrow=length(cg), ncol=length(cf), dimnames=list(cg, cf)))
          ee[rownames(Biobase::fData(eset)), colnames(Biobase::fData(eset))] <- Biobase::fData(eset)
          ee[rownames(Biobase::fData(eset.all[["GSE2034"]])), colnames(Biobase::fData(eset.all[["GSE2034"]]))] <- Biobase::fData(eset.all[["GSE2034"]])
          Biobase::fData(eset.all[["GSE2034"]]) <- ee
          ## update pheno data
          cp <- intersect(colnames(Biobase::pData(eset.all[["GSE2034"]])), colnames(Biobase::pData(eset)))
          ee <- data.frame(matrix(NA, nrow=length(cs), ncol=length(cp), dimnames=list(cs, cp)))
          ee[rownames(Biobase::pData(eset.all[["GSE2034"]])), cp] <- Biobase::pData(eset.all[["GSE2034"]])[ , cp]
          ee[rownames(Biobase::pData(eset)), cp] <- Biobase::pData(eset)[ , cp]
          Biobase::pData(eset.all[["GSE2034"]]) <- ee
          ## rename the merged dataset
          names(eset.all)[names(eset.all) == "GSE2034"] <- "GSE2034.GSE5327"
        } else {
          eset.all <- c(eset.all, eset)
          names(eset.all)[length(eset.all)] <- datasets[i, "Dataset.ID"]
        }
      },
      "GSE2109" = {
        ## select only breast tumors
        iix <- Biobase::pData(eset)[ , "Anatomical.site"] == "breast"
        Biobase::exprs(eset) <- Biobase::exprs(eset)[ , iix, drop=FALSE]
        Biobase::pData(eset) <- Biobase::pData(eset)[iix, , drop=FALSE]
        eset.all <- c(eset.all, eset)
        names(eset.all)[length(eset.all)] <- datasets[i, "Dataset.ID"]
      },
      "ISDB10278" = {
        ## select only breast tumors
        iix <- Biobase::pData(eset)[ , "tissue"] == "TUMOR"
        Biobase::exprs(eset) <- Biobase::exprs(eset)[ , iix, drop=FALSE]
        Biobase::pData(eset) <- Biobase::pData(eset)[iix, , drop=FALSE]
        eset.all <- c(eset.all, eset)
        names(eset.all)[length(eset.all)] <- datasets[i, "Dataset.ID"]
      },
      {
        ## no action is necessary
        eset.all <- c(eset.all, eset)
        names(eset.all)[length(eset.all)] <- datasets[i, "Dataset.ID"]
      })
    ## eset.all contains a list of gene-centric expression sets
  }
  
  ## log out from inSilicoDb
  try(inSilicoDb::InSilicoLogout())

  ## align clinical information
  if (verbose) {
    message("Update clinical information")
  }
  for (j in 1:length(eset.all)) {
    ## same order for all datasets
    Biobase::pData(eset.all[[j]]) <- Biobase::pData(eset.all[[j]])[ , cinfo, drop=FALSE]
    ## ensure age, survival data are numeric
    tt <- intersect(colnames(Biobase::pData(eset.all[[j]])), c("age", "size", "er", "her2", "pgr", "grade", c(paste(c("t", "e"), rep(c("rfs", "dmfs", "os"), each=2), sep="."))))
    Biobase::pData(eset.all[[j]])[ , tt] <- data.frame(apply(Biobase::pData(eset.all[[j]])[ , tt, drop=FALSE], 2, as.numeric), stringsAsFactors=FALSE)
    ## create dfs data: use rfs when available, dmfs otherwise
    surv.time <- Biobase::pData(eset.all[[j]])[ , "t.rfs"]
    surv.time[is.na(Biobase::pData(eset.all[[j]])[ , "t.rfs"])] <- Biobase::pData(eset.all[[j]])[is.na(Biobase::pData(eset.all[[j]])[ , "t.rfs"]), "t.dmfs"]
    surv.event <- Biobase::pData(eset.all[[j]])[ , "e.rfs"]
    surv.event[is.na(Biobase::pData(eset.all[[j]])[ , "e.rfs"])] <- Biobase::pData(eset.all[[j]])[is.na(Biobase::pData(eset.all[[j]])[ , "e.rfs"]), "e.dmfs"]
    Biobase::pData(eset.all[[j]]) <- cbind(Biobase::pData(eset.all[[j]]), "t.dfs"=surv.time, "e.dfs"=surv.event)
  }

  ## annotate with subtypes
  if (verbose) {
    message("Subtype classification")
  }
  for (j in 1:length(eset.all)) {
    eset.all[[j]] <- subtypeClassification(eset=eset.all[[j]], model=sbt.model)
  }

  ## merge expressionSet objects
  if (verbose) {
    message("Merge datasets")
  }
  eset.merged <- datasetMerging(esets=eset.all,  method=merging.method, standardization=merging.std, nthread=nthread)

  ## identify potential duplicated samples
  duplicates <- duplicateFinder(eset=eset.merged, topvar.genes=topvar.genes, dupl.cor=duplicates.cor, method="spearman", nthread=nthread)
  ## annotate the separate esets and the merged eset
  tt <- sapply(duplicates, paste, collapse="///")
  ## merged eset
  Biobase::pData(eset.merged) <- cbind(Biobase::pData(eset.merged), "duplicates"=NA)
  Biobase::pData(eset.merged)[names(tt), "duplicates"] <- tt
  
  ## individual esets
  if (length(tt) > 0) {
	  eset.all <- lapply(eset.all, function (x, y) {
	    Biobase::pData(x) <- cbind(Biobase::pData(x), "duplicates"=NA)
	    nn <- intersect(rownames(pData(x)), y)
	    Biobase::pData(x)[nn, "duplicates"] <- y[nn]
	    return(x)
	  }, y=tt)
  }
  ## remove duplicates
  if (remove.duplicates) {
    ## duplicates are removed by order of datasets
    ## select sample names to remove
    rmix <- duplicates
    ii <- 1
    while (length(rmix) > ii) {
      rmix <- rmix[!is.element(names(rmix), rmix[[ii]])]
      ii <- ii + 1
    }
    rmix <- unique(unlist(rmix))
    ## merged eset
    keepix <- setdiff(Biobase::sampleNames(eset.merged), rmix)
    Biobase::exprs(eset.merged) <- Biobase::exprs(eset.merged)[ , keepix, drop=FALSE]
    Biobase::pData(eset.merged) <- Biobase::pData(eset.merged)[keepix, , drop=FALSE]
    sclass <- getSubtype(eset=eset.merged, method="class")[keepix]
    sfuzzy <- getSubtype(eset=eset.merged, method="fuzzy")[keepix, , drop=FALSE]
    scrisp <- getSubtype(eset=eset.merged, method="crisp")[keepix, , drop=FALSE]
    eset.merged <- setSubtype(eset=eset.merged, subtype.class=sclass, subtype.fuzzy=sfuzzy, subtype.crisp=scrisp)
    ## individual esets
    eset.all <- lapply(eset.all, function (x, y) {
      keepix <- setdiff(Biobase::sampleNames(x), y)
      Biobase::exprs(x) <- Biobase::exprs(x)[ , keepix, drop=FALSE]
      Biobase::pData(x) <- Biobase::pData(x)[keepix, , drop=FALSE]
      sclass <- getSubtype(eset=x, method="class")[keepix]
      sfuzzy <- getSubtype(eset=x, method="fuzzy")[keepix, , drop=FALSE]
      scrisp <- getSubtype(eset=x, method="crisp")[keepix, , drop=FALSE]
      eset.merged <- setSubtype(eset=x, subtype.class=sclass, subtype.fuzzy=sfuzzy, subtype.crisp=scrisp)
      return(x)
    }, y=rmix)
  }
  
  return (list("merged"=eset.merged, "each"=eset.all))
}



