#################
## loading and changing curatedOvarianData
## Natchar February 23, 2015
#################
`getOvCaData` <- 
  function (resdir="cache", probegene.method, remove.duplicates=TRUE, topvar.genes=1000, duplicate.finder=FALSE, duplicates.cor=0.98, datasets, sbt.model=c("scmgene", "scmod2", "scmod1", "pam50", "ssp2006", "ssp2003"), merging.method=c("union", "intersection"), merging.std=c("quantile", "robust.scaling", "scaling", "none"), nthread=1, verbose=TRUE) {  
    merging.method <- match.arg(merging.method)
    merging.std <- match.arg(merging.std)
    
    ##libraries
    library(curatedOvarianData)
    library(logging)
    library(org.Hs.eg.db)

    if (duplicate.finder){
      ## use Ben's method of finding duplicates: duplicate finder method
      ## Load the filtering rules from patientselection.config
      source(system.file("extdata", "patientselection.config", package = "curatedOvarianData"))
      min.sample.size <<- 1 
      min.number.of.events <<- 1
      min.number.of.genes <<- 1
      tcga.lowcor.outliers <<- list()
      duplicates <<-list()
      remove.samples <<- c(tcga.lowcor.outliers, duplicates)
      add.surv.y <<- NULL
      source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))
      
      names(esets)<-gsub( "_.*$", "", names(esets) )
      #### change sample names of TCGA RNAseq to include .RNA at end (solves conflict between sample names of TCGA during merging)
      num.esets <-length(esets)
      RNAnames <- NULL
      for (n in 1:length(sampleNames(esets[[num.esets]]))){
        RNAnames <- c(RNAnames, sprintf("%s.RNA", sampleNames(esets[[num.esets]])[n]))
      }
      sampleNames(esets[[num.esets]]) <- RNAnames 
      RNAnames <- NULL
      for (n in 1:length(sampleNames(esets[[num.esets-1]]))){
        RNAnames <- c(RNAnames, sprintf("%s.miRNA", sampleNames(esets[[num.esets-1]])[n]))
      }
      sampleNames(esets[[num.esets-1]]) <- RNAnames  
      
      for(i in 1:length(esets)){
        ####### fData to include entrezgeneID
        entrezgene <- list()
        gs <- toTable(org.Hs.egSYMBOL)
        gs <- gs[!is.na(gs[ , "symbol"]) & !duplicated(gs[ , "symbol"]), , drop=FALSE]
        gs <- gs[fData(esets[[i]])[,"gene"], "gene_id"]
        
        fData(esets[[i]])$entrezgene <- gs
        fData(esets[[i]])$ENTREZID<- gs
        rownames(fData(esets[[i]])) <- gs
        rownames(exprs(esets[[i]])) <- gs
        names(fData(esets[[i]]))[names(fData(esets[[i]])) == "gene"] <- "SYMBOL"
        experimentData(esets[[i]])@other$cancer_type <- "ovarian"
        
        ## ovcAngiogenic (Angiogenic vs non-Angiogenic Subtyping)
        data <- t(exprs(esets[[i]]))
        annot <- fData(esets[[i]])
        colnames(data) <- rownames(annot)
        message("doing %s", names(esets)[i])
        angio <- ovcAngiogenic(data = data, annot=annot, gmap="entrezgene", do.mapping = TRUE)
        experimentData(esets[[i]])@other$class <- angio$subtype$subtype
        fuzzy <- experimentData(esets[[i]])@other$fuzzy <- data.frame("Angiogenic" = angio$subtype$Angiogenic.proba, "NonAngiogenic" = angio$subtype$nonAngiogenic.proba)
        rownames(experimentData(esets[[i]])@other$fuzzy) <- rownames(angio$subtype)
        crisp <- t(apply(fuzzy, 1, function(x){
          xx<- array(0, dim=length(x), dimnames=list(names(x)))
          xx[which.max(x)] <- 1
          return(xx)
        }))
        rownames(crisp) <- rownames(angio$subtype)
        experimentData(esets[[i]])@other$crisp <- crisp
        
      }
      eset.all <-list()
      for(e in 1:length(esets)){
        eset.all <- c(eset.all, esets[[1]])
      }
      eset.merged <- datasetMerging(esets = esets, method=merging.method, standardization=merging.std, nthread=nthread) 
      
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

      return(list("merged"=eset.merged, "each"= eset.all))
      ## end of duplicate.finder=T
    } else{
      #########################
      ### use Levi's method of finding duplicates(patientselection.config)
      #########################
    ## Load the filtering rules from patientselection.config
    source(system.file("extdata", "patientselection.config", package = "curatedOvarianData"))
    add.surv.y <- NULL
    source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))
    names(esets)<-gsub( "_.*$", "", names(esets) )
    #### change sample names of TCGA RNAseq to include .RNA at end (solves conflict between sample names of TCGA during merging)
    RNAnames <- NULL
    for (n in 1:length(sampleNames(esets[[16]]))){
      RNAnames <- c(RNAnames, sprintf("%s.RNA", sampleNames(esets[[16]])[n]))
    }
    sampleNames(esets[[16]]) <- RNAnames
    
    for(i in 1:length(esets)){
      ####### fData to include entrezgeneID
      entrezgene <- list()
      gs <- toTable(org.Hs.egSYMBOL)
      gs <- gs[!is.na(gs[ , "symbol"]) & !duplicated(gs[ , "symbol"]), , drop=FALSE]
      gs <- gs[fData(esets[[i]])[,"gene"], "gene_id"]
      
      fData(esets[[i]])$entrezgene <- gs
      fData(esets[[i]])$ENTREZID<- gs
      rownames(fData(esets[[i]])) <- gs
      rownames(exprs(esets[[i]])) <- gs
      names(fData(esets[[i]]))[names(fData(esets[[i]])) == "gene"] <- "SYMBOL"
      experimentData(esets[[i]])@other$cancer_type <- "ovarian"
      
    ## ovcAngiogenic (Angiogenic vs non-Angiogenic Subtyping)
      data <- t(exprs(esets[[i]]))
      annot <- fData(esets[[i]])
      colnames(data) <- rownames(annot)
      angio <- ovcAngiogenic(data = data, annot=annot, gmap="entrezgene", do.mapping = TRUE)
      experimentData(esets[[i]])@other$class <- angio$subtype$subtype
      fuzzy <- experimentData(esets[[i]])@other$fuzzy <- data.frame("Angiogenic" = angio$subtype$Angiogenic.proba, "NonAngiogenic" = angio$subtype$nonAngiogenic.proba)
      rownames(experimentData(esets[[i]])@other$fuzzy) <- rownames(angio$subtype)
      crisp <- t(apply(fuzzy, 1, function(x){
        xx<- array(0, dim=length(x), dimnames=list(names(x)))
        xx[which.max(x)] <- 1
        return(xx)
      }))
      rownames(crisp) <- rownames(angio$subtype)
      experimentData(esets[[i]])@other$crisp <- crisp

    }
    eset.all <-list()
    for(e in 1:length(esets)){
      eset.all <- c(eset.all, esets[[1]])
    }
    eset.merged <- datasetMerging(esets = esets, method=merging.method, standardization=merging.std, nthread=nthread) 
    return(list("merged"=eset.merged, "each"= eset.all))
    
    ## End of duplicate.finder = F
  }
## End of getOvCaData
}

