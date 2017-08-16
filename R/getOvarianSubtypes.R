#' Function to identify ovarian cancer molecular subtypes
#'
#' This function identifies the breast cancer molecular subtypes using a Subtype Clustering Model
#' @param eset an eset object
#' @param subtypeModel a string representing the desired subtyping classification model. Subtyping classification model, can be either "verhaak", "bentink", "tothill", "helland", or "konecny".
#' @param intersectThresh the fraction of genes from the eset that must be present in each verhaak subtype gene set in order for the patients in the eset to be classified according to the verhaak subtypes. Default value is 0.75 and the variable is not relevant for the other subtyping schemes.
#' @return a list with 3 elements. "subtype", which has the subtype of each patient in the eSet for the given subtyping classification model. The second element varies and is specific to the suptyping scheme, it generally provides information regarding how well each patient fit into each subtyping scheme. Annot.eset has information about the eset that was sent.
#' @export
#' @examples
#'
#' eset = loadMetaData("ovarian", "overall")[[1]]
#' subtypeList = getOvarianSubtypes(eset, "verhaak")


getOvarianSubtypes = function(eset, subtypeModel, intersectThresh = 0.75)
{
  #######general functions needed (from genefu)#########
  `rescale` <-
    function(x, na.rm=FALSE, q=0) {
      if(q == 0) {
        ma <- max(x, na.rm=na.rm)
        mi <- min(x, na.rm=na.rm)
      } else {
        ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
        mi <- quantile(x, probs=q/2, na.rm=na.rm)
      }
      xx <- (x - mi) / (ma - mi)
      attributes(xx) <- list("names"=names(x), "q1"=mi,"q2"=ma)
      return(xx)
    }

  `geneid.map` <-
    function(geneid1, data1, geneid2, data2, verbose=FALSE) {

      nn <- names(geneid1)
      geneid1 <- as.character(geneid1)
      names(geneid1) <- nn
      nn <- names(geneid2)
      geneid2 <- as.character(geneid2)
      names(geneid2) <- nn
      if(is.null(names(geneid1))) { names(geneid1) <- dimnames(data1)[[2]] }
      if(!missing(data2) && is.null(names(geneid2))) { names(geneid2) <- dimnames(data2)[[2]] }
      if(!missing(data1) && !missing(geneid1) && !missing(geneid2)) {
        ## remove probes without any measurements
        na.ix <- apply(data1, 2, function(x) { return(all(is.na(x))) })
        data1 <- data1[ , !na.ix, drop=FALSE]
        geneid1 <- geneid1[!na.ix]
      } else { stop("data1, geneid1 and geneid2 parameters are mandatory!") }
      if(!missing(data2)) {
        ## remove probes without any measurements
        na.ix <- apply(data2, 2, function(x) { return(all(is.na(x))) })
        data2 <- data2[ , !na.ix, drop=FALSE]
        geneid2 <- geneid2[!na.ix]
      } else { data2 <- NULL }

      gix1 <- !is.na(geneid1)
      gix2 <- !is.na(geneid2)

      geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
      if(length(geneid.common) == 0) {
        warning("no gene ids in common!")
        return(list("geneid1"=NA, "data1"=NA, "geneid2"=NA, "data2"=NA))
      }

      ## dataset1
      ## probes corresponding to common gene ids
      gg <- names(geneid1)[is.element(geneid1, geneid.common)]
      gid <- geneid1[is.element(geneid1, geneid.common)]
      ## duplicated gene ids
      gid.dupl <- unique(gid[duplicated(gid)])
      gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
      ## unique gene ids
      gid.uniq <- gid[!is.element(gid, gid.dupl)]
      gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
      ## data corresponding to unique gene ids
      datat <- data1[ ,gg.uniq,drop=FALSE]
      ## data for duplicated gene ids
      if(length(gid.dupl) > 0) {
        if(verbose) { message("\ndataset1 duplicates...") }
        ## compute the standard deviation with a penalization on the number of missing values
        ## this should avoid selecting the most variant probe with a lot of missing values
        pena <- apply(X=data1[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=function(x) { return(sum(is.na(x))) })
        pena <- log((nrow(data1) + 1) / (pena + 1)) + 1
        #pena <- 1
        sdr <- drop(apply(X=data1[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=sd, na.rm=TRUE)) * pena
        mysd <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "sd"=sdr)
        mysd <- mysd[order(as.numeric(mysd[ , "sd"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
        mysd <- mysd[!duplicated(mysd[ , "gid"]), , drop=FALSE]
        datat <- cbind(datat, data1[ , mysd[ , "probe"], drop=FALSE])
      }
      data1 <- datat
      geneid1 <- geneid1[dimnames(data1)[[2]]]

      #dataset2
      if(is.null(data2)) {
        #keep arbitrarily the first occurence of each duplicated geneid
        geneid2 <- geneid2[!duplicated(geneid2) & is.element(geneid2, geneid.common)]
      }
      else {
        ## probes corresponding to common gene ids
        gg <- names(geneid2)[is.element(geneid2, geneid.common)]
        gid <- geneid2[is.element(geneid2, geneid.common)]
        ## duplicated gene ids
        gid.dupl <- unique(gid[duplicated(gid)])
        gg.dupl <- names(geneid2)[is.element(geneid2, gid.dupl)]
        ## unique gene ids
        gid.uniq <- gid[!is.element(gid, gid.dupl)]
        gg.uniq <- names(geneid2)[is.element(geneid2, gid.uniq)]
        ## data corresponding to unique gene ids
        datat <- data2[ ,gg.uniq,drop=FALSE]
        ## data for duplicated gene ids
        if(length(gid.dupl) > 0) {
          if(verbose) { message("\ndataset2 duplicates...") }
          ## compute the standard deviation with a penalization on the number of missing values
          ## this should avoid selecting the most variant probe with a lotof missing values
          pena <- apply(X=data2[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=function(x) { return(sum(is.na(x))) })
          pena <- log((nrow(data2) + 1) / (pena + 1)) + 1
          #pena <- 1
          sdr <- drop(apply(X=data2[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=sd, na.rm=TRUE)) * pena
          mysd <- cbind("probe"=gg.dupl, "gid"=geneid2[gg.dupl], "sd"=sdr)
          mysd <- mysd[order(as.numeric(mysd[ , "sd"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
          mysd <- mysd[!duplicated(mysd[ , "gid"]), , drop=FALSE]
          datat <- cbind(datat, data2[ , mysd[ , "probe"], drop=FALSE])
        }
        data2 <- datat
        geneid2 <- geneid2[dimnames(data2)[[2]]]
      }

      #same order for the two datasets
      rix <- match(geneid2, geneid1)
      geneid1 <- geneid1[rix]
      data1 <- data1[ ,rix,drop=FALSE]
      return(list("geneid1"=geneid1, "data1"=data1, "geneid2"=geneid2, "data2"=data2))
    }

  ####end general functions needed#############

  subtypeModel = tolower(subtypeModel)
  if(subtypeModel == "verhaak"){
    #getVerhaakSubtypes <- function(eset, verhaakSheetOne, verhaakSheetSeven)
      data("verhaakSheetOne")
      data("verhaakSheetSeven")
      genesets <- lapply(levels(verhaakSheetSeven$CLASS), function(y) as.character(verhaakSheetSeven[verhaakSheetSeven$CLASS==y,1]))
      names(genesets) <-  levels(verhaakSheetSeven$CLASS)
      fullSetLengths = lengths(genesets)

      # For ssGSEA scores for the new samples, use the intersecting genes
      genesets <- lapply(genesets, function(x) intersect(x, fData(eset)$gene))
      intersectSetLengths = lengths(genesets)
      returnNull = FALSE
      if(sum(intersectSetLengths > intersectThresh*fullSetLengths)  == 0)
        returnNull = TRUE

      if(returnNull == FALSE){
        ## Determine the ssGSEA cutoffs for the IMR and MES subtypes
        supplementary.tcga.discovery <- verhaakSheetOne[ verhaakSheetOne$DATASET=="TCGA-discovery", c("ID", "SUBTYPE") ]
        supplementary.tcga.discovery <- supplementary.tcga.discovery[ supplementary.tcga.discovery$SUBTYPE %in% c("Mesenchymal", "Immunoreactive"), ]

        IMR.threshold <- 0.63 #min(tcga.gsva.out.with.published.subtype$IMR[ tcga.gsva.out.with.published.subtype$SUBTYPE=="Immunoreactive" ])
        MES.threshold <- 0.56 #min(tcga.gsva.out.with.published.subtype$MES[ tcga.gsva.out.with.published.subtype$SUBTYPE=="Mesenchymal" ])

        expression.matrix <- exprs(eset)
        rownames(expression.matrix) <- fData(eset)$gene
        gsva.out <- gsva(expression.matrix, genesets, method="ssgsea", tau=0.75, parallel.sz=4, mx.diff=FALSE, ssgsea.norm=FALSE, verbose = FALSE)
        gsva.out <- t(gsva.out)
        gsva.out <- apply(gsva.out, 2, function(x) ( x - min(x) ) / ( max(x) - min(x) ))

        ## Classify each sample according to the max ssGSEA subtype score, using the scheme provided in the methods.
        subclasses <- apply(gsva.out, 1, function(x) {
          if(x[which(colnames(gsva.out)=="IMR")] > IMR.threshold && x[which(colnames(gsva.out)=="MES")] > MES.threshold) {
            return(c("IMR", "MES")[which.max(x[c("IMR", "MES")])])
          }
          colnames(gsva.out)[which.max(x)]
        })

        subclasses <- factor(subclasses, levels=levels(verhaakSheetSeven$CLASS))
        ## Append a new column for Verhaak subtypes
        #eset$subtypes <- subclasses

      }else{
        #Not enough of the genes required to determine the patients subtypes are present
        subclasses <- NULL
        gsva.out = NULL
      }
      subtypeInfoList = list(subtypes = subclasses, gsva.out=gsva.out, Annotated.eset=eset)
      #return(list(Annotated.eset=eset, gsva.out=gsva.out))
  }else if(subtypeModel == "bentink"){
      ## Classify new samples
    #print("hi")
      expression.matrix <- t(exprs(eset))
      annot <- fData(eset)
      colnames(annot)[which(colnames(annot) == "EntrezGene.ID")] <- "entrezgene"

      `ovcAngiogenic` <-
        function(data, annot, hgs, gmap=c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene"), do.mapping=FALSE, verbose=FALSE) {
          gmap <- match.arg(gmap)
          if(missing(hgs)) { hgs <- rep(TRUE, nrow(data)) }
          if(do.mapping) {
            if(!is.element(gmap, colnames(annot))) { stop("gmap is not a column of annot!") }
            if(verbose) { message("the most variant probe is selected for each gene") }
            sigt <- sigOvcAngiogenic[order(abs(sigOvcAngiogenic[ ,"weight"]), decreasing=FALSE), ,drop=FALSE]
            sigt <- sigt[!duplicated(sigt[ ,gmap]), ,drop=FALSE]
            gid2 <- sigt[ ,gmap]
            names(gid2) <- rownames(sigt)
            gid1 <- annot[ ,gmap]
            names(gid1) <- colnames(data)
            rr <- geneid.map(geneid1=gid1, data1=data, geneid2=gid2)
            data <- rr$data1
            annot <- annot[colnames(data), ,drop=FALSE]
            sigt <- sigt[names(rr$geneid2), ,drop=FALSE]
            pold <- colnames(data)
            pold2 <- rownames(sigt)
            colnames(data) <- rownames(annot) <- rownames(sigt) <- paste("geneid", annot[ ,gmap], sep=".")
            mymapping <- c("mapped"=nrow(sigt), "total"=nrow(sigOvcAngiogenic))
            myprobe <- data.frame("probe"=pold, "gene.map"=annot[ ,gmap], "new.probe"=pold2)
          } else {
            gix <- intersect(rownames(sigOvcAngiogenic), colnames(data))
            if(length(gix) < 2) { stop("data do not contain enough gene from the ovcTCGA signature!") }
            data <- data[ ,gix,drop=FALSE]
            annot <- annot[gix, ,drop=FALSE]
            mymapping <- c("mapped"=length(gix), "total"=nrow(sigOvcAngiogenic))
            myprobe <- data.frame("probe"=gix, "gene.map"=annot[ ,gmap], "new.probe"=gix)
            sigt <- sigOvcAngiogenic[gix, ,drop=FALSE]
          }

          #data(modelOvcAngiogenic)
          ss <- calcSigScore(x=data.frame("probe"=colnames(data), "EntrezGene.ID"=annot[ ,gmap], "coefficient"=sigt[ ,"weight"]), data=data, annot=annot, doMapping=FALSE, signed=TRUE)$score
          ## rescale only with the high grade, late stage, serous (hgs) patients
          rr <- rescale(ss[hgs], q=0.05, na.rm=TRUE)
          ## rescale the whole dataset
          pscore <- ((ss - attributes(rr)$q1) / (attributes(rr)$q2 - attributes(rr)$q1) - 0.5) * 2
          emclust.ts <- mclust::estep(modelName="E", data=pscore, parameters=modelOvcAngiogenic)
          dimnames(emclust.ts$z) <- list(names(pscore), c("Angiogenic.proba", "nonAngiogenic.proba"))
          class.ts <- mclust::map(emclust.ts$z, warn=FALSE)
          names(class.ts) <- names(pscore)
          sbt.ts <- class.ts
          sbt.ts[class.ts == 1] <- "Angiogenic"
          sbt.ts[class.ts == 2] <- "nonAngiogenic"
          sbts <- data.frame("subtype.score"=pscore, "subtype"=sbt.ts, emclust.ts$z)
          prisk <- as.numeric(sbts[ ,"subtype"] == "Angiogenic")
          names(prisk) <- names(pscore) <- rownames(data)
          return (list("score"=pscore, "risk"=prisk, "mapping"=mymapping, "probe"=myprobe, "subtype"=sbts))
        }

      angio <- ovcAngiogenic(data = expression.matrix, annot=annot, gmap="entrezgene", do.mapping = TRUE)
      eset$subtypes <- angio$subtype$subtype
      #eset$Bentink.subtypes <- angio$subtype$subtype
      #return(list(Annotated.eset=eset, angio=angio))
      subtypeInfoList = list(subtypes = eset$subtypes, angio=angio, Annotated.eset=eset)
  }else if(subtypeModel == "helland"){

      supplementary.type.1 <- gdata::read.xls(system.file("extdata", "journal.pone.0018064.s015.XLS", package="metaGx"), sheet=1)
      supplementary.type.2 <- gdata::read.xls(system.file("extdata", "journal.pone.0018064.s015.XLS", package="metaGx"), sheet=2)
      supplementary.type.4 <- gdata::read.xls(system.file("extdata", "journal.pone.0018064.s015.XLS", package="metaGx"), sheet=3)
      supplementary.type.5 <- gdata::read.xls(system.file("extdata", "journal.pone.0018064.s015.XLS", package="metaGx"), sheet=4)
      supplementary.tables <- list(C1=supplementary.type.1, C2=supplementary.type.2, C4=supplementary.type.4, C5=supplementary.type.5)

      entrez.id.logFC.list <- lapply(supplementary.tables, function(x) {
        ## Use the supplementary table's listed probe id and gene name to determine the Entrez ID
        # If there is only one EntrezID that maps to a probe in hgu133plus2.db, use that Entrez ID.
        # If there are multiple EntrezIDs that map to a probe, then use the EntrezID (if any) that corresponds to the provided gene symbol.
        current.mapping <- suppressWarnings(AnnotationDbi::select(hgu133plus2.db, as.character(x$ID), c("ENTREZID", "SYMBOL")))
        #current.mapping <- suppressWarnings(AnnotationDbi::select(org.Hs.eg.db, as.character(x$ID), c("ENTREZID", "SYMBOL")))
        current.mapping <- current.mapping[ !is.na(current.mapping$ENTREZID), ]
        colnames(x)[1:2] <- c("PROBEID", "SYMBOL")
        mappings.with.unique.probeid <- current.mapping[ !(current.mapping$PROBEID %in% current.mapping$PROBEID[duplicated(current.mapping$PROBEID)]),]
        mappings.with.duplicate.probeid <- current.mapping[ current.mapping$PROBEID %in% current.mapping$PROBEID[duplicated(current.mapping$PROBEID)],]
        mappings.with.duplicate.probeid <- merge(x, mappings.with.duplicate.probeid, by=c("PROBEID", "SYMBOL"))[, c("PROBEID", "ENTREZID", "SYMBOL")]
        mappings.with.duplicate.probeid <- unique(mappings.with.duplicate.probeid)
        current.mapping <- rbind(mappings.with.unique.probeid, mappings.with.duplicate.probeid)
        to.return <- merge(x, current.mapping, by="PROBEID")[, c("ENTREZID", "PROBEID", "logFC")]
        return(to.return)
      })

      gene.union.in.supplementary <- Reduce(function(x,y) union(x, y), lapply(entrez.id.logFC.list, function (x) x$ENTREZID))

      intersecting.entrez.ids <- intersect(gene.union.in.supplementary, fData(eset)$EntrezGene.ID)

      # Only keep genes present in both the supplementary and this eset
      entrez.id.logFC.list <- lapply(entrez.id.logFC.list, function(x) x[x$ENTREZID %in% intersecting.entrez.ids, ])

      subtype.scores <- sapply(entrez.id.logFC.list, function(x) {
        ordered.expression.subset <- exprs(eset)[match(x$ENTREZID, fData(eset)$EntrezGene.ID),]

        return(apply(ordered.expression.subset, 2, function(y) sum((y * x$logFC))))
      })
      old.rownames <- rownames(subtype.scores)
      # Scale to mean=0, variance=1
      subtype.scores <- apply(subtype.scores, 2, scale)
      rownames(subtype.scores) <- old.rownames
      subclasses <- factor(colnames(subtype.scores)[apply(subtype.scores, 1, which.max)], levels=names(supplementary.tables))
      #eset$subtypes <- subclasses
      subtypeInfoList = list(subtypes = subclasses, subtype.scores=subtype.scores, Annotated.eset=eset)
      #return(list(Annotated.eset=eset, subtype.scores=subtype.scores))

  }else if(subtypeModel == "konecny"){

      expression.matrix <- exprs(eset)
      # Rescale per gene
      expression.matrix <- t(scale(t(expression.matrix)))
      ## take only the first entry for ambiguous probes
      pn <- sub("geneid.", "", rownames(expression.matrix))
      pn <- sapply(strsplit(pn, ","), function (x) { return (x[[1]]) })
      rownames(expression.matrix) <- pn

      ## Load centroids defined in Konecny et al., 2014
      supplementary.data <- gdata::read.xls(system.file("extdata", "jnci_JNCI_14_0249_s05.xls", package="metaGx"), sheet=4)

      ## Classify using nearest centroid with Spearman's rho

      # Subset supplementary.data to consist of centroids with intersecting genes
      # For genes with multiple probesets, take the mean
      centroids <- supplementary.data[,c(2,4:7)]
      centroids[,2:5] <- sapply(centroids[,2:5], function(x) ave(x, centroids$EntrezGeneID, FUN=mean))
      centroids <- unique(centroids)
      rownames(centroids) <- centroids$EntrezGeneID
      centroids <- centroids[,-1]

      rownames(expression.matrix) = eset@featureData@data$EntrezGene.ID
      intersecting.entrez.ids <- intersect(rownames(expression.matrix), rownames(centroids))
      centroids[rownames(centroids) %in% intersecting.entrez.ids,]
      centroids <- centroids[as.character(intersecting.entrez.ids),]
      expression.matrix <- expression.matrix[as.character(intersecting.entrez.ids),]

      expression.matrix <- as.data.frame(expression.matrix)
      spearman.cc.vals <- sapply(centroids, function(x) sapply(expression.matrix, function(y) cor(x, y , method="spearman")))

      subclasses <- apply(spearman.cc.vals, 1, function(x) as.factor(colnames(spearman.cc.vals)[which.max(x)]))

      subclasses <- factor(subclasses, levels=colnames(centroids))

      #eset$subtypes <- subclasses
      subtypeInfoList = list(subtypes = subclasses, spearman.cc.vals=spearman.cc.vals, Annotated.eset=eset)
      #return(list(Annotated.eset=eset, spearman.cc.vals=spearman.cc.vals))

  }else if(subtypeModel == "tothill"){

      gene.mapping=c("Entrez.ID")
      #gene.mapping <- match.arg(gene.mapping)

      ## Load train data with predefined class labels
      supp.table.2 <- read.table(system.file("extdata", "tothill.supptable.probes.entrez.txt", package="metaGx"), header=TRUE)
      #supp.table.2 <- read.table(system.file(file.path("extdata", "tothill.supptable.probes.entrez.txt", package="metaGx"), header=TRUE))
      supplementary.probesets <- as.character(supp.table.2$Probe.ID)
      supplementary.entrez.ids <- unique(supp.table.2[supp.table.2$Entrez.ID != "---",]$Entrez.ID)

      # Get the probeset - Entrez ID mapping for the platform used in the Tothill et al. study
      #probe.entrez.mappings <- as.list(hgu133plus2ENTREZID[mappedkeys(hgu133plus2ENTREZID)])
      #supplementary.entrez.ids <- unlist(probe.entrez.mappings[supplementary.probesets])

      ## Train a diagonal linear discriminant classifier using the Tothill data set and overlapping probesets / entrez gene ids
      # Get the expression matrix of this eset and the Tothill eset, with columns named by Entrez gene ids
      #if(gene.mapping == "Entrez.ID") {

        tothill.eset <- GSE9891
        bestProbes = getBestProbes(tothill.eset)
        tothill.expression.matrix <- exprs(tothill.eset)
        tothill.expression.matrix = tothill.expression.matrix[bestProbes, ]
        rownames(tothill.expression.matrix) <- fData(tothill.eset)$EntrezGene.ID[match(rownames(tothill.expression.matrix), rownames(fData(tothill.eset)))]
        colnames(tothill.expression.matrix) <- sub("X", "", pData(tothill.eset)$alt_sample_name)

        expression.matrix <- exprs(eset)
        bestProbes = getBestProbes(eset)
        expression.matrix = expression.matrix[bestProbes, ]
        rownames(expression.matrix) <- fData(eset)$EntrezGene.ID[match(rownames(expression.matrix), rownames(fData(eset)))]

        intersecting.entrez.ids <- intersect(supplementary.entrez.ids, intersect(rownames(expression.matrix), rownames(tothill.expression.matrix)))
        expression.matrix <- expression.matrix[rownames(expression.matrix) %in%  intersecting.entrez.ids,]
        tothill.expression.matrix <- tothill.expression.matrix[rownames(tothill.expression.matrix) %in%  intersecting.entrez.ids,]
      #} else if(gene.mapping == "Probe.ID") {
        #expression.matrix <- exprs(eset)
        #rownames(expression.matrix) <- fData(eset)$probeset[match(rownames(expression.matrix), rownames(fData(eset)))]

        #tothill.eset <- esets$GSE9891
        #tothill.expression.matrix <- exprs(tothill.eset)
        #rownames(tothill.expression.matrix) <- fData(tothill.eset)$probeset[match(rownames(tothill.expression.matrix), rownames(fData(tothill.eset)))]

        #colnames(tothill.expression.matrix) <- sub("X", "", pData(tothill.eset)$alt_sample_name)

        #intersecting.probesets <- intersect(supplementary.probesets, intersect(rownames(expression.matrix), rownames(tothill.expression.matrix)))
        #expression.matrix <- expression.matrix[rownames(expression.matrix) %in%  intersecting.probesets,]
        #tothill.expression.matrix <- tothill.expression.matrix[rownames(tothill.expression.matrix) %in%  intersecting.probesets,]
      #}
      #Transpose matrices, so each row is a sample and columns are genes
      expression.matrix <- t(expression.matrix)
      tothill.expression.matrix <- t(tothill.expression.matrix)

      ## Train and classify with diagonal LDA based on the cla
      supplementary.classes <- read.table(system.file(file.path("extdata", "tothill.supptable.1.classes.txt"), package="metaGx"), header=TRUE)
      supplementary.classes$group <- as.character(supplementary.classes$group)
      supplementary.classes <- supplementary.classes[supplementary.classes$group != "NC",]
      supplementary.classes$group <- as.factor(supplementary.classes$group)
      levels(supplementary.classes$group) <- paste0("C", levels(supplementary.classes$group))
      rownames(supplementary.classes) <- supplementary.classes$ID
      supplementary.classes <- supplementary.classes[,-1,drop=FALSE]

      tothill.train.data <- merge(supplementary.classes, tothill.expression.matrix, by="row.names")
      rownames(tothill.train.data) <- tothill.train.data$Row.names
      tothill.train.data <- tothill.train.data[,-1]
      print("Training a DLDA classifier based on common genes...")
      trained.dlda <- HiDimDA::Dlda(data=tothill.train.data[,-1], grouping=tothill.train.data$group)
      print("Finished training.")
      subclasses <- predict(trained.dlda, expression.matrix, grpcodes=levels(tothill.train.data$group))$class
      eset$subtypes <- subclasses
      #eset$Tothill.subtypes <- subclasses
      subtypeInfoList = list(subtypes = subclasses, trained.dlda=trained.dlda, Annotated.eset=eset)
  }

  #changed the return as well
  return(subtypeInfoList)
}


