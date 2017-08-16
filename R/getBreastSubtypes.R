#' Function to identify breast cancer molecular subtypes using the Subtype Clustering Model
#'
#' This function identifies the breast cancer molecular subtypes using a Subtype Clustering Model
#' @param data Matrix of gene expressions with samples in columns and probes in rows, dimnames being properly defined.
#' @param subtypeModel a string representing the desired subtyping classification model. Subtyping classification model, can be either "scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003", "intClust", "AIMS", or "claudinLow".
#' @param annot Matrix of annotations with at least one column named "EntrezGene.ID" (for ssp, scm, AIMS, and claudinLow models) or "Gene.Symbol" (for the intClust model), dimnames being properly defined.
#' @param doMapping TRUE if the mapping through Entrez Gene ids must be performed (in case of ambiguities, the most variant probe is kept for each gene), FALSE otherwise.
#' @return a list with 3 elements. "subtype", which has Subtypes identified by the subtyping classification model. "subtype.proba", which has the probabilities that a patient belongs to each subtype estimated by the subtyping classification model. "subtype.crisp", which has Crisp classes identified by the subtyping classification model.
#' @export
#' @examples
#' eset = loadMetaData("breast", "overall")[[1]]
#' expressionData = eset@assayData$exprs
#' subtypePredictions = getBreastSubtypes(data = expressionData, subtypeModel = "scmod2", annot = eset@featureData@data, doMapping = TRUE)


getBreastSubtypes = function (data, subtypeModel=c("scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003", "intclust", "aims","claudinlow"), annot, doMapping=FALSE) {

  subtypeModel = tolower(subtypeModel)
  #data(list=data(package="metaGx")[[3]][,3])
  #strEsets <- as.character(data(package="metaGx")[[3]][,3])
  #remInds = c(which(grepl("annot.", strEsets)), which(grepl("data.", strEsets)), which(grepl("demo", strEsets)))
  #strEsets = strEsets[-remInds]
  #esets <- list()
  #for (strEset in strEsets){
  #  eset <- get(strEset)
  #}

  verbose = FALSE
  data = t(data)
    #if(getRversion() >= "2.15.1")  utils::globalVariables(c("scmgene.robust","scmod2.robust","pam50.robust","ssp2006.robust","ssp2003.robust","claudinLowData"))

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
    
    medianCtr <- function(x){
      annAll <- dimnames(x)
      medians <- apply(x,1,median,na.rm=T)
      x <- t(scale(t(x),center=medians,scale=F))
      dimnames(x) <- annAll
      return(x)
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

    `subtype.cluster.predict` <-
      function(subtypeModel, data, annot, doMapping=FALSE, mapping, do.prediction.strength=FALSE, do.BIC=FALSE, plot=FALSE, verbose=FALSE) {
        #require(mclust)
        if(missing(data) || missing(annot)) { stop("data, and annot parameters must be specified") }

        sbtn <- c("ER-/HER2-", "HER2+", "ER+/HER2-")
        sbtn2 <- c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif")

        if(is.list(subtypeModel)) {
          ## retrieve model
          subtype.c <- subtypeModel[!is.element(names(subtypeModel), c("cutoff.AURKA", "mod"))]
          model.name <- subtype.c$parameters$variance$modelName
          cc <- subtypeModel$gaussian.AURKA
          mq <- subtypeModel$rescale.q
          m.mod <- subtypeModel$mod
        } else {
          ## read model file
          rr <- readLines(con=subtypeModel, n=11)[-1]
          nn <- unlist(lapply(X=rr, FUN=function(x) { x <- unlist(strsplit(x=unlist(strsplit(x=x, split=":"))[1], split=" ")); x <- x[length(x)]; return(x); }))
          m.param <- c(list(rr[[1]]), lapply(X=rr[-1], FUN=function(x) { x <- as.numeric(unlist(strsplit(x=unlist(strsplit(x=x, split=":"))[2], split=" "))); x <- x[!is.na(x)]; return(x)}))
          names(m.param) <- nn
          cn <- unlist(lapply(strsplit(nn[grep(pattern="mean", x=nn)], split="[.]"), FUN=function(x) { return(x[[2]]) }))
          m.mod <- read.m.file(subtypeModel, comment.char="#")
          #construct a fake mclust object with the parameters of the model
          subtype.c <- NULL
          tt <- m.param$pro
          names(tt) <- cn
          subtype.c$parameters$pro <- tt
          tt <- sapply(X=m.param[grep(pattern="mean", x=nn)], FUN=function(x) { return(x) })
          dimnames(tt) <- list(names(m.mod)[1:2], cn)
          subtype.c$parameters$mean <- tt
          subtype.c$parameters$variance$modelName <- model.name <- m.param$modelname
          subtype.c$parameters$variance$d <- 2
          subtype.c$parameters$variance$G <- 3
          tt <- matrix(0, ncol=2, nrow=2, dimnames=list(names(m.mod)[1:2], names(m.mod)[1:2]))
          diag(tt) <- m.param$sigma
          subtype.c$parameters$variance$sigma <- array(tt, dim=c(2,2,3), dimnames=list(names(m.mod)[1:2], names(m.mod)[1:2], cn))
          subtype.c$parameters$variance$Sigma <- tt
          subtype.c$parameters$variance$scale <- m.param$scale
          subtype.c$parameters$variance$shape <- m.param$shape
          cc <- c("mean"=m.param$gaussian.AURKA.mean, "sigma"=m.param$gaussian.AURKA.sigma)
          mq <- m.param$rescale.q
        }
        do.scale <- ifelse(is.na(mq), FALSE, TRUE)
        sbt <- rep(NA, nrow(data))
        names(sbt) <- dimnames(data)[[1]]
        sbt.proba <- matrix(NA, nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], sbtn))

        sigs.esr1 <- calcSigScore(x=m.mod$ESR1, data=data, annot=annot, doMapping=doMapping, mapping=mapping, verbose=FALSE)
        sigs.erbb2 <- calcSigScore(x=m.mod$ERBB2, data=data, annot=annot, doMapping=doMapping, mapping=mapping, verbose=FALSE)
        sigs.aurka <- calcSigScore(x=m.mod$AURKA, data=data, annot=annot, doMapping=doMapping, mapping=mapping, verbose=FALSE)
        ## signature scores
        dd <- cbind("ESR1"=sigs.esr1$score, "ERBB2"=sigs.erbb2$score, "AURKA"=sigs.aurka$score)
        ## mapping
        mymap <- list("ESR1"=sigs.esr1$probe, "ERBB2"=sigs.erbb2$probe, "AURLA"=sigs.aurka$probe)
        cln <- dimnames(subtype.c$parameters$mean)[[2]] <- as.character(1:ncol(subtype.c$parameters$mean))

        if(do.scale) {
          ## the rescaling needs a large sample size!!!
          ## necessary if we want to validate the classifier using a different dataset
          ## the estimation of survival probabilities depends on the scale of the score
          dd <- apply(dd, 2, function(x) { return((rescale(x, q=mq, na.rm=TRUE) - 0.5) * 2) })
        }
        rownames(dd) <- rownames(data)
        dd2 <- dd

        cc.ix <- complete.cases(dd[ , c("ESR1", "ERBB2"), drop=FALSE])
        if(all(!cc.ix)) {
          ps.res <- ps.res2 <- BIC.res <- NULL
          if(do.prediction.strength) {
            tt <- rep(NA, length(sbtn))
            names(tt) <- sbtn
            tt2 <- rep(NA, length(sbtn2))
            names(tt2) <- sbtn2
            ps.res <- list("ps"=NA, "ps.cluster"=tt, "ps.individual"=sbt)
            ps.res2 <- list("ps"=NA, "ps.cluster"=tt2, "ps.individual"=sbt)
          }
          if(do.BIC) {
            BIC.res <- rep(NA, 10)
            names(BIC.res) <- 1:10
          }
          return(list("subtype"=sbt, "subtype.proba"=sbt.proba, "prediction.strength"=ps.res, "BIC"=BIC.res, "subtype2"=sbt, "prediction.strength2"=ps.res2))
        }
        dd <- dd[cc.ix, , drop=FALSE]

        emclust.ts <- mclust::estep(modelName=model.name, data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], parameters=subtype.c$parameters)
        dimnames(emclust.ts$z) <- list(dimnames(dd)[[1]], cln)
        class.ts <- mclust::map(emclust.ts$z, warn=FALSE)
        names(class.ts) <- dimnames(dd)[[1]]
        uclass <- sort(unique(class.ts))
        uclass <- uclass[!is.na(uclass)]

        ps.res <- ps.res2 <- NULL
        if(do.prediction.strength) {
          if(nrow(dd) < 10) {
            warning("at least 10 observations are required to compute the prediction strength!")
            tt <- rep(NA, length(sbtn))
            names(tt) <- sbtn
            tt2 <- rep(NA, nrow(dd2))
            names(tt2) <- dimnames(dd2)[[1]]
            ps.res <- list("ps"=0, "ps.cluster"=tt, "ps.individual"=tt2)
            tt <- rep(NA, length(sbtn2))
            names(tt) <- sbtn2
            ps.res2 <- list("ps"=0, "ps.cluster"=tt, "ps.individual"=tt2)
          } else {
            ## computation of the prediction strength of the clustering
            rr3 <- mclust::Mclust(data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], modelNames=model.name, G=3)
            ## redefine classification to be coherent with subtypes
            uclass <- sort(unique(rr3$classification))
            uclass <- uclass[!is.na(uclass)]
            if(length(uclass) != 3) {
              warning("less than 3 subtypes are identified!")
              tt <- rep(NA, length(sbtn))
              names(tt) <- sbtn
              tt2 <- rep(NA, nrow(dd2))
              names(tt2) <- dimnames(dd2)[[1]]
              ps.res <- list("ps"=0, "ps.cluster"=tt, "ps.individual"=tt2)
              tt <- rep(NA, length(sbtn2))
              names(tt) <- sbtn2
              ps.res2 <- list("ps"=0, "ps.cluster"=tt, "ps.individual"=tt2)
            } else {
              mm <- NULL
              for(i in 1:length(uclass)) {
                mm <- c(mm, median(dd[rr3$classification == uclass[i],"ERBB2"], na.rm=TRUE) )
              }
              nclass <-  uclass[order(mm, decreasing=TRUE)[1]]
              mm <- NULL
              for(i in 1:length(uclass[-nclass])) {
                mm <- c(mm, median(dd[rr3$classification == uclass[-nclass][i],"ESR1"], na.rm=TRUE) )
              }
              nclass <- c(uclass[-nclass][order(mm, decreasing=TRUE)[2]], nclass, uclass[-nclass][order(mm, decreasing=TRUE)[1]])
              ## nclass contains the new order
              ncl <- rr3$classification
              for(i in 1:length(uclass)) {
                ncl[rr3$classification == nclass[i]] <- i
              }
              ## use the previously computed model to fit a new model in a supervised manner
              myclass <- mclust::unmap(ncl)
              dimnames(myclass) <-  list(dimnames(dd)[[1]], sbtn)
              mclust.tr <- mclust::mstep(modelName=model.name, data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], z=myclass)
              dimnames(mclust.tr$z) <- dimnames(myclass)
              emclust.tr <- mclust::estep(modelName=model.name, data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], parameters=mclust.tr$parameters)
              dimnames(emclust.tr$z) <- dimnames(myclass)
              class.tr <- mclust::map(emclust.tr$z, warn=FALSE)
              names(class.tr) <- dimnames(dd)[[1]]
              ## prediction strength
              ps.res <- ps.cluster(cl.tr=class.ts, cl.ts=class.tr, na.rm=TRUE)
              names(ps.res$ps.cluster) <- sbtn
              ## check for missing values in ps.individual
              tt2 <- rep(NA, nrow(dd2))
              names(tt2) <- dimnames(dd2)[[1]]
              tt2[names(ps.res$ps.individual)] <- ps.res$ps.individual
              ps.res$ps.individual <- tt2
              ## prediction strength with the separation in high and low proliferative tumors
              ## since proliferation is a continuum we fit a Gaussian using AURKA expression of the ER+/HER2- tumors
              ## refitted model
              tt <- mclust::Mclust(dd[complete.cases(class.tr, dd[ , "AURKA"]) & class.tr == 3, "AURKA"], modelNames="E", G=1)
              gauss.prolif <- c("mean"=tt$parameters$mean, "sigma"=tt$parameters$variance$sigmasq)
              class.tr2 <- class.tr
              class.tr2[class.tr == 3] <- NA
              ## probability that tumor is highly proliferative
              pprolif <- pnorm(q=dd[ , "AURKA"], mean=gauss.prolif["mean"], sd=gauss.prolif["sigma"], lower.tail=TRUE)
              ## high proliferation
              class.tr2[class.tr == 3 & pprolif >= 0.5 & complete.cases(class.tr, pprolif)] <- 3
              ## low proliferation
              class.tr2[class.tr == 3 & pprolif < 0.5 & complete.cases(class.tr, pprolif)] <- 4
              ## existing model
              tt <- mclust::Mclust(dd[complete.cases(class.ts, dd[ , "AURKA"]) & class.ts == 3, "AURKA"], modelNames="E", G=1)
              gauss.prolif <- c("mean"=tt$parameters$mean, "sigma"=tt$parameters$variance$sigmasq)
              class.ts2 <- class.ts
              class.ts2[class.ts == 3] <- NA
              ## probability that tumor is highly proliferative
              pprolif <- pnorm(q=dd[ , "AURKA"], mean=gauss.prolif["mean"], sd=gauss.prolif["sigma"], lower.tail=TRUE)
              ## high proliferation
              class.ts2[class.ts == 3 & pprolif >= 0.5 & complete.cases(class.ts, pprolif)] <- 3
              ## low proliferation
              class.ts2[class.ts == 3 & pprolif < 0.5 & complete.cases(class.ts, pprolif)] <- 4
              ## compute the prediction strength
              ps.res2 <- ps.cluster(cl.tr=class.ts2, cl.ts=class.tr2, na.rm=TRUE)
              names(ps.res2$ps.cluster) <- sbtn2
              ## check for missing values in ps.individual
              tt2 <- rep(NA, nrow(dd2))
              names(tt2) <- dimnames(dd2)[[1]]
              tt2[names(ps.res2$ps.individual)] <- ps.res2$ps.individual
              ps.res2$ps.individual <- tt2
            }
          }
        }

        BIC.res <- NULL
        if(do.BIC) {
          if(nrow(dd) >= 10) { BIC.res <- mclust::mclustBIC(data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], modelNames=c(model.name), G=1:10)[ ,model.name] } else { warning("at least 10 observations are required to compute the BIC!") }
        }

        ## subtypes
        sbt[names(class.ts)] <- sbtn[class.ts]
        sbt.proba[dimnames(emclust.ts$z)[[1]], ] <- emclust.ts$z
        ## discriminate between luminal A and B using AURKA
        gauss.prolif <- cc
        sbt2 <- sbt
        sbt2[sbt == sbtn[3]] <- NA
        ## probability that tumor is highly proliferative
        pprolif <- pnorm(q=dd2[ , "AURKA"], mean=gauss.prolif["mean"], sd=gauss.prolif["sigma"], lower.tail=TRUE)
        ## high proliferation
        sbt2[sbt == sbtn[3] & pprolif >= 0.5 & complete.cases(sbt, pprolif)] <- sbtn2[3]
        ## low proliferation
        sbt2[sbt == sbtn[3] & pprolif < 0.5 & complete.cases(sbt, pprolif)] <- sbtn2[4]
        ## subtype probabilities for luminal B and A
        sbt.proba2 <- matrix(NA, nrow(data), ncol=length(sbtn2), dimnames=list(dimnames(data)[[1]], sbtn2))
        tt <- sbt.proba[ , sbtn[3]]
        tt2 <- pprolif
        tt <- cbind(tt * tt2, tt * (1 - tt2))
        colnames(tt) <- sbtn2[3:4]
        sbt.proba2[ , sbtn2[1:2]] <- sbt.proba[ , sbtn[1:2]]
        sbt.proba2[ , sbtn2[3:4]] <- tt[ , sbtn2[3:4]]

        if(plot) {
          if(do.scale) {
            myxlim <- myylim <- c(-2, 2)
          } else {
            myxlim <- range(dd[ , "ESR1"])
            myylim <- range(dd[ , "ERBB2"])
          }
          ## plot the clusters with proliferation
          mycol <- mypch <- rep(NA, length(sbt2))
          mycol[sbt2 == sbtn2[1]] <- "darkred"
          mycol[sbt2 == sbtn2[2]] <- "darkgreen"
          mycol[sbt2 == sbtn2[3]] <- "darkorange"
          mycol[sbt2 == sbtn2[4]] <- "darkviolet"
          mypch[sbt2 == sbtn2[1]] <- 17
          mypch[sbt2 == sbtn2[2]] <- 0
          mypch[sbt2 == sbtn2[3] | sbt2 == sbtn2[4]] <- 10
          mypch <- as.numeric(mypch)
          names(mycol) <- names(mypch) <- names(sbt2)
          plot(x=dd[ , "ESR1"], y=dd[ , "ERBB2"], xlim=myxlim, ylim=myylim, xlab="ESR1", ylab="ERBB2", col=mycol[dimnames(dd)[[1]]], pch=mypch[dimnames(dd)[[1]]])
          legend(x="topleft", col=c("darkred", "darkgreen", "darkorange", "darkviolet"), legend=sbtn2, pch=c(17, 0, 10, 10), bty="n")
        }

        return(list("subtype"=sbt, "subtype.proba"=sbt.proba, "prediction.strength"=ps.res, "BIC"=BIC.res, "subtype2"=sbt2, "subtype.proba2"=sbt.proba2, "prediction.strength2"=ps.res2, "module.scores"=dd2, "mapping"=mymap))
      }
    
    `intrinsic.cluster.predict` <-
      function(sbt.model, data, annot, doMapping=FALSE, mapping, do.prediction.strength=FALSE, verbose=FALSE) {
        
        if(missing(data) || missing(annot) || missing(sbt.model)) { stop("data, annot and sbt.mod parameters must be specified") }
        if (!is.matrix(data)) { data <- as.matrix(data) }
        
        if(is.list(sbt.model)) {
          ## retrieve model
          centroids <- sbt.model$centroids
          annot.centroids <- sbt.model$centroids.map
          method.cor <- sbt.model$method.cor
          method.centroids <- sbt.model$method.centroids
          std <- sbt.model$std
          mq <- sbt.model$rescale.q
          mins <- sbt.model$mins
        } else {
          ## read model file
          ## retrieve centroids
          annot.centroids <- read.csv(sbt.model, comment.char="#", stringsAsFactors=FALSE)
          dimnames(annot.centroids)[[1]] <- annot.centroids[ , "probe"]
          ## retrieve model parameters
          rr <- readLines(con=sbt.model)[-1]
          rr <- rr[sapply(rr, function(x) { return(substr(x=x, start=1, stop=1) == "#") })]
          nn <- unlist(lapply(X=rr, FUN=function(x) { x <- unlist(strsplit(x=unlist(strsplit(x=x, split=":"))[1], split=" ")); x <- x[length(x)]; return(x); }))
          rr2 <- unlist(lapply(X=rr[is.element(nn, c("method.cor", "method.centroids", "std", "rescale.q", "mins"))], FUN=function(x) { x <- unlist(strsplit(x=unlist(strsplit(x=x, split=":"))[2], split=" ")); x <- x[length(x)]; return(x); }))
          method.cor <- rr2[nn == "method.cor"]
          method.centroids <- rr2[nn == "method.centroids"]
          std <- rr2[nn == "std"]
          mq <- as.numeric(rr2[nn == "rescale.q"])
          mins <- as.numeric(rr2[nn == "mins"])
          rr <- rr[!is.element(nn, c("method.cor", "method.centroids", "std", "rescale.q", "mins"))]
          nn <- nn[!is.element(nn, c("method.cor", "method.centroids", "std", "rescale.q", "mins"))]
          cent <- lapply(X=rr, FUN=function(x) { x <- as.numeric(unlist(strsplit(x=unlist(strsplit(x=x, split=":"))[2], split=" "))); x <- x[!is.na(x)]; return(x)})
          centroids <- NULL
          for(i in 1:length(cent)) { centroids <- cbind(centroids, cent[[i]]) }
          nn <- unlist(lapply(X=strsplit(x=nn, split="centroid."), FUN=function(x) { return(x[[2]]) }))
          dimnames(centroids) <- list(dimnames(annot.centroids)[[1]], nn)
        }
        
        number.cluster <- ncol(centroids)
        if(is.null(dimnames(centroids)[[2]])) { name.cluster <- paste("cluster", 1:ncol(centroids), sep=".") } else { name.cluster <- dimnames(centroids)[[2]] }
        
        gt <- nrow(centroids)
        #mapping
        if(doMapping) {
          #mapping through EntrezGene.ID
          #remove (arbitrarily) duplicated gene ids
          centroids.gid <- as.character(annot.centroids[ ,"EntrezGene.ID"])
          names(centroids.gid) <- as.character(annot.centroids[ , "probe"])
          myx <- !duplicated(centroids.gid) & !is.na(centroids.gid)
          centroids.gid <- centroids.gid[myx]
          annot.centroids <- annot.centroids[myx, , drop=FALSE]
          centroids <- centroids[myx, , drop=FALSE]
          gid <- as.character(annot[ ,"EntrezGene.ID"])
          names(gid) <- as.character(annot[ ,"probeset"])
          if(missing(mapping)) { ## select the most variant probes using annotations
            ## if multiple probes for one gene, keep the most variant
            rr <- geneid.map(geneid1=gid, data1=data, geneid2=centroids.gid, verbose=FALSE)
            nn <- match(rr$geneid2, centroids.gid)
            nn <- nn[!is.na(nn)]
            centroids.gid <- centroids.gid[nn]
            annot.centroids <- annot.centroids[nn, ]
            centroids <- centroids[nn, , drop=FALSE]
            data <- rr$data1
          } else { # use a predefined mapping
            nn <- as.character(mapping[ , "EntrezGene.ID"])
            # keep only centroids genes with mapped probes
            myx <- is.element(centroids.gid, nn)
            centroids.gid <- centroids.gid[myx]
            annot.centroids <- annot.centroids[myx, , drop=FALSE]
            centroids <- centroids[myx, , drop=FALSE]
            pp <- as.character(mapping[match(centroids.gid, nn), "probe"])
            myx <- is.element(pp, dimnames(data)[[2]])
            centroids.gid <- centroids.gid[myx]
            annot.centroids <- annot.centroids[myx, , drop=FALSE]
            centroids <- centroids[myx, , drop=FALSE]
            pp <- pp[myx]
            data <- data[ , pp, drop=FALSE]
          }
        }
        else {
          if(all(!is.element(dimnames(data)[[2]], dimnames(centroids)[[1]]))) { stop("no probe in common -> annot or mapping parameters are necessary for the mapping process!") }
          ## no mapping are necessary
          myx <- intersect(dimnames(centroids)[[1]], dimnames(data)[[2]])
          data <- data[ ,myx, drop=FALSE]
          centroids <- centroids[myx, , drop=FALSE]
        }
        centroids.map <- cbind("probe"=dimnames(data)[[2]], "probe.centroids"=dimnames(centroids)[[1]], "EntrezGene.ID"=as.character(annot[dimnames(data)[[2]], "EntrezGene.ID"]))
        dimnames(centroids.map)[[1]] <- dimnames(data)[[2]]
        gm <- nrow(centroids)
        if(gm == 0 || (sum(is.na(data)) / length(data)) > 0.9) { ## no mapping or too many missing values
          ncl <- rep(NA, nrow(data))
          names(ncl) <- dimnames(data)[[1]]
          nproba <- ncor <- matrix(NA, nrow=nrow(data), ncol=ncol(centroids), dimnames=list(dimnames(data)[[1]], name.cluster))
          ps.res <- NULL
          if(do.prediction.strength) { ps.res <- list("ps"=NA, "ps.cluster"=ncor[ , 1], "ps.individual"=ncl) }
          tt <- matrix(NA, ncol=nrow(centroids.map), nrow=nrow(data), dimnames=list(dimnames(data)[[1]], dimnames(centroids.map)[[1]]))
          return(list("subtype"=ncl, "subtype.proba"=nproba, "cor"=ncor, "prediction.strength"=ps.res, "centroids.map"=centroids.map, "profiles"=tt))
        }
        if(verbose) { message(sprintf("%i/%i probes are used for clustering", gm, gt)) }
        
        #standardization of the gene expressions
        switch(std,
               "scale"={
                 data <- scale(data, center=TRUE, scale=TRUE)
                 if(verbose) { message("standardization of the gene expressions") }
               }, 
               "robust"={
                 data <- apply(data, 2, function(x) { return((rescale(x, q=mq, na.rm=TRUE) - 0.5) * 2) })
                 if(verbose) { message("robust standardization of the gene expressions") }
               }, 
               "none"={ if(verbose) { message("no standardization of the gene expressions") } })
        
        ## apply the nearest centroid classifier to classify the samples again
        ncor <- t(apply(X=data, MARGIN=1, FUN=function(x, y, method.cor) { 
          rr <- array(NA, dim=ncol(y), dimnames=list(colnames(y)))
          if (sum(complete.cases(x, y)) > 3) {
            rr <- cor(x=x, y=y, method=method.cor, use="complete.obs")
          }
          return (rr)
        }, y=centroids, method.cor=method.cor))
        #nproba <- t(apply(X=ncor, MARGIN=1, FUN=function(x) { return(abs(x) / sum(abs(x), na.rm=TRUE)) }))
        ## negative correlations are truncated to zero since they have no meaning for subtypes identification
        nproba <- t(apply(X=ncor, MARGIN=1, FUN=function (x) {
          rr <- array(NA, dim=length(x), dimnames=list(names(x)))
          x[!is.na(x) & x < 0] <- 0
          if (!all(is.na(x))) {
            rr <- x / sum(x, na.rm=TRUE)
          }
          return (rr)
        }))
        dimnames(ncor) <- dimnames(nproba) <- list(dimnames(data)[[1]], name.cluster)
        ncl <- apply(X=ncor, MARGIN=1, FUN=function(x) { return(order(x, decreasing=TRUE, na.last=TRUE)[1]) })
        names(ncl) <- dimnames(data)[[1]]
        ## names of identified clusters
        ncln <- name.cluster[ncl]
        names(ncln) <- dimnames(data)[[1]]
        
        ## if one or more subtypes have not been identified, remove them for prediction strength
        myx <- sort(unique(ncl))
        myx <- myx[!is.na(myx)]
        name.cluster2 <- name.cluster[myx]
        number.cluster2 <- length(myx)
        
        ps.res <- ncl2 <- NULL
        if(do.prediction.strength) {
          ## compute the clustering and cut the dendrogram
          ## hierarchical clustering with correlation-based distance and average linkage
          hcl <- amap::hcluster(x=data, method="correlation", link="average")
          mins.ok <- stop.ok <- FALSE
          nbc <- number.cluster2
          nclust.best <- 1
          while(!mins.ok && !stop.ok) { ## until each cluster contains at least mins samples
            cl <- cutree(tree=hcl, k=nbc)
            tt <- table(cl)
            if(sum(tt >= mins) >= number.cluster2) {
              if(nbc > number.cluster2) { ## put NA for clusters with less than mins samples
                td <- names(tt)[tt < mins]
                cl[is.element(cl, td)] <- NA
                ## rename the clusters
                ucl <- sort(unique(cl))
                ucl <- ucl[!is.na(ucl)]
                cl2 <- cl
                for(i in 1:number.cluster2) { cl2[cl == ucl[i] & !is.na(cl)] <- i }
                cl <- cl2
              }
              nclust.best <- number.cluster2
              mins.ok <- TRUE
            } else {
              if(sum(tt >= mins) > nclust.best) {
                nbc.best <- nbc
                nclust.best <- sum(tt >= mins)
              }
              nbc <- nbc + 1
              if(nbc > (nrow(data) - (number.cluster2 * mins))) {
                warning(sprintf("impossible to find %i main clusters with at least %i individuals!", number.cluster2, mins))
                stop.ok <- TRUE
              }
            }
            if(stop.ok) { ## no convergence for the clustering with mininmum set of individuals
              cl <- cutree(tree=hcl, k=nbc.best)
              tt <- table(cl)
              td <- names(tt)[tt < mins]
              cl[is.element(cl, td)] <- NA
              ## rename the clusters
              ucl <- sort(unique(cl))
              ucl <- ucl[!is.na(ucl)]
              cl2 <- cl
              for(i in 1:nclust.best) { cl2[cl == ucl[i] & !is.na(cl)] <- i }
              cl <- cl2
            }
          }
          ## compute the centroids
          ## take the core samples in each cluster to compute the centroid
          ## not feasible due to low intra correlation within clusters!!!
          ## minimal pairwise cor of approx 0.3
          #cl2 <- cutree(tree=hcl, h=0.7)
          #table(cl, cl2) to detect which core cluster of samples for which cluster.
          cl.centroids <- matrix(NA, nrow=ncol(data), ncol=nclust.best, dimnames=list(dimnames(data)[[2]], paste("cluster", 1:nclust.best, sep=".")))
          for(i in 1:ncol(cl.centroids)) {
            switch(method.centroids, 
                   "mean"={ cl.centroids[ ,i] <- apply(X=data[cl == i & !is.na(cl), ,drop=FALSE], MARGIN=2, FUN=mean, na.rm=TRUE, trim=0.025) }, 
                   "median"={ cl.centroids[ ,i] <- apply(X=data[cl == i & !is.na(cl), ,drop=FALSE], MARGIN=2, FUN=median, na.rm=TRUE) }, 
                   "tukey"={ cl.centroids[ ,i] <- apply(X=data[cl == i & !is.na(cl), ,drop=FALSE], MARGIN=2, FUN=tbrm, na.rm=TRUE, C=9) })
          }
          #apply the nearest centroid classifier to classify the samples again
          ncor2 <- t(apply(X=data, MARGIN=1, FUN=function(x, y, z) { return(cor(x, y, method=z, use="complete.obs")) }, y=cl.centroids, z=method.cor))
          nproba2 <- t(apply(X=ncor2, MARGIN=1, FUN=function(x) { return(abs(x) / sum(abs(x), na.rm=TRUE)) }))
          dimnames(ncor2) <- dimnames(nproba2) <- list(dimnames(data)[[1]], dimnames(cl.centroids)[[2]])
          ncl2 <- apply(X=ncor2, MARGIN=1, FUN=function(x) { return(order(x, decreasing=TRUE)[1]) })
          names(ncl2) <- dimnames(data)[[1]]
          ## rename clusters since we do not expect to get the same id per cluster
          ## this avoids a warning in ps.cluster
          uncl <- sort(unique(ncl))
          uncl <- uncl[!is.na(uncl)]
          nclt <- ncl
          for(mm in 1:length(uncl)) {
            nclt[ncl == uncl[mm]] <- mm
          }
          uncl2 <- sort(unique(ncl2))
          uncl2 <- uncl2[!is.na(uncl2)]
          ncl2t <- ncl2
          for(mm in 1:length(uncl2)) {
            ncl2t[ncl2 == uncl2[mm]] <- mm
          }
          #prediction strength
          ps.res <- ps.cluster(cl.tr=ncl2t, cl.ts=nclt, na.rm=TRUE)
          ## put NA for clusters which are potentially not present in the dataset
          tt <- rep(NA, length(name.cluster))
          names(tt) <- name.cluster
          tt[name.cluster2] <- ps.res$ps.cluster
          ps.res$ps.cluster <- tt
        }
        
        
        return(list("subtype"=ncln, "subtype.proba"=nproba, "cor"=ncor, "prediction.strength"=ps.res, "subtype.train"=ncl2, "profiles"=data, "centroids.map"=centroids.map))
      }

    subtypeModel <- match.arg(subtypeModel)

    ## convert SCM to SSP nomenclature
    sbt.conv <- rbind(c("ER-/HER2-", "Basal"),
                      c("HER2+", "Her2"),
                      c("ER+/HER2- High Prolif", "LumB"),
                      c("ER+/HER2- Low Prolif", "LumA")
    )
    colnames(sbt.conv) <- c("SCM.nomenclature", "SSP.nomenclature")

    sbtn.ssp <- c("Basal", "Her2", "LumB", "LumA", "Normal")
    sbtn2.ssp <- c("Basal", "Her2", "Lums", "LumB", "LumA", "Normal")

    ## SCM family
    if (subtypeModel %in% c("scmgene", "scmod1", "scmod2")) {
      switch(subtypeModel,
             "scmgene" = {
               sbts <- subtype.cluster.predict(subtypeModel=scmgene.robust, data=data, annot=annot, doMapping=doMapping)[c("subtype2", "subtype.proba2")]
             },
             "scmod1" = {
               sbts <- subtype.cluster.predict(subtypeModel=scmod1.robust, data=data, annot=annot, doMapping=doMapping)[c("subtype2", "subtype.proba2")]
             },
             "scmod2" = {
               sbts <- subtype.cluster.predict(subtypeModel=scmod2.robust, data=data, annot=annot, doMapping=doMapping)[c("subtype2", "subtype.proba2")]
             }
      )
      names(sbts) <- c("subtype", "subtype.proba")
      ## compute crisp classification
      sbts$subtype.crisp <- t(apply(sbts$subtype.proba, 1, function (x) {
        xx <- array(0, dim=length(x), dimnames=list(names(x)))
        xx[which.max(x)] <- 1
        return (xx)
      }))

      ## reorder columns
      #ss <- sbtn2.ssp[is.element(sbtn2.ssp, colnames(sbts$subtype.proba))]
      #sbts$subtype.proba <- sbts$subtype.proba[ , ss, drop=FALSE]
      #sbts$subtype.crisp <- sbts$subtype.crisp[ , ss, drop=FALSE]

      ## set the proper names
      names(sbts$subtype) <- rownames(sbts$subtype.proba) <- rownames(sbts$subtype.crisp)<- rownames(data)
    }

    ## SSP family
    if (subtypeModel %in% c("ssp2003", "ssp2006", "pam50")) {
      switch(subtypeModel,
             "pam50" = {
               sbts <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=data, annot=annot, doMapping=doMapping)[c("subtype", "subtype.proba")]
             },
             "ssp2006" = {
               sbts <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=data, annot=annot, doMapping=doMapping)[c("subtype", "subtype.proba")]
             },
             "ssp2003" = {
               sbts <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=data, annot=annot, doMapping=doMapping)[c("subtype", "subtype.proba")]
             }
      )
      sbts$subtype <- factor(as.character(sbts$subtype), levels=sbtn.ssp)
      ## compute crisp classification
      sbts$subtype.crisp <- t(apply(sbts$subtype.proba, 1, function (x) {
        xx <- array(0, dim=length(x), dimnames=list(names(x)))
        xx[which.max(x)] <- 1
        return (xx)
      }))

      ## merge LumA and LumB: #sum the probability for LumA and LumB to get the probability for Luminals in general
      #lums.proba <- apply(sbts$subtype.proba[ , c("LumB", "LumA"), drop=FALSE], 1, sum, na.rm=TRUE)
      #sbts$subtype.proba <- cbind(sbts$subtype.proba, "Lums"=lums.proba)
      #lums.crisp <- as.numeric(is.element(sbts$subtype, c("LumA", "LumB")))
      #sbts$subtype.crisp <- cbind(sbts$subtype.crisp, "Lums"=lums.crisp)

      ## reorder columns
      #ss <- sbtn2.ssp[is.element(sbtn2.ssp, colnames(sbts$subtype.proba))]
      #sbts$subtype.proba <- sbts$subtype.proba[ , ss, drop=FALSE]
      #sbts$subtype.crisp <- sbts$subtype.crisp[ , ss, drop=FALSE]

      ## set the proper names
      names(sbts$subtype) <- rownames(sbts$subtype.proba) <- rownames(sbts$subtype.crisp)<- rownames(data)
    }

    ## IntClust family
    if (subtypeModel %in% c("intclust")) {
      #message("Note: Need a Gene.Symbol column in the annotation object")
      sbts<-NULL
      myx <- !is.na(annot[ , "gene"]) & !duplicated(annot[ , "gene"])
      dd <- t(data[ , myx, drop=FALSE])
      rownames(dd) <- annot[myx, "gene"]
      ## remove patients with more than 80% missing values
      rix <- apply(dd, 2, function (x, y) { return ((sum(is.na(x) / length(x))) > y) }, y=0.8)
      cix <- apply(dd, 2, function (x, y) { return ((sum(is.na(x) / length(x))) > y) }, y=0.8)
      dd <- dd[!rix, !cix, drop=FALSE]
      features <- iC10::matchFeatures(Exp=dd, Exp.by.feat="gene")
      features <- iC10::normalizeFeatures(features, method="scale")
      res <- iC10::iC10(features)
      ## compute crisp classification
      crisp <- t(apply(res$posterior, 1, function (x) {
        xx <- array(0, dim=length(x), dimnames=list(names(x)))
        xx[which.max(x)] <- 1
        return (xx)
      }))
      sbts$subtype <- array(NA, dim=nrow(data), dimnames=list(rownames(data)))
      sbts$subtype[!rix] <- res$class
      sbts$subtype.proba <- array(NA, dim=c(nrow(data), ncol(res$posterior)), dimnames=list(rownames(data), colnames(res$posterior)))
      sbts$subtype.proba[!rix, ] <- res$posterior
      sbts$subtype.crisp <- t(apply(sbts$subtype.proba, 1, function (x) {
        xx <- array(0, dim=length(x), dimnames=list(names(x)))
        xx[which.max(x)] <- 1
        return (xx)
      }))
      ## set the proper colnames
      colnames(sbts$subtype.proba) <- colnames(sbts$subtype.crisp) <- paste("iC", colnames(sbts$subtype.proba), sep="")
      sbts$subtype <- paste("iC", sbts$subtype, sep="")
      sbts$subtype <- factor(sbts$subtype, levels=colnames(sbts$subtype.proba))
      ## set the proper rownames
      names(sbts$subtype) <- rownames(sbts$subtype.proba) <- rownames(sbts$subtype.crisp)<- rownames(data)
    }

    ## AIMS classifier
    if (subtypeModel %in% c("aims")) {
      sbts <- AIMS::applyAIMS(eset=t(data), EntrezID=annot[ , "EntrezGene.ID"])[c("cl", "all.probs")]
      sbts$subtype <- sbts$cl
      sbts$subtype.proba <- matrix(unlist(sbts$all.probs$`20`), ncol = 5, byrow = TRUE)
      colnames(sbts$subtype.proba) <- colnames(sbts$all.probs$`20`)
      rownames(sbts$subtype.proba) <- rownames(sbts$subtype)

      ## compute crisp classification
      sbts$subtype.crisp <- t(
        apply(sbts$subtype.proba, 1, function (x) {
          xx <- array(0, dim=length(x), dimnames=list(names(x)))
          xx[which.max(x)] <- 1
          return (xx)
        })
      )
      sbts<-sbts[- which(names(sbts) %in% c("cl","all.probs"))]
    }

    ## CLAUDIN-LOW classifier
    if (subtypeModel %in% c("claudinlow")) {
      train<-claudinLowData
      train$xd<- medianCtr(train$xd)

      if(doMapping) {
        gid1 <- as.numeric(rownames(train$xd))
        names(gid1) <- rownames(train$xd)
        gid2 <- as.numeric(as.character(annot[ ,"EntrezGene.ID"]))
        names(gid2) <- dimnames(annot)[[1]]

        ## remove missing and duplicated geneids from the gene list
        rm.ix <- is.na(gid1) | duplicated(gid1)
        gid1 <- gid1[!rm.ix]
        rr <- geneid.map(geneid1=gid2, data1=data, geneid2=gid1, verbose=FALSE)
        gt <- length(rr$geneid2)
        if(is.na(rr$geneid1[1])) {
          gm <- 0
          #no gene ids in common
          res <- rep(NA, nrow(data))
          names(res) <- dimnames(data)[[1]]
          gf <- c("mapped"=0, "total"=gt)
          if(verbose) { message(sprintf("probe candidates: 0/%i", gt)) }
          return(list("score"=res, "risk"=res, "mapping"=gf, "probe"=NA))
        }
        gid1 <- rr$geneid2
        gid2 <- rr$geneid1
        data <- rr$data1
        #mymapping <- c("mapped"=gm, "total"=gt)
        myprobe <- cbind("probe"=names(gid1), "EntrezGene.ID"=gid1, "new.probe"=names(gid2))
        ## change the names of probes in the data
        dimnames(data)[[2]] <- names(gid2) <- names(gid1)
      }

      test<- medianCtr(t(data)) #probes as rows, median-centered
      
      claudinLow <- function(x, classes="", y, nGenes="", priors="equal", std=FALSE, distm="euclidean", centroids=FALSE){
        
        dataMatrix <- x
        features <- dim(x)[1]
        samples <- dim(x)[2]
        sampleNames <- dimnames(x)[[2]]
        featureNames <- dimnames(x)[[1]]
        
        #parse the test file - same as train file but no rows of classes
        tdataMatrix <- y
        tfeatures <- dim(y)[1]
        tsamples <- dim(y)[2]
        tsampleNames <- dimnames(y)[[2]]
        tfeatureNames <- dimnames(y)[[1]]
        
        overlapSets<-function(x,y){
          
          # subset the two lists to have a commonly ordered gene list
          x<-x[dimnames(x)[[1]] %in% dimnames(y)[[1]],]
          y<-y[dimnames(y)[[1]] %in% dimnames(x)[[1]],]
          
          #and sort such that thing are in the correct order
          x<-x[sort.list(row.names(x)),]
          y<-y[sort.list(row.names(y)),]
          
          return(list(x=x,y=y))
        }
        
        #dimnames(tdataMatrix)[[2]] <- paste("x",seq(1,471))
        temp <- overlapSets(dataMatrix,tdataMatrix)
        dataMatrix <- temp$x
        tdataMatrix <- temp$y
        sfeatureNames <- row.names(dataMatrix)
        
        # standardize both sets
        if(std){
          dataMatrix <- standardize(dataMatrix)
          tdataMatrix <- standardize(tdataMatrix)
        }
        
        if(!centroids){
          thisClass <- as.vector(classes[,1])
          nClasses <- nlevels(as.factor(thisClass))
          classLevels <- levels(as.factor(thisClass))
          for(j in 1:nClasses){
            thisClass[thisClass==classLevels[j]] <- j
          }
          thisClass <- as.numeric(thisClass)
          dataMatrix <- dataMatrix[,!(is.na(thisClass))]
          thisClass <- thisClass[!(is.na(thisClass))]
          
          scores <- apply(dataMatrix,1,bwss,thisClass)
          trainscores <- vector()	
          for(j in 1:dim(dataMatrix)[1]){			
            trainscores[j] <- scores[[row.names(dataMatrix)[j]]]$bss / scores[[row.names(dataMatrix)[j]]]$wss
          }
          
          dataMatrix <- dataMatrix[sort.list(trainscores,decreasing=T),]
          tdataMatrix <- tdataMatrix[sort.list(trainscores,decreasing=T),]	
          
          if(nGenes==""){
            nGenes <- dim(dataMatrix)[1]
          }
          print(paste("Number of genes used:",nGenes))
          
          dataMatrix <- dataMatrix[1:nGenes,]
          tdataMatrix <- tdataMatrix[1:nGenes,]
          
          centroids <- matrix(,nrow=nGenes,ncol=nClasses)
          for(j in 1:nClasses){
            centroids[,j] <- apply(dataMatrix[,thisClass==j],1,mean)
          }
          dimnames(centroids) <- list(row.names(dataMatrix),NULL)
          
        }else{
          nGenes <- dim(dataMatrix)[1]
          print(paste("Number of genes used:",nGenes))
          centroids <- dataMatrix
          nClasses <- dim(centroids)[2]
          classLevels <- dimnames(centroids)[[2]]
        }
        
        distances <- matrix(ncol=nClasses,nrow=dim(tdataMatrix)[2])
        for(j in 1:nClasses){
          if(distm=="euclidean"){
            distances[,j] <- dist(t(cbind(centroids[,j],tdataMatrix)))[1:(dim(tdataMatrix)[2])]
          }
          if(distm=="correlation" | distm=="pearson"){
            distances[,j] <- t(-1*cor(cbind(centroids[,j],tdataMatrix),use="pairwise.complete.obs"))[2:(dim(tdataMatrix)[2]+1)]
          }
          if(distm=="spearman"){
            distances[,j] <- t(-1*cor(cbind(centroids[,j],tdataMatrix),method="spearman",use="pairwise.complete.obs"))[2:(dim(tdataMatrix)[2]+1)]
          }
          colnames(distances) <- c("euclidian distance to Claudin-low", "euclidian distance to Others")
          rownames(distances) <- tsampleNames
          
        }
        
        scores <- apply(distances,1,min)
        prediction <- vector(length=tsamples)
        for(i in 1:tsamples){
          prediction[i] <- classLevels[match(scores[i],distances[i,])]
        }
        names(prediction) <- tsampleNames
        prediction <- data.frame(Samples=tsampleNames, prediction)
        colnames(prediction) <- c("Samples", "Call")
        
        return(list(predictions=prediction,testData=tdataMatrix,distances=distances,centroids=centroids))
      }
      
      #Run Classifier Call
      #is returning everyone to be others and none to be claudinLow subtype.
      #likely due to euclidean distance smalled for group others always, look into later
      predout<-claudinLow(train$xd, as.matrix(train$classes$Group,ncol=1), test)
      sbts<-NULL
      sbts$subtype<-factor(as.character(predout$predictions$Call))
      colnames(predout$centroids)<-c("Claudin","Others")

      ## apply the nearest centroid classifier to classify the samples again
      ncor <- t(apply(X=data, MARGIN=1, FUN=function(x, y) {
        rr <- array(NA, dim=ncol(y), dimnames=list(colnames(y)))
        if (sum(complete.cases(x, y)) > 3) {
          rr <- cor(x=x, y=y, method="spearman", use="complete.obs")
        }
        return (rr)
      }, y=predout$centroids))

      #Calculate posterior probability based on the correlationss
      # nproba <- t(apply(X=ncor, MARGIN=1, FUN=function(x) { return(abs(x) / sum(abs(x), na.rm=TRUE)) }))

      # negative correlations are truncated to zero since they have no meaning for subtypes identification
      nproba <- t(apply(X=ncor, MARGIN=1, FUN=function (x) {
        rr <- array(NA, dim=length(x), dimnames=list(names(x)))
        x[!is.na(x) & x < 0] <- 0
        if (!all(is.na(x))) {
          rr <- x / sum(x, na.rm=TRUE)
        }
        return (rr)
      }))

      sbts$subtype.proba<-nproba

      ## compute crisp classification - in this case, really based on the binary call from the CL classifier
      #     sbts$subtype.crisp <- t(
      #       apply(sbts$subtype.proba, 1, function (x) {
      #         xx <- array(0, dim=length(x), dimnames=list(names(x)))
      #         xx[which.max(x)] <- 1
      #         return (xx)
      #       })
      #     )
      #     colnames(sbts$subtype.crisp)<-c("Claudin","Others")

      # In this case, really based on the binary call from the CL classifier. Use that for accuracy
      CLsubtypes<-c("Claudin","Others")
      sbts$subtype.crisp <- matrix(0, nrow=nrow(predout$predictions), ncol=2,dimnames=list(rownames(predout$predictions),CLsubtypes))
      for(count in 1:nrow(predout$predictions))
      {
        if(predout$predictions$Call[count]=="Others")
          sbts$subtype.crisp[count,2]<-1
        else sbts$subtype.crisp[count,1]<-1
      }
    }
    return (sbts)
  }
