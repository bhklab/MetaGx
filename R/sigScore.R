`sigScore` <- 
function (eset, sig, method=c("principal.component", "weighted.average"), scaling=c("none", "std", "robust")) {
  #Take a specify gene set and create a single metagene base it's the selected method
  #
  # Args:
  #   eset: expression set object
  #   sig: vector of probe/gene identifiers or data.frame with two columns reporting the probe/gene identifier (character) and coefficiwents/weigths (numeric), respectively
  #   method: summarization method for signature score computation, either PCA or simple weighted average
  #   scaling: shoudl the signature score be scaled, std stands for standard z score and robust use quantile-based scaline
  #
  # Returns:     
  #     Signature scores
  
  method <- match.arg(method)
  scaling <- match.arg(scaling)
  if (!is.vector(sig) && !is.data.frame(sig)) {
    stop ("sig must be a vector or a data frame")
  }
  if (is.vector(sig)) {
    sig <- data.frame(as.character(sig), rep(1, length(sig)))
  }
  sig <- sig[ , c(1, 2), drop=FALSE]
  dimnames(sig) <- list(sig[ , 1], c("feature", "coefficient"))
  signn <- nrow(sig)
  sig <- sig[complete.cases(sig), , drop=FALSE]
  if (nrow(sig) == 0) {
    stop ("All values are missing in the signature")
  }
  sig <- sig[intersect(sig[ , "feature"], rownames(Biobase::fData(eset))), , drop=FALSE]
  if (nrow(sig) == 0) {
    warning ("Feature identifers from the signature do not map with those of the expressionSet object")
    sigscore <- array(NA, dim=nrow(Biobase::pData(eset)), dimnames=list(rownames(Biobase::pData(eset))))
    return (sigscore)
  }
  if (nrow(sig) < signn) {
    warning(sprintf("%i/%i genes were present in the expressionSet object", nrow(sig), signn))
  }

  if (nrow(sig) > 1) {
    switch(method,
     "principal.component" = {
       sigscore <- array(NA, dim=nrow(Biobase::pData(eset)), dimnames=list(rownames(Biobase::pData(eset))))
       if(nrow(sig) >= 3) {
         ## compute the first principal component
         rr <- Biobase::exprs(eset)[sig[ , "feature"], , drop = FALSE]
         rr <- rr[ , complete.cases(t(rr)), drop=FALSE]
         prm <- prcomp(x=t(rr), na.action="na.omit")
         sigscore[rownames(prm$x)] <- prm$x[ , 1]
         if (lsa::cosine(x=prm$rotation[ , 1], y=sig[ , "coefficient"]) < 0) {
           ## inversion of the "direction" of the principal component
           sigscore <- - sigscore
         }        
       } else {
         stop (sprintf("Not enough genes (%i) in the signature to perform principal component analsysis", nrow(sig)))
       }
     },
     "weighted.average" = {
        sigp <- sig[sig[ , "coefficient"] >= 0, , drop=FALSE]
        if (nrow(sigp) > 0) {
          ssp <- apply(X=Biobase::exprs(eset)[sigp[ , "feature"], , drop=FALSE], MARGIN=2, FUN=function (x, y) {
            nix <- is.na(x)
            return(sum(x * y, na.rm=TRUE) / sum(y[!nix]))
          }, y=abs(sigp[ , "coefficient"]))
        } else {
          ssp <- 0
        }
        sign <- sig[sig[ , "coefficient"] < 0, , drop=FALSE]
        if (nrow(sign) > 0) {
          ssn <- apply(X=Biobase::exprs(eset)[sign[ , "feature"], , drop=FALSE], MARGIN=2, FUN=function (x, y) {
            nix <- is.na(x)
            return(sum(x * y, na.rm=TRUE) / sum(y[!nix]))
          }, y=abs(sign[ , "coefficient"]))
        } else {
          ssn <- 0
        }
        sigcr <- sum(abs(sigp[ , "coefficient"]), na.rm=TRUE) / sum(abs(sig[ , "coefficient"]), na.rm=TRUE)
        sigscore <- (sigcr * ssp - (1 - sigcr) * ssn) / 2
      }
    )
  } else {
    sigscore <- Biobase::exprs(eset)[which(is.element(rownames(Biobase::fData(eset)), sig[1, 1]))[1], ] * sig[1, 2]
  }
  ## scaling
  switch(scaling,
    "std" = {
      sigscore <- scale(sigscore, center=TRUE, scale=TRUE)
    },
    "robust" = {
      sigscore <- genefu::rescale(x=sigscore, q=0.05, na.rm=TRUE)
    },
    "none" = {
    }
  )
  
  return(sigscore)
}


## end

