#########
## Natchar Ratanasirigulchai
## March 6, 2015
#########

forestplot <- function(esets, surv.type, entrezgene, formula=y~entrezgene,
                       mlab="Overall", rma.method="FE", at=NULL,xlab="Hazard Ratio",...) {
  entrezgene <-  paste("geneid", as.character(entrezgene ), sep=".") 
  require(metafor)
   esets <- esets[sapply(esets, function(x) entrezgene %in% featureNames(x))]
   coefs <- sapply(1:length(esets), function(i) {

     switch(surv.type,
            "os" = {
              stime <- pData(esets[[i]])[, "days_to_death"]
            },
            "rfs" = {
              stime <- pData(esets[[i]]) [,"days_to_tumor_recurrence"]
            }, 
            "dmfs"= {
              stime <- pData(esets[[i]]) [, "t.dmfs"]
            })

     #   sevent <- Biobase::pData(eset)[ , sprintf("e.%s", surv.type)]
     switch(surv.type,
            "os" = {
              sevent <-as.numeric(pData(esets[[i]])[, "vital_status"] == "deceased")
            },
            "rfs" = {
              sevent <- pData(esets[[i]]) [,"recurrence_status"]
            },
            "dmfs" = {
              sevent <- pData(esets[[i]])[, "e.dmfs"]
            })
     
     tmp <- as(phenoData(esets[[i]]), "data.frame")
#      tmp$y <- esets[[i]][[y]]
     tmp$y <- Surv(stime, sevent)
     tmp$entrezgene <- exprs(esets[[i]])[entrezgene,]
    
       summary(coxph(formula,data=tmp))$coefficients[1,c(1,3)]
     })
     res.rma <- metafor::rma(yi = coefs[1,], sei = coefs[2,],
                               method=rma.method)
      if (is.null(at)) at <- log(c(0.25,1,4,20))
    forest.rma(res.rma, xlab=xlab, slab=gsub("_eset$","",names(esets)),
                 atransf=exp, at=at, mlab=mlab,...)
    return(res.rma)
    }
