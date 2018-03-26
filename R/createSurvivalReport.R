
#' Function that assesses the prognostic value of a gene signature and returns a pdf of the results
#'
#' This function runs a survival analysis on breast or ovarian cancer patients for a gene signature using overall survival or relapse freee survival data
#' @param geneSigList A list of character vectors containing the ensemble IDs, entrez IDs, or gene symbols for the genes signatures/genes to conduct a survival analysis on. Entrez IDs are
#'  recommended as if gene symbols or ensemble Ids cant be mapped to entrez Ids then the gene will be omitted from the analysis. Note that The names of the
#' elements in the list will correspond to the names of the signatures in the report, and default names will be provided if the list elements are not named.
#' @param geneDirecList A list of numeric vectors composed of +1 and - 1 indicating the direction of association for each vector of genes supplied in the geneSigList vectors.
#' +1 implies that the expression of that gene in a patient will be added to the patients score and -1 implies the expression will be subtracted from their score.
#' If one is looking for high scores to be associated with good survival, than genes that high expression is believed to lead to good prognosis should be given +1 and genes that high expression is believed
#' to lead to bad prognosis should be given a minus 1.
#' @param cancerType A string representing the type of cancer one would like to test the gene signature on. Options are currently "ovarian" and "breast"
#' @param subtype a string representing the subtyping scheme the patients will be classified under. For breast cancer, the options are "scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003", "intClust", "AIMS", or "claudinLow".
#' For ovarian cancer the options are "verhaak", "bentink", "tothill", "helland", and "konecny".
#' @param survivalMetric a string specifying what type of survival data to use in the analysis. If the cancerType is set to "ovarian"
#' then enter "overall" for overall survival and "relapse" for relapse free survival. If the cancerType is set to "breast" then the options
#' are "relapse" for relapse free survival, "overallMeta" for overall survival on the metabric study data, "overallTcga" for
#' overall survival on the TCGA study data, "overall" for overall survival on the remaining studies (No Metabric and TCGA data), and metastasis for distant metastasis free survival.
#' Note that for all cancer types the option "hierarchy" is available and this survivalMetric parameter checks for relapse free, then distant metastasis, and then overall survival data
#' in each dataset and keeps the first survival data it finds. Note that the ovarian datasets do not have distant metastasis free survival info and for the breast analyses metabric and tcga are removed.
#' @param numGroups an integer specifying the number of groups the patients will be split into when generating survival curves and calculating the log rank p value. The value is 2 by default and cannot exceed 5.
#' @param dataNames a character vector specifying the names of the datasets to include in the analysis. The names must match the names in the column "Data Name" from the dataframe returned from the
#' obtainDataInfo function. Use loadMetaData followed by obtainDataInfo with the cancerType and survivalMetric of interest to get the appropriate data names in the obtainDataInfo table for your analysis. By default
#' all datasets are included.
#' @param patientNames a character vector specifying which patients from the datasets should dbe used in the analysis. The names should corresponsd with the IDs of the patients from the esets obtained via the loadMetaData function
#' @param patientScoresList a list of numeric vectors specifying the patient scores in each signature for the patients provided in the patientNames variables
#' @param soloGeneAnalysis a boolean specifying whether each unique gene in geneSigList should have an survival analysis conducted on it in order to assess its
#' prognostic value independent of its gene signature. Results will show up in a section called Individual Gene Survival Analysis and high scores in the survival curves correspond to high expression as the direction is defaulted to 1
#' to allow for easy comparison amongst all the genes. Default value is FALSE
#' @param removeMid a number greater than 0 and less than .5 specifying the fraction of the patients with scores/risk predictions in the middle of the patients scores to be removed when generating survival results
#' The default is 0, all patients used in the analysis
#' @param censorTime a number specifying the point in time (years) at which the survival data must be censored. The default is 10 years
#' @param addBenchmarks a boolean specifying whether to add known signatures from literature to the analysis for comparison to the provided signatures. Default value is FALSE
#' @param genesRequired a fraction between 0 and 1 (1 inclusive) specifying what fraction of genes from a signature must be present in a dataset for the patients in that dataset to be used in the analysis of that signature (default is 0, use all datasets)
#' @param docTitle a title for the output document with the results. The default is metaGxReport.pdf
#' @return A pdf called MetaGxReport in R's working directory at the time the function was called
#' @importFrom matrixStats rowIQRs iqr
#' @importFrom XML htmlParse xpathApply xmlGetAttr
#' @importFrom knitr knit2pdf
#' @importFrom gplots heatmap.2
#' @importFrom httr GET
#' @importFrom GSVA gsva
#' @importFrom ranger ranger
#' @importFrom survcomp censor.time combine.est D.index km.coxph.plot
#' @importFrom survival Surv
#' @importFrom forestplot forestplot fpTxtGp fpColors
#' @importFrom grid gpar unit
#' @export
#' @examples
#' geneSigList = list(c("KLK14","RHOX8","ADAMTS20","IDO1"), c("ADAMTS20","IDO1", "FAP"))
#' names(geneSigList) = c("4 Gene Signature", "3 Gene Signature")
#' geneDirecList = list(c(1, 1, -1, -1), c(-1, -1, 1))
#' createSurvivalReport(geneSigList, geneDirecList, cancerType = "breast", subtype = "scmod2", survivalMetric = "relapse", addBenchmarks = TRUE)
#'

createSurvivalReport = function(geneSigList, geneDirecList, cancerType, subtype, survivalMetric, numGroups = 2, dataNames = NULL, patientNames = NULL, patientScores = NULL, soloGeneAnalysis = FALSE, removeMid = 0, censorTime = 10, addBenchmarks = FALSE, genesRequired = 0, docTitle = "metaGxReport.pdf")
{
  #install.packages('knitr', dependencies = TRUE)
  #Sys.setenv(JAVA_HOME='C:\\Users\\micha\\Documents\\jre1.8.0_111')
  #library(metaGx)
  #library(gdata)
  #library(forestplot)
  #library(xlsx)
  #library(GSVA)
  #library(mclust)
  #library(org.Hs.eg.db)
  #library(hgu133plus2.db)
  #library(HiDimDA)
  #library(matrixStats)
  #library(survival)
  #library(survcomp)
  #library(MetaGxOvarian)
  #library(MetaGxBreast)
  #library(gplots)
  #library(knitr)
  #library(XML)
  #library(httr)
  #library(Biobase)
  #library(ranger)
  #below 3 needed for breast subtype code
  #library(AIMS)
  #library(iC10)
  #library(limma)
  #source(system.file("extdata", "patientselection.config", package="MetaGxOvarian"))
  #source(system.file("extdata", "createEsetList.R", package="MetaGxOvarian"))

  #if(mord == FALSE)
  #  setwd("C:\\Users\\micha\\Documents\\PMH Research\\Yale Alberto Meta Analysis\\ProjectTest\\Sweave Project MetaGx Fresh\\R Code")
  #source("testGenesPrognostics.R")
  #source("getSurvPlots.R")
  # @importFrom MetaGxOvarian createEsetList

  if(is.list(geneSigList) == FALSE)
    stop("The geneSigList variable is not a list")

  if(is.list(geneDirecList) == FALSE)
    stop("The geneDirecList variable is not a list")

  #geneFrame = metaGx::setupSurvivalReport(geneIds, geneDirecs, genesInSig)
  #geneIds = genesOfInt$Gene.Symbol
  #geneDirecs = genesOfInt$t
  #genesInSig = genesOfInt$geneInSig
  if(numGroups > 5)
    stop("The numGroups variable cannot be greater than 5")

  if(length(geneSigList) != length(geneDirecList))
    stop("There is not a 1 to 1 relationship between the geneSigList supplied and the geneDirecList supplied (lengths of each list are not equal)")

  if(is.null(patientScoresList) == FALSE){
    for(i in 1:length(patientScoresList))
      if((sum(length(patientNames) != length(patientScoresList[[i]]))))
        stop("patientScoresList elements must contain a single score for each patient in patientNames")
    }

  if(!is.null(patientScoresList) & is.null(patientNames))
    stop("Cannot provide patientScores without patientNames")

  for(i in 1:length(geneSigList))
    if(length(geneSigList[[i]]) != length(geneDirecList[[i]]))
      stop(paste0("There is not a 1 to 1 relationship between the gene IDs in geneSigList and the directions in geneDirecList (length of geneSigList[[", i,"]] is not equal to the length of geneDirecList[[",i,"]])"))

  #file is placed in the directory location the knitrenv is in, so currently place file in directory at time of function call
  #setwd(origDirec)
  knitrenv <- new.env()
  assign("geneSigList", geneSigList, knitrenv)
  assign("geneDirecList", geneDirecList, knitrenv)
  assign("cancerType", cancerType, knitrenv)
  assign("subtype", subtype, knitrenv)
  assign("survivalMetric", survivalMetric, knitrenv)
  assign("numGroups", numGroups, knitrenv)
  assign("patientNames", patientNames, knitrenv)
  assign("patientScoresList", patientScoresList, knitrenv)
  assign("dataNames", dataNames, knitrenv)
  assign("soloGeneAnalysis", soloGeneAnalysis, knitrenv)
  assign("removeMid", removeMid, knitrenv)
  assign("censorTime", censorTime, knitrenv)
  assign("addBenchmarks", addBenchmarks, knitrenv)
  assign("genesRequired", genesRequired, knitrenv)
  #setwd(system.file("latex", package = "metaGx"))
  knitInput = system.file("latex", "MetaGxReport.Rnw", package = "metaGx")
  #knit2pdf(input = knitInput, output = knitOutput)
  knit2pdf(input = knitInput, envir=knitrenv)
  #2:38 start
}

#TO DO
#4. make sure report works for each subtype (check that plots are good for case with > 4 subtypes)
#5. add options to report functions where needed, eg censortime, normalizeEsets, and combineD model parameter
#5. create vignette
#6. add functions that haibe-kains wants
#7. write paper+publish package

#library(knitr)
#setwd("C://Users//micha//Documents//PMH Research//metaGx//inst//latex")
#knit2pdf("MetaGxReport.Rnw")

#origDirec = getwd()
#library(gdata)
#setwd("C:\\Users\\micha\\Documents\\PMH Research\\Yale Alberto Meta Analysis\\ProjectTest\\Sweave Project MetaGx Fresh\\R Code")
#fileName = "Paul Gene Signature Small"
#fileId = paste0(fileName,".xls")
#genesOfInt = read.xls(fileId, perl = "C:\\Perl64\\bin\\perl.exe")
#setwd(origDirec)
#geneIds = genesOfInt$Gene.Symbol
#geneDirecs = genesOfInt$t
#genesInSig = genesOfInt$geneInSig
#cancerType = "ovarian"
#numGroups = 2
#subtype = "consensusOv"
#survivalMetric = "overall"
#soloGeneAnalysis = FALSE
#censorTime = 10
#removeMid = 0
#addBenchmarks = FALSE
#dataNames = NULL
#genesRequired = 0.750
#includeAll = FALSE
#patientNames = NULL
#patientScores = NULL
#geneSigList = geneSigListOvarian
#geneDirecList = geneDirecListOvarian
#ovarianSigs = obtainCancerSigs("ovarian")
#geneSigList = ovarianSigs$geneSigList
#geneDirecList = ovarianSigs$geneDirecList

#library("devtools")
#library(roxygen2)
#setwd("C:\\Users\\micha\\Documents\\PMH Research\\metaGx")
#document()
#setwd("..")
#install("metaGx")
#library(metaGx)
#createSurvivalReport(geneIds, geneDirecs, "ovarian", "verhaak", "overall", numGroups = 2, genesInSig = genesInSig, removeMid = 1)

#Paul Analysis
#genesSigList = list(c("43847", "434768", "80070", "3620", "51561", "50616"))
#names(genesSigList) = c("Gene Signature")
#geneDirecList = list(c(1, 1, 1, 1, 1, 1))
#createSurvivalReport(genesSigList, geneDirecList, cancerType = "ovarian", subtype = "verhaak", survivalMetric = "overall", soloGeneAnalysis = TRUE, removeMid = 0.33333, censorTime = 5)

#library("devtools")
#library(roxygen2)
#setwd("C:\\Users\\micha\\Documents\\PMH Research\\MetaGxBreast")
#document()
#setwd("..")
#install("MetaGxBreast")
#library(MetaGxBreast)

#library("devtools")
#library(roxygen2)
#setwd("C:\\Users\\micha\\Documents\\PMH Research\\MetaGxOvarian_0.99.0\\MetaGxOvarian")
#document()
#setwd("..")
#install("MetaGxOvarian")
#library(MetaGxOvarian)

#ovarian sigs

#Paper: 10-gene biomarker panel: a new hope for ovarian cancer?
#below decent with relapse free survival, censored at 5 years
#geneIds10Sig = c("AEBP1", "COL11A1", "COL5A1", "COL6A2", "LOX", "POSTN", "SNAI2", "THBS2", "TIMP3", "VCAN")
#geneDirecs10Sig = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

#Paper: A Gene Signature Predicting for Survival in Suboptimally Debulked Patients with Ovarian Cancer
#below decent with relapse free survival, censored at 5 years
#setwd("C:\\Users\\micha\\Documents\\PMH Research\\metaGx Paper\\Ovarian Signatures")
#geneSig = read.xlsx(file = "Bioinformatics Sig Paper Signature.xlsx", sheetIndex = 1)
#geneSig = geneSig[!duplicated(geneSig$Entrez.ID), ]
#geneIds = geneSig$Entrez.ID
#geneDirec = sign(as.numeric(as.character(gsub("?^'", "-", geneSig$Cox.score.coefficient))))

#Paper:  Integrated genomic analyses of ovarian carcinoma
#very prognostic
#setwd("C:\\Users\\micha\\Documents\\PMH Research\\metaGx Paper\\Ovarian Signatures")
#geneSig = read.xlsx(file = "Verhaak Paper Signature.xlsx", sheetIndex = 1)
#geneDirec = as.character(geneSig$Gene.set)
#geneIds = as.character(geneSig$Name)
#geneDirec[geneDirec == "poor"] = 1
#geneDirec[geneDirec == "good"] = -1
#names(geneDirec) = NULL
#geneDirec = as.numeric(as.character(geneDirec))

#Paper: A DNA Repair Pathway-Focused Score for Prediction of Outcomes in Ovarian Cancer Treated With Platinum-Based Chemotherapy
#not prognostic, weird method section , different score calculation then calcSigScore
#setwd("C:\\Users\\micha\\Documents\\PMH Research\\metaGx Paper\\Ovarian Signatures")
#geneSig = read.xlsx(file = "Kang Paper Signature.xlsx", sheetIndex = 1)
#geneIds = c("ATM", "H2AFX",    "MDC1",     "RNF8",     "TOP2A",    "BRCA2" ,    "C17orf70", "FANCB",    "FANCE" ,   "FANCF" ,   "FANCG",
#             "FANCI" ,   "PALB2",    "MUS81" ,   "NBN"   ,   "SHFM1" ,   "DDB1" ,    "ERCC8" ,   "RAD23A" ,  "XPA"  ,    "MAD2L2" ,  "POLH" , "UBE2I")
#geneDirec = c(1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, -1 ,1 , 1, -1 , 1, -1, 1, 1, 1)

#Paper: Tumor-Infiltrating Plasma Cells Are Associated with Tertiary Lymphoid Structures, Cytolytic T-Cell Responses, and Superior Prognosis in Ovarian Cancer
#Not prognostic? but could be that signatur eisnt used in same context
#setwd("C:\\Users\\micha\\Documents\\PMH Research\\metaGx Paper\\Ovarian Signatures")
#geneSig = read.xlsx(file = "Plasma Cells Paper Signature.xlsx", sheetIndex = 1)
#geneIds = as.character(geneSig$PCvsno.PC.B.cells.gene)
#geneDirec = rep(1, length(geneIds))

#Does not appear to be prognostic at all?
#Paper: The prognostic signi???cance of speci???c HOX gene expressionpatterns in ovarian cancer
#geneIds = c("HOXC13", "HOXB6", "HOXA13", "HOXD13", "HOXD1")
#geneDirec = c(1, 1, 1, 1, 1)

#setwd("C:\\Users\\micha\\Documents\\PMH Research")
#geneSigList = list(geneIds)
#names(geneSigList) = c("Gene Signature")
#geneDirecList = list(geneDirec)
#createSurvivalReport(geneSigList, geneDirecList, cancerType = "ovarian", subtype = "verhaak", survivalMetric = "relapse", censorTime = 5)



