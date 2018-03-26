#' Function to identify which probes in an esets expression data that should be kept amongst those with the same IDs
#'
#' This function returns the rows indices of the probes in eset@assayData$exprs that should be kept amongst probes with the same IDs. Probes with the highest IQR amongst those with the same IDs
#' are kept.
#' @param eset an eset object
#' @return a vector containing the row indices of the probes in the eset@assayData$exprs expresson matrix that have the highest IQR for each probe ID. The names of the bestProbes vector are the entrez IDs
#' @export
#' @examples
#'
#' eset = loadMetaData("ovarian", "overall")[[1]]
#' bestProbeRows = getBestProbes(eset)
#' expressMat = eset@assayData$exprs[bestProbeRows, ]
#'

getBestProbes = function(eset)
{
  
  #dataValsOrig = eset@assayData$exprs;
  #Result is much worse if the below line is included
  #dataVals <- normalizeBetweenArrays(dataVals,method="quantile")
  #dataInfo = eset@featureData@data;
  #dataInfoEntrezGene.ID = dataInfo$EntrezGene.ID
  
  #below seems uneccessary now that createEsetList.R has been fixed
  #dataInfoEntrezGene.ID = as.character(dataInfo$EntrezGene.ID)
  ##some missing ids that are "///"
  #slashStarts = which(substr(dataInfoEntrezGene.ID, 1, 1) == "/")
  #slashStLeng = nchar(dataInfoEntrezGene.ID[slashStarts])
  #slashRem = substr(dataInfoEntrezGene.ID[slashStarts], 4, slashStLeng)
  #dataInfoEntrezGene.ID[slashStarts] = slashRem
  #missingEnts = which(dataInfoEntrezGene.ID == "")
  ##set the missing ents to -100 so the probe selection works
  ##dont end up using any of them anyways
  #dataInfoEntrezGene.ID[missingEnts]  = -100
  
  ##multiple entrez ids have ent1///ent2///....///entN
  #multEnts = which(grepl("/", dataInfoEntrezGene.ID))
  #multStrings = dataInfoEntrezGene.ID[multEnts]
  #firstSlashes = rapply(gregexpr(pattern = "/", multStrings), function(x) head(x, 1))
  #firstEnts = substr(multStrings, 1, firstSlashes - 1)
  #dataInfoEntrezGene.ID[multEnts] = firstEnts
  
  #multEnts = which(grepl(",", dataInfoEntrezGene.ID))
  #multStrings = dataInfoEntrezGene.ID[multEnts]
  #firstSlashes = rapply(gregexpr(pattern = ",", multStrings), function(x) head(x, 1))
  #firstEnts = substr(multStrings, 1, firstSlashes - 1)
  #dataInfoEntrezGene.ID[multEnts] = firstEnts
  
  #dataInfoEntrezGene.ID = as.numeric(dataInfoEntrezGene.ID)
  
  #dataInfoProbeset = as.vector(dataInfo$probeset)
  #dataInfoGene = as.vector(dataInfo$gene)
  
  #repeatEnt = duplicated(dataInfoEntrezGene.ID)
  #repeatEnt = which(repeatEnt == TRUE)
  
  #rmList = c()
  #for(i in 1:length(repeatEnt))
  #{
  #  repInd = repeatEnt[i]
  #  entId = dataInfoEntrezGene.ID[repInd]
  #  probesWithId = which(dataInfoEntrezGene.ID == entId)
  #  probes = dataValsOrig[probesWithId,]
  #  iqrs = rowIQRs(probes);
  #  winnerInd = which(iqrs == max(iqrs))
  #  keeper = probesWithId[winnerInd]
  #  probesRm = probesWithId[-winnerInd]
  #  rmList = c(rmList, probesRm)
  #  bestProbes = c(bestProbes, keeper)
  #}
  #above code checks more then once per repeat id, remove non unique
  #rmList = unique(rmList)
  #bestProbes is now a just a list of all the rows to keep in the data set
  #bestProbes = c(1:length(dataInfoEntrezGene.ID))
  #bestProbes = bestProbes[-rmList]
  
  #return(bestProbes)
  
  #0. name ent id vec according to position in vec
  #1. find ents with ///
  #2. get individual ent ids
  #3. rename that row in expression frame to first ent
  #4. append new rows to expression data with same values but rowname the other ent ids
  #   make sure name in ent id vec is original index
  #5. sort data and ent vec by numeric id
  #6. run through while loop but select best probe by name in ent id vec
  
  dataVals = eset@assayData$exprs;
  dataInfo = eset@featureData@data;
  dataInfoEntrezGene.ID = dataInfo$EntrezGene.ID
  
  dataInfoProbeset = as.vector(dataInfo$probeset)
  dataInfoGene = as.vector(dataInfo$gene)
  
  #old method below that is faster, should probably switch back to this and add a reordering
  #deal with / and , in entrez id list
  dataInfoEntrezGene.ID = as.character(dataInfo$EntrezGene.ID)
  dataInfoEntrezGene.ID = gsub(" ", "", dataInfoEntrezGene.ID)
  dataInfoEntrezGene.ID = gsub(",", "///", dataInfoEntrezGene.ID)
  dataInfoEntrezGene.ID = gsub("//////", "///", dataInfoEntrezGene.ID)
  names(dataInfoEntrezGene.ID) = seq(1, length(dataInfoEntrezGene.ID))
  
  missingEnts = which(dataInfoEntrezGene.ID == "///")
  #one weird dataset where ents are missing, set to not a real ent so 
  #wont be selected later by user
  if(length(missingEnts) > 0) dataInfoEntrezGene.ID[missingEnts] = "0"
  
  #code that follows relys on their not being a /// at start of entrez
  #so remove cases where there is. Now /// separates ents only
  startSlashes = which(grepl("///", substr(dataInfoEntrezGene.ID, 1, 3)))
  if(length(startSlashes) > 0)
    dataInfoEntrezGene.ID[startSlashes] = substr(dataInfoEntrezGene.ID[startSlashes], nchar("///")+1, nchar(dataInfoEntrezGene.ID[startSlashes]))
  endSlashes = which(grepl("///", substr(dataInfoEntrezGene.ID, nchar(dataInfoEntrezGene.ID)-2, nchar(dataInfoEntrezGene.ID))))
  if(length(endSlashes) > 0)
    dataInfoEntrezGene.ID[endSlashes] = substr(dataInfoEntrezGene.ID[endSlashes], 1, nchar(dataInfoEntrezGene.ID[endSlashes]) - 3)
  
  multNameInds = which(grepl("///", dataInfoEntrezGene.ID))
  
  dataInfoEntrezGene.IDNew = dataInfoEntrezGene.ID
  #time consuming process due to need to compare ents with ///
  #and to make code run faster by ordering ents by number
  if(length(multNameInds) > 0)
  {
    for(i in 1:length(multNameInds))
    {
      slashInds = gregexpr("///", dataInfoEntrezGene.ID[multNameInds[i]])[[1]]
      slashInds = c(slashInds, nchar(dataInfoEntrezGene.ID[multNameInds[i]]) + 1)
      startInd = 0
      for(j in 1:length(slashInds))
      {
        entName = substr(dataInfoEntrezGene.ID[multNameInds[i]], startInd, slashInds[j]-1)
        if(is.na(entName)) stop()
        if(j == 1){
          dataInfoEntrezGene.IDNew[multNameInds[i]] = entName
          #rownames(dataVals)[multNameInds[i]] = entName
        }else{
          dataInfoEntrezGene.IDNew[length(dataInfoEntrezGene.IDNew)+1] = entName
          names(dataInfoEntrezGene.IDNew)[length(dataInfoEntrezGene.IDNew)] = names(entName)
          dataVals = rbind(dataVals, dataVals[multNameInds[i], ])
          rownames(dataVals)[nrow(dataVals)] = entName
          #print(geneName)
        }
        startInd = slashInds[j] + 3
      }
    }
  }
  
  #multEnts = which(grepl("/", dataInfoEntrezGene.ID))
  #multStrings = dataInfoEntrezGene.ID[multEnts]
  #firstSlashes = rapply(gregexpr(pattern = "/", multStrings), function(x) head(x, 1))
  #firstEnts = substr(multStrings, 1, firstSlashes - 1)
  #dataInfoEntrezGene.ID[multEnts] = firstEnts
  
  
  #dataInfoEntrezGene.ID = paste0("/", dataInfoEntrezGene.ID, "/")
  
  
  #naRows = which(rowSums((is.na(dataVals))) > 0)
  #if(length(naRows) > 0)
  #{
  #  dataInfoEntrezGene.ID = dataInfoEntrezGene.ID[-naRows]
  #  dataVals = dataVals[-naRows, ];
  #}
  dataValsOrig = dataVals
  namesEntrezNew = names(dataInfoEntrezGene.IDNew)
  dataInfoEntrezGene.IDNew = as.numeric(dataInfoEntrezGene.IDNew)
  names(dataInfoEntrezGene.IDNew) = namesEntrezNew
  sortInd = sort(dataInfoEntrezGene.IDNew, decreasing = FALSE, index.return=TRUE)$ix;
  dataInfoEntrezGene.IDPreSort = dataInfoEntrezGene.ID
  dataInfoEntrezGene.ID = dataInfoEntrezGene.ID[sortInd];
  dataInfoEntrezGene.IDNew = dataInfoEntrezGene.IDNew[sortInd]
  dataVals = dataValsOrig[sortInd, ];
  dataInfoProbeset = dataInfoProbeset[sortInd];
  dataInfoGene = dataInfoGene[sortInd];
  
  
  bestProbes = c()
  i = 1;
  
  while(i < nrow(dataVals))
  {
    entrezId = dataInfoEntrezGene.IDNew[i];
    #slashLoc = gregexpr(pattern = "/", entrezId)[[1]][1]
    #if(slashLoc > 0)
    #  entrezId = substr(entrezId, 1, slashLoc-1)
    origInd = i;
    probesWithId = c();
    while(dataInfoEntrezGene.IDNew[i] == entrezId)
    {
      probesWithId = c(probesWithId, i);
      i = i + 1;
      if(i == dim(dataVals)[1])
      {
        if(dataInfoEntrezGene.IDNew[i] == entrezId)
          probesWithId = c(probesWithId, i)
        if(dataInfoEntrezGene.IDNew[i] != entrezId){
          bestProbes = c(bestProbes, names(dataInfoEntrezGene.IDNew)[i])
          #names(bestProbes)[length(bestProbes)] = as.character(dataInfoEntrezGene.ID)[as.numeric(names(dataInfoEntrezGene.IDNew)[i])]
        }
        
        break
      }
    }
    if(length(probesWithId) > 1)
    {
      #print(origInd)
      probes = dataVals[probesWithId,]
      iqrs = sapply(1:nrow(probes), function(x) matrixStats::iqr(probes[x, ]))
      #iqrs = rowIQRs(probes[!rowSums(is.na(probes)) == ncol(probes), !colSums(is.na(probes)) == nrow(probes)]);
      keepProbe = which(iqrs == max(iqrs)) + (origInd - 1);
      bestProbes = c(bestProbes, names(dataInfoEntrezGene.IDNew)[keepProbe])
      #names(bestProbes)[length(bestProbes)] = as.character(dataInfoEntrezGene.ID)[as.numeric(names(dataInfoEntrezGene.IDNew)[keepProbe])]
    }
    if(length(probesWithId) == 1)
    {
      bestProbes = c(bestProbes, names(dataInfoEntrezGene.IDNew)[origInd])
      #names(bestProbes)[length(bestProbes)] = as.character(dataInfoEntrezGene.ID)[as.numeric(names(dataInfoEntrezGene.IDNew)[origInd])]
      
    }
  }
  names(bestProbes) = dataInfo$EntrezGene.ID[as.numeric(bestProbes)]
  names(bestProbes) = gsub(" ", "", names(bestProbes))
  names(bestProbes) = gsub(",", "///", names(bestProbes))
  #potential duplicates from splitting up /// ent rows into multiple rows with same values
  #and then using nmes(dataInfoEnts) to only keep the row index of the original row
  bestProbes = bestProbes[!duplicated(bestProbes)]
  
  #get the rownames of the probes to keep and then match them with the rows from the unsorted data
  #keepRows = rownames(dataVals)[bestProbes]
  #bestProbes = match(keepRows, rownames(dataValsOrig))
  return(bestProbes)
}

