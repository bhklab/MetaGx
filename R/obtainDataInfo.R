#' Function to obtain information about various esets
#'
#' This function returns a data frame which describes the number of events (deceased or tumor recurrence), median time to an event, the number of patients, the number of genes, the platform of an eset, and the publication date (NA if no internet connection)
#' @param dataList a list containing the datasets one would like info about
#' @param survivalMetric a string specifying whether to use overall survival or relapse free survival as the event of interest. Note that
#' an error will be thrown if the esets in dataList do not contain survival information for the requested event. Enter "overall" for
#' overall survival and "relapse" for relapse free survival.
#' @return a data frame with information about the esets supplied
#' @export
#' @examples
#'
#' dataList = loadMetaData(cancerType = "ovarian", survivalMetric = "overall")
#' dataInfoTab = obtainDataInfo(dataList, survivalMetric = "overall")
#' head(dataInfoTab)

obtainDataInfo = function(dataList, survivalMetric)
{
  #dont censor data for table (Haibe-Kains)
  dataInfoFrame = as.data.frame(NULL)
  survEventList = getSurvEventData(dataList, survivalMetric)
  
  if(grepl("overall", tolower(survivalMetric)))
  {
    numEvStr = "# Deceased"
    evTimeStr = "Median Survival Time (Years)"
  }
  if(grepl("relapse", tolower(survivalMetric)))
  {
    numEvStr = "# Tumor Recurrences"
    evTimeStr = "Median Recurrence Time (Years)"
  }
  if(grepl("metastasis", tolower(survivalMetric)))
  {
    numEvStr = "# Distant Metastases"
    evTimeStr = "Median Appearence Time (Years)"
  }
  if(grepl("hierarchy", tolower(survivalMetric)))
  {
    numEvStr = "# Events (Metastases/Recurrence/Deceased)"
    evTimeStr = "Median Event Time (Years)"
  }
  
  
  rowInfo = c()
  for(i in 1:length(dataList))
  {

    evTime = survEventList$eventTimeList[[i]]
    events = survEventList$eventList[[i]]
    if((sum(is.na(evTime)) == length(evTime)) | (sum(is.na(events)) == length(events)))
      stop(paste("esets in the data list do not have", survivalMetric, "survival data"))
    survInf = survfit(formula = Surv(evTime, events) ~ 1)
    medEvTime = as.numeric(summary(survInf)$table["median"])
    numPatients = as.integer(summary(survInf)$table["records"])
    numEvents = as.integer(summary(survInf)$table["events"])
    #rowInfo = c(data name, number patients, number of genes, number of events, median event time, platform)
    platformStr = as.character(dataList[[i]]@experimentData@other$platform_shorttitle)
    if(length(platformStr) == 0)
      platformStr = "Unknown"
    
    has_internet = function(){
      !is.null(curl::nslookup("r-project.org", error = FALSE))
    }

    if(Biobase::testBioCConnection()){
      #print(i)
      dataUrl = dataList[[i]]@experimentData@contact
      if(dataUrl == "" & dataList[[i]]@experimentData@pubMedIds != "")
        dataUrl = paste0("https://www.ncbi.nlm.nih.gov/pubmed/?term=", dataList[[i]]@experimentData@pubMedIds)
      colonInd = regexpr(";", dataUrl)[1]
      if(colonInd != -1)
        dataUrl = substr(dataUrl, 1, colonInd - 1)
      #print(dataUrl)
      if(dataUrl != "")
        doc <- htmlParse(rawToChar(GET(gsub(" ", "", dataUrl))$content))
      
      if(dataUrl != "")
      {
        if(length(xpathApply(doc, '//h1')) > 1)
        {
          doc.text = xpathApply(doc, '//meta')
          dateTextInd = which(lapply(1:length(doc.text), function(x) xmlGetAttr(doc.text[[x]], "name")) == "description")
          if(length(dateTextInd) > 0)
          {
            dateText = xmlGetAttr(doc.text[[dateTextInd]], "content")
            colonInd = regexpr(";", dateText)[1]
            dateStr = substr(dateText, 1, colonInd - 1) 
          }else{
            dateTextInd = which(lapply(1:length(doc.text), function(x) xmlGetAttr(doc.text[[x]], "name")) == "citation_date")
            if(length(dateTextInd) > 0){
              dateStr = xmlGetAttr(doc.text[[dateTextInd]], "content")
            }else{
              dateStr = "Unknown"
            }
          }
          
        }else{
          dateStr = "Unknown"
        }
      }else{
        dateStr = "Unknown"
      }
      #print(dateStr)
    }else{
      dateStr = "NA"
    }

    rowInfo = c(names(dataList)[i], numPatients, length(unique(dataList[[i]]@featureData@data$EntrezGene.ID)), numEvents, sprintf("%.3g",medEvTime), platformStr, dateStr)
    dataInfoFrame = rbind(dataInfoFrame, rowInfo, stringsAsFactors = FALSE)
  }
  colnames(dataInfoFrame) = c("Data Name", "# Patients", "# Genes", numEvStr, evTimeStr, "Platform", "Date")
  #colNameVec = c("Data Name", "# Patients", "# Genes", numEvStr, evTimeStr, "platform")

  return(dataInfoFrame)
}

