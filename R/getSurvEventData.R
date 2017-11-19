#' Function to obtain events and event times for a given survival metric in each dataset
#'
#' This function returns the event times and event statuses for the patients in the datasets
#' @param dataList a list containing the datasets to be used in the analysis
#' @param survivalMetric a string specifying what type of survival data to use in the analysis. If the cancerType is set to "ovarian"
#' then enter "overall" for overall survival and "relapse" for relapse free survival. If the cancerType is set to "breast" then the options
#' are "relapse" for relapse free survival, "overallMeta" for overall survival on the metabric study data, "overallTcga" for
#' overall survival on the TCGA study data, and "overall" for overall survival on the remaining studies (No Metabric and TCGA data).
#' @return a list with two lists. One called eventList that indicates whether the event occurred or not for the patients of each dataset and the
#' second list, called eventTimeList, indicates the time (in years) that the events in eventList were observed relative to the patients first checkup 
#' @export
#' @examples
#'
#' dataList = loadMetaData(cancerType = "breast", survivalMetric = "overall")
#' survEventList = getSurvEventData(dataList, survivalMetric = "overall")

getSurvEventData = function(dataList, survivalMetric)
{
 
  getEsetEvents = function(eset, survivalMetric)
  {
    if(grepl("overall", survivalMetric, ignore.case = TRUE)){
      evTime = eset@phenoData@data$days_to_death/365.25;
      events = as.character(eset@phenoData@data$vital_status)
      events[events == "living"] = 0
      events[events == "deceased"] = 1
    }else if(grepl("relapse", survivalMetric, ignore.case = TRUE)){
      evTime = eset@phenoData@data$days_to_tumor_recurrence/365.25;
      events = as.character(eset@phenoData@data$recurrence_status)
      events[events == "norecurrence"] = 0
      events[events == "recurrence"] = 1
    }else if(grepl("metastasis", survivalMetric, ignore.case = TRUE)){
      evTime = eset@phenoData@data$dmfs_days/365.25;
      events = as.character(eset@phenoData@data$dmfs_status)
      events[events == "norecurrence"] = 0
      events[events == "recurrence"] = 1
    }else{
      stop("The survivalMetric parameter must be either overall, relapse, metastisis, or hierarchy")
    }
    events = as.numeric(events)
    return(list(events, evTime))
  }
  eventList = list()
  eventTimeList = list()
  for(i in 1:length(dataList))
  {
    eset = dataList[[i]]
    if(survivalMetric == "hierarchy"){
      eventData = getEsetEvents(eset, "relapse")
      
      if((sum(is.na(eventData[[1]])) == length(eventData[[1]])) | (sum(is.na(eventData[[2]])) == length(eventData[[2]])))
        eventData = getEsetEvents(eset, "metastasis")
      
      if((sum(is.na(eventData[[1]])) == length(eventData[[1]])) | (sum(is.na(eventData[[2]])) == length(eventData[[2]])))
        eventData = getEsetEvents(eset, "overall")
      
    }else{
      eventData = getEsetEvents(eset, survivalMetric)
    }
    eventList[[i]] = eventData[[1]]
    eventTimeList[[i]] = eventData[[2]]
  }
  names(eventList) = names(dataList)
  names(eventTimeList) = names(dataList)
  
  survEventList = list(eventList, eventTimeList)
  names(survEventList) = c("eventList", "eventTimeList")
  return(survEventList)
   
}