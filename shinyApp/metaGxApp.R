#konecny crashed but verhaak didnt crash at same point

############ OUTLINE ##############
#first page: select data, subtype, survival type etc and display survival plots/generate pdf
#second page: show analysis information
#third page: show survival plots and forest plots
#fourth page: show table of most prognostic genes, allow user to select individual gene and see survival plots for it

########### LIBRARIES ############
library(shiny)
library(metaGx)
library(xlsx)
library(survival)
library(d3heatmap)
library(gdata)
library(forestplot)
library(xlsx)
library(GSVA)
library(mclust)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(HiDimDA)
library(matrixStats)
library(survival)
library(survcomp)
library(MetaGxOvarian)
library(MetaGxBreast)
library(gplots)
library(knitr)
library(XML)
library(httr)
library(Biobase)
#below 3 needed for breast subtype code
library(AIMS)
library(iC10)
library(limma)


## (C) Wail
# Function used to link to ensembl database
#createLink <- function(val) {
#  sprintf('<b><a style="color: blue; background: Transparent; border: none;outline:none;" href="http://useast.ensembl.org/Homo_sapiens/Gene/Idhistory?g=%s" target="_blank" class="btn btn-primary">Link</a></b>',val)
#}

########### SOURCE FUNTIONS & LOAD DATA #############
#need to generate metaGxReport in shiny directory in order for user to download it from app directory later
setwd("C:\\Users\\micha\\Documents\\PMH Research\\metaGx Shiny")
breastSubtypes = c("scmod2", "scmod1", "scmgene", "pam50", "ssp2006", "ssp2003", "intClust", "AIMS", "claudinLow")
breastSurvTypes = c("Overall Survival", "Relapse Free Survival", "Overall Metabric Only", "Overall TCGA Only")
ovarianSubtypes = c("verhaak", "bentink", "tothill", "helland", "konecny")
ovarianSurvTypes = c("Overall Survival", "Relapse Free Survival")

userInput = NULL
reportDone = NULL
# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #need to adjust so report can only be downloaded after its generated

  output$dataTable <- renderUI({
    
    dataList = loadMetaData(cancerType(), survivalMetric = survivalMetric())
    dataTab = obtainDataInfo(dataList, survivalMetric = survivalMetric())

    #dateInfo = dataTab$Date
    #print(output$dataTable)
    #print(input$show_vars)
    
    checkboxGroupInput(inputId='show_vars', 'Datasets', paste(dataTab$`Data Name`, dataTab$Date),
                       selected = paste(dataTab$`Data Name`, dataTab$Date))
  }) 
  
  output$dlButton <- renderUI({
    if(length(sigResults()) > 0)
      downloadButton('survivalReport', label = "Download Report")
    
  })
  
  output$genRepButton <- renderUI({
    if(length(input$datafile) > 0)
      actionButton("genReport", "Generate Report")
    
  })
  
  selectedData <- reactive({
    input$show_vars
  })
    
  output$subtype <- renderUI({
    if(input$cancerType == "Breast"){
      subtypeVec = breastSubtypes
    }
    else if(input$cancerType == "Ovarian"){
      subtypeVec = ovarianSubtypes
    }
    #selectInput('subtype', 'Subtype', subtypeVec, multiple=FALSE, selectize=TRUE)
    radioButtons(inputId="subtype", label="Subtype", choices = subtypeVec)
  }) 
  
  output$survivalMetric <- renderUI({
    if(input$cancerType == "Breast"){
      survTypeVec = breastSurvTypes
    }
    else if(input$cancerType == "Ovarian"){
      survTypeVec = ovarianSurvTypes
    }
    radioButtons(inputId="survivalMetric", label="Survival Type", choices = survTypeVec)
  }) 
  
  
  output$survivalReport = downloadHandler(
    filename = 'myreport.pdf',
    
    content = function(file) {
      file.copy('MetaGxReport.pdf', file, overwrite = TRUE)
    },
    
    contentType = 'application/pdf'
  )
  #values = reactiveValues()
  
  #dataset<-reactive({ 
  #  inFile <- input$uploadFile 
  #  dat<-read.xlsx(inFile$datapath, 1)
  #  return(dat)
  #})
  
  #output$summary <- renderText({summary(dataset())})
  
  
  #read csv file
  output$contents <- renderTable({
    # input$datafile will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile <- input$datafile
    
    if(is.null(inFile))
      return(NULL)
    
    read.xlsx(inFile$datapath, 1)
    
  })
  
  
  #infile <- reactive({
  #  infile <- input$datafile
  #  if (is.null(infile)) {
  #    # User has not uploaded a file yet
  #    return(NULL)
  #  }
  #  objectsLoaded <- load(input$datafile$name) 
  #  # the above returns a char vector with names of objects loaded
  #  df <- eval(parse(text=objectsLoaded[1])) 
  #  # the above finds the first object and returns it
  #  return(df)
  #})
  
  #myData <- reactive({
  #  df<-infile()
  #  if (is.null(df)) return(NULL)
  #  return(df)
  #})
  
  #output$value2 <- renderPrint({
  #  names(myData())
  #})
  
  #createReport = reactive({
  #  cancerType = tolower(input$cancerType)
  #  subtype = input$subtype
  #
  #})
  
  #observeEvent(input$createReport, {
  #  session$sendCustomMessage(type = 'testmessage', message = 'Thank you for clicking')
  #})
  
  subtype = reactive({
    shinySub = tolower(input$subtype)
    shinySub
  })
  
  cancerType = reactive({
    shinyCanc = tolower(input$cancerType)
    shinyCanc
  })
  
  numGroups = reactive({
    shinyGroups = as.numeric(input$numGroups)
    shinyGroups
  })
  
  censorTime = reactive({
    shinyCens = as.numeric(input$censorTime)
    shinyCens
  })
  
  removeMid = reactive({
    shinyMidRem = as.numeric(input$removeMid)
    shinyMidRem
  })
  
  survivalMetric = reactive({
    shinySurvMet = as.character(input$survivalMetric)
    shinySurvMet = sub(" ", "", shinySurvMet)
    shinySurvMet
  })
  
  geneSig = reactive({
    inFile <- input$datafile 
    data <- read.xlsx(inFile$datapath, 1)
    shinySig = as.vector(as.character(data[,1]))
    shinySig
  })
  
  geneDirec = reactive({
    inFile <- input$datafile 
    data <- read.xlsx(inFile$datapath, 1)
    shinyDirecs = as.vector(data[,2])
    shinyDirecs
  })
  
  output$text1 <- renderText({
    #inFile <- input$datafile 
    #data <- read.xlsx(inFile$datapath, 1)
    
    paste(paste(subtype()), paste(cancerType()))
    #paste(geneSig())
  })
  
  dataInfoTab = reactive({
    dataTab = obtainDataInfo(sigResults()$datasets, survivalMetric = survivalMetric())
  })
  
  output$dataInfo = DT::renderDataTable(
    
    #addCheckboxButtons <- paste0('<input type="checkbox" name="row', mymtcars$id, '" value="', mymtcars$id, '">',"")
    DT::datatable(dataInfoTab(), escape = FALSE, 
                  options = list(
                    lengthChange = FALSE), rownames = FALSE) 
  )
  
  geneInfoTab = reactive({
    geneInfoFrame = getGeneInfo(geneSig(), geneDirec())
  })
  
  output$geneInfo = DT::renderDataTable(
    DT::datatable(geneInfoTab(), escape = FALSE, 
                  options = list(
                    lengthChange = FALSE), rownames = FALSE) 
  )
  
  output$sessionInfo = renderText(paste(capture.output(sessionInfo()), collapse="<br>"))
  
  
#  createSurvRep <- eventReactive(input$genReport, {
#    
#    print(selectedData())
#    dataNames = c()
#    for(i in 1:length(selectedData()))
#      dataNames = c(dataNames, substr(selectedData()[i], 1, regexpr(" ", selectedData()[i])-1))
#    })
  
  sigResults <- eventReactive(input$genReport, {
    #print(selectedData())
    dataNames = c()
    for(i in 1:length(selectedData()))
      dataNames = c(dataNames, substr(selectedData()[i], 1, regexpr(" ", selectedData()[i])-1))
    #dataNames = NULL
    addBenchmarks = FALSE
    soloGeneAnalysis = TRUE
    
    print(cancerType())
    print(subtype())
    print(survivalMetric())
    if(survivalMetric() == "Overall Survival"){
      
    }
    print(survivalMetric())
    print(subtype())
    print(cancerType())
    print(dataNames)
    geneInfoLists = getGenesProgValue(list(geneSig()), list(geneDirec()), cancerType(), subtype(), survivalMetric(), numGroups(), dataNames, soloGeneAnalysis, removeMid(), censorTime(), addBenchmarks)
    
    
    #load("geneInfoLists.RData")
    knitrenv <- new.env()
    assign("geneSigList", list(geneSig()), knitrenv)
    assign("geneDirecList", list(geneDirec()), knitrenv)
    assign("cancerType", cancerType(), knitrenv)
    assign("subtype", subtype(), knitrenv)
    assign("survivalMetric", survivalMetric(), knitrenv)
    assign("numGroups", numGroups(), knitrenv)
    assign("dataNames", dataNames, knitrenv)
    assign("soloGeneAnalysis", soloGeneAnalysis, knitrenv)
    assign("removeMid", removeMid(), knitrenv)
    assign("censorTime", censorTime(), knitrenv)
    assign("addBenchmarks", addBenchmarks, knitrenv)
    assign("geneInfoLists", geneInfoLists, knitrenv)
    #setwd(system.file("latex", package = "metaGx"))
    knitInput = system.file("latex", "MetaGxAppReport.Rnw", package = "metaGx")
    #knit2pdf(input = knitInput, output = knitOutput)
    knit2pdf(input = knitInput, envir=knitrenv)

    geneInfoLists
    
  })
  
  
  
  output$subNames <- renderUI({
    subNameVec = names(sigResults()$patientSurvData[[1]])
    selectInput('subNames', 'Patient Subtype', subNameVec, multiple=FALSE, selectize=TRUE)
    #radioButtons(inputId="subNames", label="Patient Subtype", choices = subNameVec)
  })
  
  subNames = reactive({
    shinySubnames = tolower(input$subNames)
    shinySubnames
  })
  
#  dataInfoTab = reactive({
#    geneSigInfo = getGenesProgValue(list(geneSig()), list(geneDirec()), cancerType(), subtype(), survivalMetric = "overall")
#    
#  })
  
  output$survPlot <- renderPlot({
    ## TODO:: Fix the range to reasonable values.
    #print(subNames())
    plotInd = which(tolower(names(sigResults()$patientSurvData[[1]])) == tolower(subNames()))
    #print(plotInd)
    if(tolower(subNames()) == "all patients"){
      titleStr = paste("Survival Curve \nUsing All Patients")
    }else{
      titleStr = paste("Survival Curve \nUsing ", toupper(subNames()), "Patients")
    }
    makeSurvivalPlot(sigResults()$patientSurvData[[1]][[plotInd]], as.numeric(numGroups()), normalizeEsetScores = TRUE, titleStr)
  })
  
  output$subNamesForest <- renderUI({
    subNameVec = names(sigResults()$patientSurvData[[1]])
    selectInput('subNamesForest', 'Patient Subtype', subNameVec, multiple=FALSE, selectize=TRUE)
    #radioButtons(inputId="subNames", label="Patient Subtype", choices = subNameVec)
  })
  
  subNamesForest = reactive({
    shinySubnames = tolower(input$subNamesForest)
    shinySubnames
  })
  
  output$forestPlot <- renderPlot({
    plotInd = which(tolower(names(sigResults()$patientSurvData[[1]])) == tolower(subNamesForest()))
    makeForestPlot(sigResults()$patientSurvData[[1]][[plotInd]])
  })
  
  output$corPlot <- renderD3heatmap({
    #make method an input parameter
    geneInfoFrame = getGeneInfo(geneSig(), geneDirec())
    entrezIds = geneInfoFrame$`Entrez ID`
    naGenes = which(is.na(entrezIds))
    if(length(naGenes) > 0)
      entrezIds = entrezIds[-naGenes]
    corMatrix = correlateGenes(entrezIds, cancerType(), survivalMetric = "overall", method = "pearson")
    corPlot = d3heatmap(corMatrix, scale="column", colors="Blues", dendrogram = "none")
    corPlot
  })
  
  output$geneCorPlot <- renderUI({
    d3heatmap::d3heatmapOutput('corPlot', width = "75%", height = "500px")
  })
  
  output$geneSymbFor <- renderUI({
    geneSymbs = names(geneSigInfo$patientSurvData)
    selectInput('geneSymbFor', 'Gene', geneSymbs, multiple=FALSE, selectize=TRUE)
  })
  
  geneSymbFor = reactive({
    shinySymbFor = tolower(input$geneSymbFor)
    shinySymbFor
  })
  
  output$geneForestPlot <- renderPlot({
    print(geneSymbFor())
    geneInd = which(tolower(names(geneSigInfo$patientSurvData)) == tolower(geneSymbFor()))
    print(geneInd)
    makeForestPlot(sigResults()$patientSurvData[[geneInd]][[1]])
  })
  
  
  output$geneSymbSurv <- renderUI({
    geneSymbs = names(geneSigInfo$patientSurvData)
    selectInput('geneSymbSurv', 'Gene', geneSymbs, multiple=FALSE, selectize=TRUE)
  })
  
  geneSymbSurv = reactive({
    shinySymbSurv = tolower(input$geneSymbSurv)
    shinySymbSurv
  })
  
  output$geneSurvPlot <- renderPlot({
    geneInd = which(tolower(names(geneSigInfo$patientSurvData)) == tolower(geneSymbSurv()))
    ## TODO:: Fix the range to reasonable values.
    #print(subNames())
    #plotInd = which(tolower(names(sigResults()$patientSurvData[[1]])) == tolower(subNames()))
    #print(plotInd)
    #if(tolower(subNames()) == "all patients"){
    #  titleStr = paste("Survival Curve \nUsing All Patients")
    #}else{
    #  titleStr = paste("Survival Curve \nUsing ", toupper(subNames()), "Patients")
    #}
    makeSurvivalPlot(sigResults()$patientSurvData[[geneInd]][[1]], as.numeric(numGroups()), normalizeEsetScores = TRUE, "titleStr")
  })
  
  

}

#note tested changes to drug names and verified results remained the same, aka pkiSenseSigs and ovarianPki are in the same order

# Define UI for application that draws a histogram
UI <- navbarPage(title = "Gene Signature Analysis", 
                 # Application title
                 tabPanel(title = "Signature Data Input",
                          
                          # Sidebar with a slider input for number of bins 
                          #sidebarLayout(
                            sidebarPanel(
                              fileInput("datafile", "XLSX file"),
                              tags$hr(),
                              checkboxInput("header", "Header", TRUE)
                            ),
                            
                            sidebarPanel(
                              selectInput('cancerType', 'Cancer', c("Ovarian", "Breast"), multiple=FALSE, selectize=TRUE)
                            ),
                            
                          
                          uiOutput("subtype"),
                          uiOutput("survivalMetric"),
                          uiOutput("dataTable"),
                          

                          sidebarPanel(
                            sliderInput("numGroups", "Number of Patient/Survival Groups", min=2, max=4, value=2)
                          ),
                          
                          sidebarPanel(
                            sliderInput("removeMid", "Percent of Middle Score Patients to Remove", min=0, max=0.5, value=0)
                          ),
                          
                          sidebarPanel(
                            sliderInput("censorTime", "Year to Censor the Survival Data", min=5, max=15, value=10)
                          ),

                            # Show the contents of the file
                            #Need to adjust later so user scrolls through pages to see genes, so table does not blow up
                            #Add cancer type, survival metric, and subtype inputs
                            mainPanel(textOutput("text1"),
                              tableOutput("contents")
                            ),
                          
                          
                          #actionButton("genReport", "Generate Report"),
                          uiOutput("genRepButton"),
                          uiOutput("dlButton")
                          #need to adjust so report can only be downloaded after its generated, button only appears after generated?
                          #downloadButton('survivalReport', label = "Download Report")
                          
                          #)
                 ),
                 navbarMenu(title = "Analysis Information",
                            tabPanel(title = "Dataset Information",
                                     mainPanel(
                                       DT::dataTableOutput("dataInfo")
                                     )
                            ),
                            tabPanel(title = "Gene Information",
                                     mainPanel(
                                       DT::dataTableOutput("geneInfo")
                                     )
                            ),
                            tabPanel(title = "R Session Information",
                                     mainPanel(
                                       htmlOutput("sessionInfo")
                                     )
                            )
                          
                 ),
                 navbarMenu(title = "Signature Results",
                            tabPanel(title = "Survival Plots",
                                     uiOutput("subNames"),
                                     mainPanel(
                                       plotOutput("survPlot")
                                     )
                            ),
                            tabPanel(title = "Forest Plots",
                                     uiOutput("subNamesForest"),
                                     mainPanel(
                                       plotOutput("forestPlot")
                                     )
                            ),
                            tabPanel(title = "Signature Gene Correlations",
                                     mainPanel(
                                       shiny::uiOutput("geneCorPlot")
                                     )
                            )
                 ),
                 navbarMenu(title = "Individual Gene Results",
                            tabPanel(title = "Survival Plots",
                                     uiOutput("geneSymbSurv"),
                                     mainPanel(
                                       plotOutput("geneSurvPlot")
                                     )
                            ),
                            tabPanel(title = "Forest Plots",
                                     uiOutput("geneSymbFor"),
                                     mainPanel(
                                       plotOutput("geneForestPlot")
                                     )
                            )
                 )
)

shinyApp(ui = UI, server = server)



