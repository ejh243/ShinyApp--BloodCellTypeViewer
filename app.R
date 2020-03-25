# EH's editted version of GM's orginal code

library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(dplyr)
library(dbplyr)
library(RSQLite)

# Define UI for application
ui<- dashboardPage(
  dashboardHeader(title = "Quantifying the contribution of specific blood cell types to whole blood DNA methylation profiles", titleWidth = 1300),
  dashboardSidebar(disable=TRUE),
  dashboardBody(
    fluidRow(
      column(width = 4,
             box(
               title = "Description", width = NULL,
               "This tool profiles how variation in five different blood cell types influences variation in DNA methylation measured in whole blood.It also provides co-variation statistics for three peripheral tissues: whole blood, buccal and saliva. This information for a single CpG can be visulaised in the plots on the left or for a list of multiple CpGs a batch query can be performed and summary statistics downloaded." , 
                br(),br(),
                "The DNA methylation data behind this tool were generated with Illumina EPIC array data from 30 individuals. Each participant provided a whole blood, buccal epithelial and nasal epithelial sample. Whole blood samples were processed using fluorescence-activated cell sorting (FACS) to purify and isolate five cell-types (monocytes, granulocytes, CD4 T-cells, CD8 T-cells, and B-cells). These data are described in our manuscript which is currently under review.",
                br(),br(),
                textInput("probe", "EPIC array CpG ID:", value="cg10142008"),
				textOutput("downloadtext"),
              HTML("<br/>"),
              downloadButton("downloaddata1","Download summary statistics")
  
             ),
      
             box(
               title = "Batch Query", width = NULL,
               
				fileInput("file1", "Choose CSV file with list of CpG identifiers",
                multiple = FALSE,
                accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
              textOutput("downloadtext"),
              HTML("<br/>"),
              downloadButton("downloaddata2","Download summary statistics")
            )
          ),

      column(width = 8,
             box(
               title=uiOutput("boxplottitle"), width = NULL, height=370,
               withSpinner(plotOutput("boxplot"))
             ),
             
             box(
               title=uiOutput("scattertitle1"), width=NULL, height=500,
               withSpinner(plotOutput("scatterplots1"))
             ),
	     box(
               title=uiOutput("scattertitle2"), width=NULL, height=250,
               withSpinner(plotOutput("scatterplots2"))
             )
	     
          )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
	my_db_file <- "Data/BloodCellType.sqlite"
	my_db <- dbConnect(RSQLite::SQLite(), my_db_file)
	#probeID <- "cg00000158"
	
	celltypes<-c("Buccal", "Nasal", "Whole.Blood", "B.cells","CD4.T.cells","CD8.T.cells","Granulocytes","Monocytes")
	officialnames<-c("Buccal","Nasal","Whole Blood", "B-cells","CD4T-cells","CD8T-cells","Granulocytes","Monocytes")
	names(officialnames)<-celltypes
  
  findprobebetas<-reactive({
    dat<-unlist(dbGetQuery(my_db, paste0("SELECT * FROM dnam WHERE row_names == ?"), params = c(input$probe))[-1])
	probebetas<-matrix(dat, ncol = length(celltypes))
	colnames(probebetas)<-c("Buccal","B.cells","CD4.T.cells","CD8.T.cells","Monocytes","Whole.Blood","Granulocytes","Nasal")
	probebetas<-as.data.frame(probebetas[,celltypes])
    return(probebetas)
  }) 
  
  summarydata<-reactive({
    tmp<-dbGetQuery(my_db, 'SELECT * FROM sumStats WHERE row_names == ?', params = c(input$probe))
    return(tmp)
  })
  
  
  cols<-rainbow(8) 
  names(cols)<-celltypes
  
  output$boxplot <- renderPlot({
    par(mar=c(2,4,1,1))
    probebetas<-findprobebetas()
    boxplot(probebetas, use.cols=FALSE, names=officialnames, ylab="DNA Methylation", main="", col=cols)
  
    }, height=300, width=800)
  
  output$boxplottitle<-renderUI({paste("Distribution of DNA Methylation levels at ",input$probe,"across cell & tissue types")})
  
 
}

# Run the application 
shinyApp(ui = ui, server = server)

