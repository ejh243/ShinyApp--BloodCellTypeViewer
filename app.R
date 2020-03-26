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
                textInput("probe", "EPIC array CpG ID:", value="cg00000158"),
				textOutput("downloadtext1"),
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
		
		checkboxInput("header","Header", value = FALSE, width = NULL),
		checkboxInput("rownames","Rownames", value = FALSE, width = NULL),
                downloadButton("downloaddata2","Download summary statistics for uploaded file")
            ),

             box(
               title = "Summary Stats Output", width = NULL,
		tableOutput('tbl')

            )

          ),

      column(width = 8,
             box(
               title=uiOutput("boxplottitle"), width = NULL, height=370,
               plotOutput("boxplot", height = "320px")
             ) ,            
             box(
               title=uiOutput("scattertitle1"), width=NULL, height=550,
               withSpinner(plotOutput("scatterplots1", height = "500px"))
             ) ,
	     box(
               title=uiOutput("scattertitle2"), width=NULL, height=300,
               withSpinner(plotOutput("scatterplots2", height = "250px"))
             )

            	     
          )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
	colInfo<-read.csv("Data/columnHeaders.csv", header = FALSE)
	colnames(colInfo)<-c("Header", "Contents")
	my_db_file <- "Data/BloodCellType.sqlite"
	my_db <- dbConnect(RSQLite::SQLite(), my_db_file)
	
	celltypes<-c("Buccal", "Nasal", "Whole.Blood", "B.cells","CD4.T.cells","CD8.T.cells","Granulocytes","Monocytes")
	officialnames<-c("Buccal","Nasal","Whole Blood", "B-cells","CD4T-cells","CD8T-cells","Granulocytes","Monocytes")
	names(officialnames)<-celltypes
  
  findprobebetas<-reactive({
   dat<-unlist(dbGetQuery(my_db, "SELECT * FROM dnam WHERE row_names == ?", params = c(input$probe))[-1])
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
    boxplot(probebetas, ylab = "DNA methylation", col = cols, names=officialnames)
  
    })
  
  output$boxplottitle<-renderUI({paste("Distribution of DNA Methylation levels at ",input$probe,"across cell & tissue types")})
  
 output$scatterplots1 <- renderPlot({
    par(mfrow=c(2,3))
    probebetas<-findprobebetas()
    for (each in c("B.cells","CD4.T.cells","CD8.T.cells","Granulocytes","Monocytes")){
		r<-signif(cor(probebetas[,"Whole.Blood"], probebetas[,each], use="complete"),3)
		plot(probebetas[,"Whole.Blood"], probebetas[,each],  xlab="Whole Blood", ylab = each, pch=20, col=cols[each], 
           cex.lab=1.4, cex.main=1.2, main=paste("r =",r))
	}
	})
	
  output$scattertitle1<-renderUI({paste("Co-variation between whole blood and blood cell types at ",input$probe)})
 
 output$scatterplots2 <- renderPlot({
    par(mfrow=c(1,3))
    probebetas<-findprobebetas()
    for (each in c("Buccal","Nasal")){
                r<-signif(cor(probebetas[,"Whole.Blood"], probebetas[,each], use="complete"),3)
                plot(probebetas[,"Whole.Blood"], probebetas[,each],  xlab="Whole Blood", ylab = each, pch=20, col=cols[each],
           cex.lab=1.4, cex.main=1.2, main=paste("r =",r))
        }
        })

  output$scattertitle2<-renderUI({paste0("Co-variation between peripheral tissues at ",input$probe)})

output$downloadtext1 <- renderText({paste("Download summary statistics for CpG:", input$probe)})
  
  output$downloaddata1 <- downloadHandler(
    filename = function(){
      paste(input$probe, ".csv", sep="")},
    
    content = function(file){
      write.csv(summarydata(), file, row.names = FALSE)},
    
    contentType = "text/csv"
    
  )

batchquery<-reactive({
if(input$rownames == TRUE){
	rowN=1
} else {
	rowN=NULL
}
df <- read.csv(input$file1$datapath,header = input$header, row.names = rowN)
        tmp<-dbGetQuery(my_db, 'SELECT * FROM sumStats WHERE row_names == ?',
        params = list(df[,1]))
	return(tmp)
})

  output$downloaddata2 <- downloadHandler(
	
   filename = function(){
      paste("sumstats", input$file1$name, sep="_")},

    content = function(file){
      write.csv(batchquery(),file, row.names = FALSE)},

    contentType = "text/csv"

  )
output$tbl <- renderTable({ colInfo })

}

# Run the application 
shinyApp(ui = ui, server = server)

