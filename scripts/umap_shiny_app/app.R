#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
rm(list = ls())
library("shiny")
library("tidyverse")
library("org.Hs.eg.db")
library("CORElearn")
library("e1071")

joined_tpm        <- readRDS("./data/joined_tpm_tibble.rds")
joined_meta       <- readRDS("./data/joined_meta.rds")
joined_tpm_numeric_indx <- unlist(lapply(joined_tpm,is.numeric))
joined_tpm_numeric <- joined_tpm[,joined_tpm_numeric_indx]
joined_tpm_numeric_variance <- apply(joined_tpm_numeric,2,var)


#entrez_to_ensembl   <- mapIds(org.Hs.eg.db,as.character(tx),"ENSEMBL","SYMBOL")

ui <- fluidPage(
  
  # Application title
  titlePanel("Met500: Feature Selection and UMAP"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectInput("dis",
                  label = "Select Disease:",
                  choices = sort(unique(joined_meta$disease)) ),
      
      
      selectInput("dataset",
                  label = "Select Dataset:",
                  choices = unique(joined_meta$dataset) )
      
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Summary Table", tableOutput("main_table")),
                  tabPanel("VolcanoPlot")
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  
  dis_filter <- reactive({joined_meta$disease == input$dis})
  dataset_filter <- reactive({joined_meta$dataset == input$dataset})
  
  fitered_tpm <- reactive({joined_tpm_numeric[dis_filter() & dataset_filter(),]})
  filtered_meta <- reactive({joined_meta[dis_filter() & dataset_filter(),]})
  anova_list <- reactive({})
  

   p_vals <- reactive({rpeat_anova(filtered_meta,c(1,2))})
   
  output$main_table <- renderTable({
    filtered_meta()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

