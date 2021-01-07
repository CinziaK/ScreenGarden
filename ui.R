#ui
library(tidyverse)
library(lubridate)
library(rlang)
library(ggplot2)
library(Cairo)
library(gghighlight)
library(shiny)
library(shinythemes)
library(png)
# Define UI
ui <- fluidPage(theme = shinytheme("flatly"),
                navbarPage(
                  # theme = "cerulean",  # <--- To use a theme, uncomment this
                  "ScreenGarden",
                  tabPanel("Home",
                           fluidRow(
                        
                               img(src = "SG2.png", height = "100%" , width = "100%")
                             ), 
                           
                           hr(),
                           fluidRow(h5("ScreenGarden is designed to facilitate the analysis of plate-based
                                       high-throughput growth assay with one or two different 
                                       controls in 96,384 or 1536 colony format. For a detailed 
                                       description of each analysis step, please download the 
                                       instructions using the download button below.
                                       Enjoy the screengardening!")),
                           
                           fluidRow(downloadButton("downloadpdf", "Download Instructions"))
                           ),
                                     
                  tabPanel("CalculateLGRs",
                           sidebarPanel(
                             textInput("query", "Query Name:", ""),
                             textInput("control", "Control Name:", ""),
                             
                             # Horizontal line ----
                             tags$hr(),
                             
                             # Input: Select a file ----
                             fileInput("file1", "Choose colonyAreas.txt File",
                                       multiple = TRUE,
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv")),
                             fileInput("file2", "Choose Keyfile",
                                       multiple = TRUE,
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv")),
                             tags$hr(),
                             # Horizontal line ----
                             # Input: Checkbox if file has header ----
                             #checkboxInput("header", "Header", FALSE),
                             
                             # Input: Select separator ----
                             #radioButtons("sep", "Separator",
                             #choices = c( Tab = "\t"),
                             #selected = "\t"),
                             
                             radioButtons("replicates", "Replicates",
                                          choices = c( "4" = 4,
                                                       "16" = 16),
                                          selected = 4),
                             
                             radioButtons("array", "Plate array",
                                          choices = c("384" = 384,
                                                      "1536" = 1536),
                                          selected = 1536),
                             
                             checkboxInput("smooth", "Smoothing", TRUE),
                             
                             
                             # Input: Select number of rows to display ----
                             radioButtons("disp", "Display",
                                          choices = c(Means = "Means",
                                                      Replicates = "Replicates"),
                                          selected = "Means"),
                             
                             
                             #Download Button
                             downloadButton("downloadData", "Download")
                             
                           ),
                           
                           # Main panel for displaying outputs ----
                           mainPanel(
                             
                             # Output: Data file ----
                             tableOutput("contents")
                             
                           )
                  ),
                  
                  
                  # Navbar 2, tabPanel
                  tabPanel("Combine2controls",sidebarPanel(
                    
                    fileInput("file3", "Choose CTR1 File",
                              multiple = TRUE,
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
                    fileInput("file4", "Choose CTR2 File",
                              multiple = TRUE,
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
                    tags$hr(),
                    
                    #Download Button
                    downloadButton("downloadData2", "Download")
                    ,
                    
                    # Main panel for displaying outputs ----
                    mainPanel(
                      
                      # Output: Data file ----
                      tableOutput("contents2")
                      
                    ))
                  ),
                  tabPanel("Plots", sidebarPanel(
                    
                    fileInput("file5", "Choose CSV File",
                              multiple = TRUE,
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
                    tags$hr(),
                    # "Empty inputs" - they will be updated after the data is uploaded
                    selectInput('xcol', 'X Variable Plot1', ""),
                    selectInput('ycol', 'Y Variable Plot1', "", selected = ""),
                    
                    tags$hr(),
                    tags$hr(),
                    
                    
                    sliderInput(inputId = "bins",
                                label = "Number of bins in histogram of data distribution:",
                                min = 1,
                                max = 100,
                                value = 50)
                    
                  ),
                  mainPanel(
                    fluidRow(
                      column(width = 10, class = "well",
                             h4("Data plot"),
                             plotOutput("plot1", height = 500,
                                        dblclick = "plot1_dblclick",
                                        brush = brushOpts(
                                          id = "plot1_brush",
                                          resetOnNew = TRUE
                                        )
                             )
                      )
                    ),
                    downloadButton(outputId = "down", label = "Download the plot"),
                    plotOutput(outputId = "distPlot")
                  )
                  ),
                  tabPanel("Mixture Model", sidebarPanel(
                    
                    fileInput("file6", "Choose CSV File",
                              multiple = TRUE,
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
                    tags$hr(),
                    downloadButton("downloadData3", "Download")
                  ),
                  
                  mainPanel(
                    fluidRow(
                      splitLayout(cellWidths = c("50%", "50%"), plotOutput("fitplot"), plotOutput("componentplot")),
                      tableOutput("contents3")
                      
                    )
                  )),
                  tabPanel("SPIalyser", sidebarPanel(
                    
                    fileInput("file7", "Choose CSV File",
                              multiple = TRUE,
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
                    tags$hr(),
                    downloadButton("downloadData4", "Download")
                  ),
                  
                  mainPanel(
                    fluidRow(
                      column(width = 12, class = "well",
                             h4("Cellular localisation of SPIs (LGR > 0.4)"),
                             plotOutput("donut", height = 400
                                        )
                             )
                      )
                    )
                  )
                
                ) # navbarPage
) # fluidPage

