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
                  "ScreenGarden",
                  tabPanel("Home", #First tab is homepage, with Logo and description
                           fluidRow(
                        
                               img(src = "SG2.png", height = "100%" , width = "100%") 
                             ), 
                           
                           hr(),
                           fluidRow(h5("ScreenGarden is designed to facilitate the analysis of plate-based
                                       high-throughput growth assay with one or two different 
                                       controls in 384 or 1536 colony format. For a detailed 
                                       description of each analysis step, please download the 
                                       instructions using the download button below.
                                       Enjoy the screengardening!")),
                           
                           fluidRow(downloadButton("downloadpdf", "Download Instructions"), # download the instructions to run ScreenGarden
                                    downloadButton("downloadzip", "Download R scripts"), #download the scripts, for user who wants to run ScreenGarden directly from R studio
                                    downloadButton("downloadzip2", "Download example files")) #download the scripts, for user who wants to run ScreenGarden directly from R studio
                                    
                           ),
                                     
                  tabPanel("CalculateLGRs", # second tab to calculate LGRs compared to one control
                           sidebarPanel(
                             textInput("query", "Query Name:", ""), # field to enter query name
                             textInput("control", "Control Name:", ""), # field to enter control name
                             
                             # Horizontal line ----
                             tags$hr(),
                             
                             # Input: Select a file ----
                             fileInput("file1", "Choose Colony Size File", # field to enter colony size input file such as cm engine output (check format if error message)
                                       multiple = TRUE,
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv")),
                             fileInput("file2", "Choose Keyfile", # field to enter keyfile (check format if error message)
                                       multiple = TRUE,
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv")),
                             tags$hr(),
                             # Horizontal line ----
                             
                             radioButtons("filetype", "Which software was used to measure colony size?", # select input file format
                                          choices = c( "CM engine" = 1, 
                                                       "other" = 2),
                                          selected = 1),
                             radioButtons("replicates", "Replicates", # select number of replicates 
                                          choices = c( "4" = 4,
                                                       "16" = 16,
                                                       "1" = 1),
                                          selected = 4),
                             
                             radioButtons("array", "Plate array", # select number of colonies on plate
                                          choices = c("384" = 384,
                                                      "1536" = 1536),
                                          selected = 1536),
                             
                             radioButtons("correct", "Plate correction method", # select normalisation method (based on median or specific positive controls)
                                          choices = c("Plate median" = 7,
                                                      "Mean Positive Control" = 8),
                                          selected = 7),
                             
                             checkboxInput("smooth", "Smoothing", TRUE), # smoothing should be selected for screens that use median-correction as correction method, when most colonies are not deficient in growth
                             
                             
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
                  
                  
                  # Navbar 2, tabPanel for combining the analyses of two independent controls
                  tabPanel("Combine2controls",sidebarPanel(
                    
                    fileInput("file3", "Choose CTR1 File", # select output files from previous tab 
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
                  # plot data for easy and quick analysis 
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
                    
                    #plot histogram to show data distribution
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
                  # Mixture model tab to define q-value thresholds for screens that can be fitted into a bimodal distribution
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
                  ))
                
                ) # navbarPage
) # fluidPage

