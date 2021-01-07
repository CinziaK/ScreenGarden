library(tidyverse)
library(lubridate)
library(rlang)
library(ggplot2)
library(Cairo)
library(gghighlight)
library(shiny)
library(shinythemes)
library(mclust)

l <- read.csv("data/loc.csv")
# Define server function  
shinyServer(function(input, output, session) {
  
  output$downloadpdf <- downloadHandler(
    filename = "SG.pdf",
    content = function(file) {
      file.copy("www/SG.pdf", file)
  })
  
  getData <- reactive({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    req(input$file2)
    req(query <- unlist(strsplit(input$query, ",")))
    req(control <- unlist(strsplit(input$control, ",")))
    
    dat <- read.csv(input$file1$datapath,
                     header = F,
                     sep = "\t")
    
    if (input$array == 1536) {
      (rn = 32) & (cn = 48)
    } else {
      if (input$array == 384) {
        (rn = 16) & (cn = 24)
      } else {
        stop('invalid array size')}}
    
    
    data = dat[c(1)]
    i =  nrow(data)/(1536+1)
    data = as.data.frame(split(data, 1:(1536+1)))
    rownames(data) <- data[,1]
    data = data[c(2:(1536+1))]
    
    
    if (input$array == 1536){
      data2 = as.tibble(t(data))
    }else{
      data2 = as.tibble(t(sdata))
    }

    Row = rep(c(1:(rn)),times=(cn))
    Colum = rep(c(1:cn), each =((rn)))
    
    data3 = mutate(data2, Row, Colum)
    data3 = as.tibble(data3)
    for ( col in 1:ncol(data3)){
      colnames(data3)[col] <-  sub(",.tif*", "", colnames(data3)[col])
    }
    
    data4 = gather(data3[c(1:i)], key = Plasmid, value = colonysize, na.rm = FALSE, convert = FALSE)
    
    d = nrow(data4)/(1536/cn)
    Row = rep(c(1:(1536/cn)),times=d)
    Column = rep(c(1:cn), each =(32))
    Column = rep(Column, times = i)
    
    data5 = mutate(data4, Row, Column)
    data6 = separate( data5, col = 1, into = c("Plasmid", "Plate"), sep = ",")
    data6$Plate = as.numeric(as.character(data6$Plate))
    
    
    key = read.csv(input$file2$datapath,
                   header = T,
                   sep = "\t")
    
    key
    key2 = mutate(key, Row = (key$Row * 2-1), Column = (key$Column * 2-1))
    key3 = mutate(key2, Row2 = (key$Row * 2), Column2 = (key$Column * 2))
    key4 = mutate(key3, Row3 = (key$Row), Column3 = (key$Column * 2))
    key5 = mutate(key4, Row4 = (key$Row * 2), Column4 = (key$Column))
    key6 = key3 %>% gather(n, Row, Row,Row2, na.rm = FALSE, convert = FALSE)
    key7 = key6 %>% gather(n, Column, Column,Column2,na.rm = FALSE, convert = FALSE)
    key7 = key7[c(1:3,5)]
    
    key8 = mutate(key7, Row = (key7$Row * 2-1), Column = (key7$Column * 2-1))
    key9 = mutate(key8, Row2 = (key7$Row * 2), Column2 = (key7$Column * 2))
    key10 = key9 %>% gather(n, Row, Row,Row2, na.rm = FALSE, convert = FALSE)
    key11 = key10 %>% gather(n, Column, Column,Column2,na.rm = FALSE, convert = FALSE)
    key12 = key11[c(1:3,5)]
    
    if (input$replicates == 4){
      colAreas = inner_join(data6, key7, by = c("Plate", "Row", "Column"))
    }else{
      if (input$replicates == 16){
        colAreas = inner_join(data6, key12, by = c("Plate", "Row", "Column"))
      }
    } 
    
    colAreas$colonysize = as.numeric(as.character(colAreas$colonysize))
    colAreas = colAreas[c(2,4,5,1,3,6)]
    
    x2 = mutate(colAreas, colonysize + 1) # no 0 Colonysizes, important for calculations
    x3 <- select (x2,-c(colonysize)) # filters Colony size column out
    #xnospace = apply(x3[1],2,function(x)gsub('\\s+', '',x))
    #x4 = x3$Plasmid <- xnospace
    #x5 = transform(x3, x4)
    
    x4 = reactive({
      if ((input$query == "" %in% x3$Plasmid) & (input$control == "" %in% x3$Plasmid)){
        return (x4 = filter(x3, Plasmid %in% c(input$query, input$control)))
      } else {
        stop('invalid name of query or control')}
    })# filters plasmids which are not used for comparison! Thus, results.tab can contain multiple plasmids
    x4 = reactive({
      if ((input$control == "" %in% x3$Plasmid) & (input$query == "" %in% x3$Plasmid)){
        return (x4 = filter(x3, Plasmid %in% c(input$control, input$query)))
      } else {
        stop('invalid name of query or control')}
    })# filters plasmids which are not used for comparison! Thus, results.tab can contain multiple plasmids
    x4 = filter(x3, Plasmid %in% c(input$query, input$control))
    
    xrearranged <- spread (x4,Plasmid,"colonysize + 1")
    colnames(xrearranged) = gsub(" ","",colnames(xrearranged))
    colnames(xrearranged)[which(colnames(xrearranged)==query)] = "query"
    colnames(xrearranged)[which(colnames(xrearranged)==control)] = "control" 
    LGR_GBP = c() # necessary later (don't know why its here actually but it works)
    
    platenumber = c(max(xrearranged$Plate))
    querymedian = c()
    controlmedian = c()     #make vector functions for medians -> advantage is that table isnt affected
    for(j in 1:platenumber){
      tempMedian = median(as.numeric(filter(xrearranged, Plate==j, !ORF== "BLANK")$query))
      querymedian = c(querymedian, rep(tempMedian, nrow(filter(xrearranged, Plate==j))))
      tempMedian = median(as.numeric(filter(xrearranged, Plate==j, !ORF== "BLANK")$control))
      controlmedian = c(controlmedian, rep(tempMedian, nrow(filter(xrearranged, Plate==j))))
    }
    #calculate median based on different plates (filter for 1:12) and exclude blank ORFs as they would impact media, claculate for all plasmids
    
    xmedian = mutate(xrearranged, querymedian, controlmedian)  #insert medians in new columns with mutate function
    xmedianfiltered = filter(xmedian, !ORF== "BLANK") # filter median for blank ORFs! controls are non-growing, thus median and smoothing would be incorrect if included
    mediancorrected <- transform(xmedianfiltered, norm_query = query/querymedian,  norm_control = control/controlmedian)
    # divide coloniesizes by plate median
    
    #################################################################################
    # claculate logarithm
    log = mutate(mediancorrected, log_norm_query = log(norm_query), log_norm_control = log(norm_control))
    
    LGR = c()
    LGR <- (log$log_norm_control- log$log_norm_query)
    xLGR = mutate(log, LGR)
   ###############
       if(input$smooth) {
      xLGR_S = xLGR
      
      xLGR_S = xLGR_S[order(xLGR_S$Plate),]
      
      platemedian = c()
      for(j in 1:platenumber){
        tempMedian2 = median(as.numeric(filter(xLGR_S, Plate==j)$LGR))
        platemedian = c(platemedian, rep(tempMedian2, nrow(filter(xLGR_S, Plate==j))))
      }
        
      xLGR_S = mutate(xLGR_S, platemedian)
        
      xLGR_S = xLGR_S[order(xLGR_S$Plate, xLGR_S$Row),]
      rowmean = c()
      for(j in 1:platenumber) { for (w in 1:32) {
        tempmean2 = mean(as.numeric(filter(xLGR_S, Plate==j, Row==w)$LGR))
        rowmean = c(rowmean, rep(tempmean2, nrow(filter(xLGR_S, Plate==j, Row==w))))
      }}
      
      xLGR_rsmoothed = mutate(xLGR_S, rowmean)
      
      xLGR_rsmoothed2  <- transform(xLGR_rsmoothed, LGRs = (LGR - (rowmean - platemedian)))
      
      platemedian2 = c()
      for(j in 1:platenumber){
        tempMedian2 = median(as.numeric(filter(xLGR_rsmoothed2, Plate==j)$LGRs))
        platemedian2 = c(platemedian2, rep(tempMedian2, nrow(filter(xLGR_rsmoothed2, Plate==j))))
      }
      xLGR_rsmoothed2 = mutate(xLGR_rsmoothed2, platemedian2)
      
      xLGR_rsmoothed2 = xLGR_rsmoothed2[order(xLGR_rsmoothed2$Plate, xLGR_rsmoothed2$Column),]
      
      columnmean = c()
      for(j in 1:platenumber) { for (h in 1:48) {
        tempmean3 = mean(as.numeric(filter(xLGR_rsmoothed2, Plate==j, Column==h)$LGRs))
        columnmean = c(columnmean, rep(tempmean3, nrow(filter(xLGR_rsmoothed2, Plate==j, Column==h))))
      }}
      xLGR_rcsmoothed2 = mutate(xLGR_rsmoothed2, columnmean)
      xLGR_rcsmoothed2  <- transform(xLGR_rcsmoothed2, LGRs = (LGRs - (columnmean - platemedian2)))
      
      xLGR_rcsmoothed2 = xLGR_rcsmoothed2[order(xLGR_rcsmoothed2$Plate, xLGR_rcsmoothed2$Row),]
      
      
      ############################## Levelplates #################################
      # plates are leveled to 0. Median scores are calculated for a given plate and adjusted to 0.
      
      platemedian3 = c()
      for(j in 1:platenumber){
        tempMedian3 = median(as.numeric(filter(xLGR_rcsmoothed2, Plate==j)$LGRs))
        platemedian3 = c(platemedian3, rep(tempMedian3, nrow(filter(xLGR_rcsmoothed2, Plate==j))))
      }
      xLGR_rcsmoothed2 = mutate(xLGR_rcsmoothed2, platemedian3)
      xLGR <- transform(xLGR_rcsmoothed2, LGRs = (LGRs - platemedian3))
      xLGR = select (xLGR,-c(rowmean, columnmean, platemedian, platemedian2, platemedian3))
    }
    
    ################################################################
    # calculate meanRow and meanColumn to combine replicates
    xLGR2 <- as.tibble(xLGR [,1:3])
    
    
    if (input$replicates == 4){
      xLGR2 <- mutate(xLGR2, meanRow = ceiling(as.numeric(Row)/2), meanColumn = ceiling(as.numeric(Column)/2))
    }else{ 
      if(input$replicates == 16){
        xLGR2 <- mutate(xLGR2, meanRow = ceiling(as.numeric(Row)/4), meanColumn = ceiling(as.numeric(Column)/4))}
    }
    
    
    xLGR3 <- inner_join(xLGR, xLGR2, by = c("Plate", "Row", "Column"))
    xLGR3 = transform(xLGR3, Plate = as.numeric(Plate))
    xLGR3 = transform(xLGR3, Plate = as.numeric(Plate))
    xLGR3 = arrange(xLGR3, Plate, meanRow, meanColumn)
    ###############################################################
    # organising the single replictates in rows
    if (input$replicates == 4) {
      t = rep(c(1,2,3,4),times=(nrow(xLGR3)/4))
    }else{
      if (input$replicates == 16) {
        t = rep(c(1:16),times=(nrow(xLGR3)/16))
      }
    }
    
    xt = mutate (xLGR3,t)
    
    if(input$smooth){
      xtreduced <- xt [c(1,4:6,9:17)]
    }else{
      xtreduced <- xt [c(1,4:6,9:16)]}
    
    if(input$smooth){
      xtreduced1 = gather(xtreduced,  'x', value, query, control, norm_query, norm_control, log_norm_query, log_norm_control, LGR, LGRs)
    }else{
      xtreduced1 = gather(xtreduced,  'x', value, query, control, norm_query, norm_control, log_norm_query, log_norm_control, LGR)
    }
    
    xtreduced2 = mutate(xtreduced1, Plasmid_replicate_stats = paste(x, t))
    xtreduced3 = xtreduced2 [c(1:4,7:8)]
    y = spread (xtreduced3, Plasmid_replicate_stats, value)
    
    #################################################################
    # Z-scores for each replicate
    #LGR1 = "LGR 1"
    #LGR2 = "LGR 2"
    #LGR3 = "LGR 3"
    #LGR4 = "LGR 4"
    #which(colnames(y)==LGR1)
    #y1 = y[c(which(colnames(y)==LGR1),which(colnames(y)==LGR2),which(colnames(y)==LGR3),which(colnames(y)==LGR4))]
    #y1mean = mean(colMeans(y1))
    #y1sd = mean(apply(y1,2,sd))
    
    #Z_score_1 = c()
    #Z_score_2 = c()
    #Z_score_3 = c()
    #Z_score_4 = c()
    #Z_score_1 = (y$"LGR 1"-y1mean)/y1sd
    #Z_score_2 = (y$"LGR 2"-y1mean)/y1sd
    #Z_score_3 = (y$"LGR 3"-y1mean)/y1sd
    #Z_score_4 = (y$"LGR 4"-y1mean)/y1sd
    
    #yz = mutate(y,Z_score_1,Z_score_2,Z_score_3,Z_score_4)
    
    ####################
    ####calculate mean and standard deviation of replictaes####
    if(input$smooth){
      aggregatedmean <-  aggregate(xLGR3 [,5:16], by = list(xLGR3$Plate, xLGR3$meanRow, xLGR3$meanColumn), FUN=mean, na.rm=TRUE)
    }else{
      aggregatedmean <-  aggregate(xLGR3 [,5:15], by = list(xLGR3$Plate, xLGR3$meanRow, xLGR3$meanColumn), FUN=mean, na.rm=TRUE)
    }
    if(input$smooth){
      aggregatedmean = rename(aggregatedmean, mean_usLGR = LGR, Mean_LGR = LGRs, Plate = Group.1, Row = Group.2, Column = Group.3)
    }else{
      aggregatedmean = rename(aggregatedmean, Mean_LGR = LGR, Plate = Group.1, Row = Group.2, Column = Group.3)
    }
    
    if(input$smooth){
      aggregatedsd <- aggregate(xLGR3 [,5:16], by = list(xLGR3$Plate, xLGR3$meanRow, xLGR3$meanColumn), FUN=sd, na.rm=TRUE)
    }else{
      aggregatedsd <- aggregate(xLGR3 [,5:15], by = list(xLGR3$Plate, xLGR3$meanRow, xLGR3$meanColumn), FUN=sd, na.rm=TRUE)
    }
    if(input$smooth){
      aggregatedsd = rename(aggregatedsd, sd_usLGR = LGR, sd_LGR = LGRs, Plate = Group.1, Row = Group.2, Column = Group.3)
    }else{
      aggregatedsd = rename(aggregatedsd, sd_LGR = LGR, Plate = Group.1, Row = Group.2, Column = Group.3)
    }
    
    if(input$smooth){
      aggregatedsd <- aggregatedsd [c(1:3,12:13)]
    }else{
      aggregatedsd <- aggregatedsd [c(1:3,12)]
    }
    
    aggregated <- inner_join(aggregatedmean, aggregatedsd, by = c("Plate", "Row", "Column"))
    
    if(input$smooth){
      aggregated <- aggregated [c(1:13,16:17)]
    }else{
      aggregated <- aggregated [c(1:12,15)]
    }
    aggregated1 <- transform(aggregated, as.numeric(as.character(Plate)))
    
    ####### calculate Z-score over mean of replicates
    mean_Z_score = c()
    mean_Z_score = as.vector(scale(aggregated1$Mean_LGR))
    aggregated2 = mutate(aggregated1, mean_Z_score)
    
    ########## if you ever reading this: you're almost done! so have a glass of prosecco and celebrate! rn it's friday, half past 5 and i want to go home...
    #keyfile = read_delim("keyfile.txt", delim="\t")
    #system(paste("perl", "keyfile_letters_to_numbers.pl"))
    
    
    key = read.csv(input$file2$datapath,
                   header = T,
                   sep = "\t")
    key = rename(key, Plate = "Plate")
    key = transform(key, Plate = as.numeric(Plate))
    
    mean_file <- inner_join(aggregated2, key, by = c("Plate", "Row", "Column"))
    
    if(input$smooth){
      mean_file = mean_file[c(1:5,8:17)]
    }else{
      mean_file = mean_file[c(1:5,8:15)]
    }
    
    
    
    if(input$disp == "Means") {
      return(mean_file)
    } else {
      return(y)
    }
  }
  
  )
  output$contents <- renderTable(
    
    getData()
    
  )
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    
    filename = function() { 
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    
    content = function(file) {
      
      write.csv(getData(), file)
    }
  )
  
  getData2 <- reactive({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file3)
    req(input$file4)
    
    ctr1 <- read.csv(input$file3$datapath,
                     header = T)
    ctr2 <- read.csv(input$file4$datapath,
                     header = T)
    
    Merge.file <- inner_join(ctr1, ctr2, by = c("Plate", "Row", "Column", "ORF"), suffix = c(".1",".2"))
    Mean_LGR = c() # calculates the mean LGR of both comparisons
    Mean_LGR = rowMeans(Merge.file[c("Mean_LGR.1", "Mean_LGR.2")], na.rm=TRUE)
    Mean_Z = c() # calculates the mean Z score of both comparisons
    Mean_Z = rowMeans(Merge.file[c("mean_Z_score.1", "mean_Z_score.2")], na.rm=TRUE)
    
    Merge.File = mutate( Merge.file, Mean_LGR, Mean_Z)
    
    
    return(Merge.File)
  })
  output$contents2 <- renderTable(
    
    getData2()
    
  ) 
  
  # Downloadable csv of selected dataset ----
  output$downloadData2 <- downloadHandler(
    
    filename = function() { 
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    
    content = function(file) {
      
      write.csv(getData2(), file)
    }
  )
  
  
  data <- reactive({
    
    req(input$file5)
    
    
    A <- read.csv(input$file5$datapath,
                  header = T)
    A = as.data.frame(arrange(A, desc(Mean_LGR)))
    n = c(1:nrow(A))
    A = mutate (A, n)
    
    updateSelectInput(session, inputId = 'xcol', label = 'X Variable Plot1',
                      choices = names(A), selected = names(A))
    updateSelectInput(session, inputId = 'ycol', label = 'Y Variable Plot1',
                      choices = names(A), selected = names(A)[2])
    return(A)
    
  })
  
  
  plot1 <- reactive({
    ranges <- reactiveValues(x = NULL, y = NULL)
    # output$plot1 <- renderPlot({
    p  = data()
    
    #if(is.null(plot.x)) return()
    #if(plot.x$xcol == "" | plot.x$ycol =="") return()
    plot1 <- ggplot(p,
                    aes_string(
                      x =  input$xcol,
                      y = input$ycol
                    ))+ geom_point()+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y ,  expand = TRUE)+
      gghighlight(Mean_LGR > 0.4) + geom_hline(yintercept = c(0.4,-0.4),colour="red", linetype="dashed")
    
    print(plot1)
  })
  
  output$plot1 <- renderPlot({
    plot1()
  })
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  #observeEvent(input$plot1_dblclick, {
  # brush <- input$plot1_brush
  #if (!is.null(brush)) {
  # ranges$x <- c(brush$xmin, brush$xmax)
  #ranges$y <- c(brush$ymin, brush$ymax)
  
  #} else {
  # ranges$x <- NULL
  #ranges$y <- NULL
  
  
  # }
  
  #})
  
  output$down <- downloadHandler(
    filename = function(){paste(input$data, '.pdf', sep = '')},
    content = function(file){
      ggsave(file, plot1())
    },
    contentType = "application/pdf"
  )
  
  
  output$distPlot <- renderPlot({
    B =  data()
    x    <- B$Mean_LGR
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    hist(x, breaks = bins, col = "#75AADB", border = "white",
         xlab = "Mean LGR", 
         ylab = "Number of strains",
         main = "Distribution of LGRs")
  })
  MMdata <- reactive({
    
    req(input$file6)
    
    MM <- read.csv(input$file6$datapath,
                   header = T)
    
    select = dplyr::select
    MM2 = MM %>% select(ORF, Mean_LGR)
    
    modG2 = densityMclust(MM2$Mean_LGR, G = 2)
    # The order of the fit models is not always such that m_1<m_2 so define ordering based on this
    Order = sort(modG2$parameters$mean, index.return = T)$ix
    maxG = 2
    minG = 1
    means = modG2$parameters$mean[Order]
    pros = modG2$parameters$pro[Order]
    sds = sqrt(modG2$parameters$variance$sigmasq[Order])
    
    # Define whether the model fitting is successful based on whether the means are within 1.5 sd of each other
    failure = 0
    if((means[2]- means[1])<1.5*sds[1]){
      failure = 1
    }
    print(failure)
    if(failure==0){
      
      # Define distributions: f1 is mixture model, f2 is hit peak component
      f1 = function(x) modG2$parameters$pro[minG]*dnorm(x, mean = modG2$parameters$mean[minG], sd = sqrt(modG2$parameters$variance$sigmasq[minG]))+modG2$parameters$pro[maxG]*dnorm(x, mean = modG2$parameters$mean[maxG], sd = sqrt(modG2$parameters$variance$sigmasq[maxG]))
      f2 = function(x) modG2$parameters$pro[maxG]*dnorm(x, mean = modG2$parameters$mean[maxG], sd = sqrt(modG2$parameters$variance$sigmasq[maxG]))
      
      # Work out q(x) and write data
      
      q = function (x) f2(x)/f1(x)
      
      MMdata = mutate(MM2, qval = q(Mean_LGR))
      
      return(MMdata)
    } else {
      renderPrint("Mixture model fitting failed")}
    
  })
  
  
  output$contents3 <- renderTable(
    MMdata())
  
  output$downloadData3 <- downloadHandler(
    
    filename = function() { 
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    
    content = function(file) {
      
      write.csv(MMdata(), file)
    }
  )
  
  #Plots
  fitplot <- reactive({
    
    MM = MMdata()
    modG2 = densityMclust(MM$Mean_LGR, G = 2)
    # The order of the fit models is not always such that m_1<m_2 so define ordering based on this
    Order = sort(modG2$parameters$mean, index.return = T)$ix
    maxG = 2
    minG = 1
    means = modG2$parameters$mean[Order]
    pros = modG2$parameters$pro[Order]
    sds = sqrt(modG2$parameters$variance$sigmasq[Order])
    # Define distributions: f1 is mixture model, f2 is hit peak component
    f1 = function(x) modG2$parameters$pro[minG]*dnorm(x, mean = modG2$parameters$mean[minG], sd = sqrt(modG2$parameters$variance$sigmasq[minG]))+modG2$parameters$pro[maxG]*dnorm(x, mean = modG2$parameters$mean[maxG], sd = sqrt(modG2$parameters$variance$sigmasq[maxG]))
    f2 = function(x) modG2$parameters$pro[maxG]*dnorm(x, mean = modG2$parameters$mean[maxG], sd = sqrt(modG2$parameters$variance$sigmasq[maxG]))
    
    # Work out q(x) and write data
    
    q = function (x) f2(x)/f1(x)
    
    x = seq(min(MM$Mean_LGR),max(MM$Mean_LGR),length = 1000)
    dist = f1(x)
    data = hist(MM$Mean_LGR, breaks = 150, plot = FALSE)
    distVals = tibble(X = x, Y = dist, Distribution = "Mixture model")
    dataVals = tibble(X = data$mids, Y = data$density, Distribution = "Screen data")
    Vals = rbind(distVals, dataVals)
    #png(filename = paste0(name,"MMfit.png"), width = 600)
    fitplot <- ggplot(Vals, aes(x= X, y = Y, colour = Distribution))+ geom_line(linetype = 1) +  theme_bw() +scale_color_brewer(palette="Set1")+labs(x="LGR",y="density")+theme(text = element_text(size = 12,family = "Helvetica")) 
    print(fitplot)
  })
  
  componentplot <- reactive({
    MM = MMdata()
    modG2 = densityMclust(MM$Mean_LGR, G = 2)
    # The order of the fit models is not always such that m_1<m_2 so define ordering based on this
    Order = sort(modG2$parameters$mean, index.return = T)$ix
    maxG = 2
    minG = 1
    means = modG2$parameters$mean[Order]
    pros = modG2$parameters$pro[Order]
    sds = sqrt(modG2$parameters$variance$sigmasq[Order])
    # Define distributions: f1 is mixture model, f2 is hit peak component
    f1 = function(x) modG2$parameters$pro[minG]*dnorm(x, mean = modG2$parameters$mean[minG], sd = sqrt(modG2$parameters$variance$sigmasq[minG]))+modG2$parameters$pro[maxG]*dnorm(x, mean = modG2$parameters$mean[maxG], sd = sqrt(modG2$parameters$variance$sigmasq[maxG]))
    f2 = function(x) modG2$parameters$pro[maxG]*dnorm(x, mean = modG2$parameters$mean[maxG], sd = sqrt(modG2$parameters$variance$sigmasq[maxG]))
    
    # Work out q(x) and write data
    
    q = function (x) f2(x)/f1(x)
    x = seq(min(MM$Mean_LGR),max(MM$Mean_LGR),length = 1000)
    dist = f1(x)
    comp = tibble(X =x, Y = modG2$parameters$pro[1]*dnorm(x, mean = modG2$parameters$mean[1], sd = sqrt(modG2$parameters$variance$sigmasq[1])), Distribution = "Distribution 1")
    int= tibble(X=x, Y = modG2$parameters$pro[2]*dnorm(x, mean = modG2$parameters$mean[2], sd = sqrt(modG2$parameters$variance$sigmasq[2])), Distribution = paste("Distribution", 2))
    comp = rbind(comp, int)
    comp = rbind(comp, tibble(X=x, Y=dist, Distribution = "Mixture model"))
    # png(filename = paste0(name,"MMcomponents.png"), width = 600)
    componentplot <- ggplot(comp, aes(x = X, y = Y, colour = Distribution)) + geom_line() +  theme_bw() +scale_color_brewer(palette="Set1")+labs(x="LGR",y="density")+theme(text = element_text(size = 12,family = "Helvetica"),legend.background = element_rect(fill="transparent"))
    print(componentplot)
  })
  
  output$fitplot <- renderPlot({
    fitplot()
  })
  output$componentplot <- renderPlot({
    componentplot()
  })
  
  donutdata <- reactive({
    
    req(input$file7)
    
    
    D <- read.csv(input$file7$datapath,
                  header = T)
  
    x0 = aggregate(l$sum, by=list(loc=l$loc), FUN=sum)
    x0 = filter(x0, loc != "n/a")
    
    
    tx = inner_join(D,l, by = "ORF")
    tx = tx[!duplicated(tx$ORF), ]
    x2 = filter(tx, Mean_LGR > 0.4)
    x2 = aggregate(x2$sum, by=list(loc=x2$loc), FUN=sum)
    x2
    
    x2 = filter(x2, loc != "n/a")
    x3 = inner_join(x2,x0, by = "loc")
    v0 = as.vector(x3$x.y)
    v2 = as.vector(x3$x.x)
    sumv0 = sum(v0)
    sumv2 = sum(v2)
    rv0 = v0/sumv0
    rv2 = v2/sumv2
    rich = rv2/rv0
    rich2 = format(round(rich,2), nsmall = 2)
    rich2
    x2 = mutate (x2, rich2)
    x2$fraction = x2$x / sum(x2$x)
    
    # Compute the cumulative percentages (top of each rectangle)
    x2$ymax = cumsum(x2$fraction)
    x2$ymin = c(0, head(x2$ymax, n=-1))
    
    d = x2  %>%  unite("Localisation", c(loc, rich2), sep = ", x", remove = FALSE)
    
    return(d)
    
  })
  
  
  donut <- reactive({
    
    d  = donutdata()
    
   donut = ggplot(d, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Localisation)) +
      geom_rect(stat="identity", color="white") + 
      scale_fill_brewer(palette="Set3") +
      scale_color_brewer(palette="Set3") +
      coord_polar(theta="y") + 
      xlim(c(2, 4)) +
      theme_void()
    
    print(donut)
  })
  
  output$donut <- renderPlot({
    donut()
  })
})


# Create Shiny object
#shinyApp(ui = ui, server = server)
