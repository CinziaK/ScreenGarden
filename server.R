### ScreenGarden script for plate-based large scale data analysis ###

# This script compares colony sizes of experimental plates to one or two independent controls and can be used to plot data.
# The user can upload colony size files derived from ScreenMill's CM Engine or other software to calculate log growth ratios. 
# written by Cinzia Klemm and Rowan Howell in 2019/2020

# if ScreenGarden is run locally from RStudio, following packages have to be installed prior running the script:
# tidyverse, lubridate, rlang, ggplot2, Cairo, gghighlight, shiny, shinythemes, mclust
library(tidyverse) # load package
library(lubridate) # load package
library(rlang) # load package
library(ggplot2) # load package
library(Cairo) # load package
library(gghighlight) # load package
library(shiny) # load package
library(shinythemes) # load package
library(mclust) # load package

# Define server function  
shinyServer(function(input, output, session) {
   # attach ScreenGarden manual which can be downloaded from the homepage, a step-by-step guide how to run ScreenGarden analysis 
  output$downloadpdf <- downloadHandler(
    filename = "SG.pdf",
    content = function(file) {
      file.copy("www/SG.pdf", file)
  })
  # attach ScreenGarden R scripts to homepage, for users who want to run ScreenGarden locally from R studio (offline)
 output$downloadzip <- downloadHandler(
   filename = "SGcode.zip",
   content = function(file) {
   file.copy("www/SGcode.zip", file)
    })
  
 # attach ScreenGarden R scripts to homepage, for users who want to run ScreenGarden locally from R studio (offline)
 output$downloadzip2 <- downloadHandler(
   filename = "ExampleFiles.zip",
   content = function(file) {
     file.copy("www/ExampleFiles.zip", file)
   })
 
  getData <- reactive({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1) # insert colony size file
    req(input$file2) # insert keyfile defining the ID at each position
    req(query <- unlist(strsplit(input$query, ",")))
    req(control <- unlist(strsplit(input$control, ",")))
    
    # if CM Engine input file is used (filetype 1), the file has to be transformed to the right format
    # if other software is used to measure colony sizes (filetype 2), the data has to be manually pre-formatted (see SG manual)
    if (input$filetype == 1){
    dat <- read.csv(input$file1$datapath,
                     header = F,
                     sep = "\t")
    
    ###################
    ###error is wrong array size selected
    ldat = nrow(dat)
   if (ldat %% 1537 == 0 | ldat %% 385 == 0) {T
    }else{
    validate('ERROR: The CM engine output (colonyAreas.txt) file is not in the 1536 or 384 colony format. Please check your input file!
             
    If you are not using the CM engine output for data analysis, please select "other" in the left panel'
      )}
    ####################
    
    if (input$array == 1536) {
      (rn = 32) & (cn = 48)
    } else {
      if (input$array == 384) {
        (rn = 16) & (cn = 24)
      } else {
        stop('invalid array size')}}
    
    if (input$array == 1536){
    data = dat[c(1)]
    i =  nrow(data)/(1536+1)
    data = as.data.frame(split(data, 1:(1536+1)))
    rownames(data) <- data[,1]
    data = data[c(2:(1536+1))]
    }else{
      if (input$array == 384){
        data = dat[c(1)]
        i =  nrow(data)/(384+1)
        data = as.data.frame(split(data, 1:(384+1)))
        rownames(data) <- data[,1]
        sdata = data[c(2:(384+1))]
      }
    }
    
    
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
    
    data4 = gather(data3[c(1:i)], key = Label, value = colonysize, na.rm = FALSE, convert = FALSE)
    
    if (input$array == 1536){
      d = nrow(data4)/(1536/cn)
      Row = rep(c(1:(1536/cn)),times=d)
      Column = rep(c(1:cn), each =(32))
      Column = rep(Column, times = i)
    }else{
      if (input$array == 384){
        d = nrow(data4)/(384/cn)
        Row = rep(c(1:(384/cn)),times=d)
        Column = rep(c(1:cn), each =(16))
        Column = rep(Column, times = i)
      }
    }
    
    data5 = mutate(data4, Row, Column)
    data6 = separate( data5, col = 1, into = c("Label", "Plate"), sep = ",")
    data6$Plate = as.numeric(as.character(data6$Plate))
    }else{
      if (input$filetype == 2){
      data6 <- read.csv(input$file1$datapath,
                      header = T,
                      sep = "\t")
      }
    }
    #############################    
    #define correct structure of file
    Correct_Columns <- c('Label', 'Plate', 'colonysize', 'Row', 'Column')
    #read current structure of table
    Current_Columns <- colnames(data6)
    
    #Check whether CSV was correctly imported from Source
    shiny::validate(
      need(all(Correct_Columns %in% Current_Columns), "ERROR: Imported colonysize file has incorrect columnheaders, please correct headers to: 
           Label    Plate   colonysize    Row   Column
           This script is case sensitive.")
    )
    #############################
    # upload keyfile, a file that defines the ID of each specific plate position. 
    # The file does NOT need to have replicates, as the file is automatically transformed in the following steps.
    key = read.csv(input$file2$datapath,
                   header = T,
                   sep = "\t")
    ############################    
    #define correct structure of file
    Correct_Columns2 <- c('Plate', 'Row', 'Column', 'ID')
    #read current structure of table
    Current_Columns2 <- colnames(key)
    
    shiny::validate(
      need(all(Correct_Columns2 %in% Current_Columns2), "ERROR: Imported keyfile has incorrect columnheaders, please correct headers to: 
           Plate    Row   Column    ID
           This script is case sensitive.")
    )
    #############################
    # transform keyfile depending on the number of replicates
    key # one replicate
    key2 = mutate(key, Row = (key$Row * 2-1), Column = (key$Column * 2-1))
    key3 = mutate(key2, Row2 = (key$Row * 2), Column2 = (key$Column * 2))
    key4 = mutate(key3, Row3 = (key$Row), Column3 = (key$Column * 2))
    key5 = mutate(key4, Row4 = (key$Row * 2), Column4 = (key$Column))
    key6 = key3 %>% gather(n, Row, Row,Row2, na.rm = FALSE, convert = FALSE)
    key7 = key6 %>% gather(n, Column, Column,Column2,na.rm = FALSE, convert = FALSE)
    key7 = key7[c(1:3,5)] # four replicates
    
    key8 = mutate(key7, Row = (key7$Row * 2-1), Column = (key7$Column * 2-1))
    key9 = mutate(key8, Row2 = (key7$Row * 2), Column2 = (key7$Column * 2))
    key10 = key9 %>% gather(n, Row, Row,Row2, na.rm = FALSE, convert = FALSE)
    key11 = key10 %>% gather(n, Column, Column,Column2,na.rm = FALSE, convert = FALSE)
    key12 = key11[c(1:3,5)] # sixteen replicates
    # select keyfile depending on selected number of replicates
    if (input$replicates == 4){
      colAreas = inner_join(data6, key7, by = c("Plate", "Row", "Column"))
    }else 
      if (input$replicates == 16){
        colAreas = inner_join(data6, key12, by = c("Plate", "Row", "Column")) 
      } else 
      if (input$replicates == 1){
        colAreas = inner_join(data6, key, by = c("Plate", "Row", "Column"))
      }
      
    colAreas$colonysize = as.numeric(as.character(colAreas$colonysize))
    colAreas = colAreas[c(2,4,5,1,3,6)]
    # add +1 to colony size so colony size does not equal zero
    x2 = mutate(colAreas, colonysize + 1) 
    x3 <- select (x2,-c(colonysize)) # filters Colony size column out
    
    # check if query and control names are correct
    if (input$query %in% x3$Label){
    }else{
      validate('ERROR: Invalid name of query. Please correct!')}
    
    if (input$control %in% x3$Label){
    }else{
      validate('ERROR: Invalid name of control. Please correct!')}
    # select query and control based on input names in user interface
    # even if input files contain multiple experimental and control plates, only the ones defined as query and control are selected for following analysis
    x4 = reactive({
      if ((input$query == "" %in% x3$Label) & (input$control == "" %in% x3$Label)){
        return (x4 = filter(x3, Label %in% c(input$query, input$control)))
      }else{
        validate('invalid name of query')}
    })# filters plasmids which are not used for comparison! 
    x4 = reactive({
      if ((input$control == "" %in% x3$Label) & (input$query == "" %in% x3$Label)){
        return (x4 = filter(x3, Label %in% c(input$control, input$query)))
      }else{
        validate('invalid name of control')}
    })# filters plasmids which are not used for comparison! 
    x4 = filter(x3, Label %in% c(input$query, input$control))
    
    xrearranged <- spread (x4,Label,"colonysize + 1")
    colnames(xrearranged) = gsub(" ","",colnames(xrearranged))
    colnames(xrearranged)[which(colnames(xrearranged)==query)] = "query"
    colnames(xrearranged)[which(colnames(xrearranged)==control)] = "control" 
    LGR_GBP = c() 
    
    ### normalisation/correction based on plate median (7) or mean Control colony size (8)
    platenumber = c(max(xrearranged$Plate)) # define maximal plate number
    querymedian = c()
    controlmedian = c() #make vector functions for medians -> advantage is that table isnt affected
    querymean = c()
    controlmean = c() 
    
    if(input$correct == 7) {
    for(j in 1:platenumber){
        tempMedian = median(as.numeric(filter(xrearranged, Plate==j, !ID== "BLANK")$query))
        querymedian = c(querymedian, rep(tempMedian, nrow(filter(xrearranged, Plate==j))))
        tempMedian = median(as.numeric(filter(xrearranged, Plate==j, !ID== "BLANK")$control))
        controlmedian = c(controlmedian, rep(tempMedian, nrow(filter(xrearranged, Plate==j))))
    }
      } else {
        if(input$correct == 8) {
        for(j in 1:platenumber){
          tempMean = mean(as.numeric(filter(xrearranged, Plate==j, ID== "Control")$query))
          querymean = c(querymean, rep(tempMean, nrow(filter(xrearranged, Plate==j))))
          tempMean = mean(as.numeric(filter(xrearranged, Plate==j, ID== "Control")$control))
          controlmean = c(controlmean, rep(tempMean, nrow(filter(xrearranged, Plate==j))))
           }
        }
      }
    # filter BLANKs and add calculate median/controlmean-corrected colony sizes
    if(input$correct == 7) {
      xmedian = mutate(xrearranged, querymedian, controlmedian)
      xmedianfiltered = filter(xmedian, !ID== "BLANK") # filter median for blank IDs! controls are non-growing, thus median and smoothing would be incorrect if included
      mediancorrected <- transform(xmedianfiltered, norm_query = query/querymedian,  norm_control = control/controlmedian)
      # divide coloniesizes by plate median
    } else 
      if(input$correct == 8) {
        xmean = mutate(xrearranged, querymean, controlmean)
        xmeanfiltered = filter(xmean, !ID== "BLANK") 
        meancorrected <- transform(xmeanfiltered, norm_query = query/querymean,  norm_control = control/controlmean)
      }
    # divide coloniesizes by colony size mean
    # claculate logarithm
    if (input$correct == 7) {
      log = mutate(mediancorrected, log_norm_query = log(norm_query), log_norm_control = log(norm_control))
    } else {
      if(input$correct == 8) {
        log = mutate(meancorrected, log_norm_query = log(norm_query), log_norm_control = log(norm_control))
      }
    }
    #calculate LGRs  
    LGR = c()
    LGR <- (log$log_norm_control- log$log_norm_query)
    xLGR = mutate(log, LGR)
   ### smoothing based on colony position (Row, Column) on plate ###
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
      ### Levelplates ###
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
    
    ### calculate meanRow and meanColumn to combine replicates ###
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
    ### organising the single replictates in rows for replicates file ###
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
    ### calculate mean and standard deviation of replictaes for mean file ###
    if(input$smooth){
      aggregatedmean <-  aggregate(xLGR3 [,5:16], by = list(xLGR3$Plate, xLGR3$meanRow, xLGR3$meanColumn), FUN=mean, na.rm=TRUE)
    }else{
      aggregatedmean <-  aggregate(xLGR3 [,5:15], by = list(xLGR3$Plate, xLGR3$meanRow, xLGR3$meanColumn), FUN=mean, na.rm=TRUE)
    }
    if(input$smooth){
      aggregatedmean = dplyr::rename(aggregatedmean, mean_usLGR = LGR, Mean_LGR = LGRs, Plate = Group.1, Row = Group.2, Column = Group.3)
    }else{
      aggregatedmean = dplyr::rename(aggregatedmean, Mean_LGR = LGR, Plate = Group.1, Row = Group.2, Column = Group.3)
    }
    
    if(input$smooth){
      aggregatedsd <- aggregate(xLGR3 [,5:16], by = list(xLGR3$Plate, xLGR3$meanRow, xLGR3$meanColumn), FUN=sd, na.rm=TRUE)
    }else{
      aggregatedsd <- aggregate(xLGR3 [,5:15], by = list(xLGR3$Plate, xLGR3$meanRow, xLGR3$meanColumn), FUN=sd, na.rm=TRUE)
    }
    if(input$smooth){
      aggregatedsd = dplyr::rename(aggregatedsd, sd_usLGR = LGR, sd_LGR = LGRs, Plate = Group.1, Row = Group.2, Column = Group.3)
    }else{
      aggregatedsd = dplyr::rename(aggregatedsd, sd_LGR = LGR, Plate = Group.1, Row = Group.2, Column = Group.3)
    }
    
    if(input$smooth){
      aggregatedsd <- aggregatedsd [c(1:3,12:13)]
    }else{
      aggregatedsd <- aggregatedsd [c(1:3,12)]
    }
    
    #  Create Paired data.
    t.testx <-xLGR3[,12]-xLGR3[,11]
    
    #  Create new t-test so that it doesn't crash upon the first instance of an error.  
    my_t.test<-function(x){
      A<-try(t.test(x), silent=TRUE)
      if (is(A, "try-error")) return(NA) else return(A$p.value)
    }
    
    #  Use aggregate with new t-test.  
    aggregate.ttest = aggregate(t.testx, by=list(xLGR3$Plate, xLGR3$meanRow, xLGR3$meanColumn),FUN=my_t.test)
    aggregate.ttest = dplyr::rename(aggregate.ttest, p.value = x, Plate = Group.1, Row = Group.2, Column = Group.3)
    ############
    
    aggregated <- inner_join(aggregatedmean, aggregatedsd, by = c("Plate", "Row", "Column"))
    aggregated <- inner_join(aggregated, aggregate.ttest, by = c("Plate", "Row", "Column"))
    
    aggregated1 <- transform(aggregated, as.numeric(as.character(Plate)))
    ### adjust p.value with fdr correction (Benjamini&Hochberg) ###
    q.value = p.adjust(aggregated1$p.value, method = "fdr", n = length(aggregated1$p.value))
    aggregated2 = mutate (aggregated1, q.value)
    ### calculate -ln(p.value) & -ln(q.value) ###
    aggregated2 <- mutate(aggregated2, negLOG.p.value = -log(p.value))
    aggregated2 <- mutate(aggregated2, negLOG.q.value = -log(q.value))
    
    ### calculate Z-score over mean of replicates ###
    mean_Z_score = c()
    mean_Z_score = as.vector(scale(aggregated2$Mean_LGR))
    aggregated2 = mutate(aggregated2, mean_Z_score)
    
    ### add keyfile again to merge mean LGRs with ID information ###
    key = read.csv(input$file2$datapath,
                   header = T,
                   sep = "\t")
    key = dplyr::rename(key, Plate = "Plate")
    key = transform(key, Plate = as.numeric(Plate))
    
    mean_file <- inner_join(aggregated2, key, by = c("Plate", "Row", "Column"))
    
   # if(input$smooth){
    #  mean_file = mean_file[c(1:5,8:21)]
    #}else{
     # mean_file = mean_file[c(1:5,8:19)]
    #}
    ### print selected data (mean data or repliactes data) ###
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
      if(input$disp == "Means") {
        paste("mean_file-", Sys.Date(), ".csv", sep="")
        }else{
          paste("replicates_file-", Sys.Date(), ".csv", sep="")
        }
      },
    
    content = function(file) {
      
      write.csv(getData(), file)
    }
  )
  
  ###########################################################################################
  ### Script for combining to independent controls (second tab on ScreenGarden, optional) ###
  ###########################################################################################

  getData2 <- reactive({
    
    ### input$file1 will be NULL initially. After the user selects
    ### and uploads a file, head of that data file by default,
    ### or all rows if selected, will be shown.
    
    req(input$file3)
    req(input$file4)
    
    ctr1 <- read.csv(input$file3$datapath,
                     header = T)
    ctr2 <- read.csv(input$file4$datapath,
                     header = T)
    
    Merge.file <- inner_join(ctr1, ctr2, by = c("Plate", "Row", "Column", "ID"), suffix = c(".1",".2"))
    Mean_LGR = c() # calculates the mean LGR of both comparisons
    Mean_LGR = rowMeans(Merge.file[c("Mean_LGR.1", "Mean_LGR.2")], na.rm=TRUE)
    Mean_Z = c() # calculates the mean Z score of both comparisons
    Mean_Z = rowMeans(Merge.file[c("mean_Z_score.1", "mean_Z_score.2")], na.rm=TRUE)
    max_p = c() # calculates the max pvalues of both comparisons
    max_p =  pmax(Merge.file$p.value.1, Merge.file$p.value.2)
    max_q = c() # calculates the max qvalues of both comparisons
    max_q = pmax(Merge.file$q.value.1, Merge.file$q.value.2)
    max_negLOG.p = c() # calculates the max neglogP score of both comparisons
    max_negLOG.p = -log(pmax(Merge.file$p.value.1, Merge.file$p.value.2))
    max_negLOG.q = c() # calculates the max neglogQ score of both comparisons
    max_negLOG.q = -log(pmax(Merge.file$q.value.1, Merge.file$q.value.2))
    
    Merge.File = mutate( Merge.file, Mean_LGR, Mean_Z, max_p, max_q, max_negLOG.p, max_negLOG.q)
    
    return(Merge.File)
  })
  output$contents2 <- renderTable(
    
    getData2()
    
  ) 
  ### Downloadable csv of selected dataset ----
  output$downloadData2 <- downloadHandler(
    filename = function() { 
      paste("merge_file-", Sys.Date(), ".csv", sep="")
    },
    
    content = function(file) {
      
      write.csv(getData2(), file)
    }
  )
  
  ###################################################################################################
  ### Script for plotting columns of mean file or merge file(third tab on ScreenGarden, optional) ###
  ###################################################################################################
  
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
    p  = data()
    plot1 <- ggplot(p,
                    aes_string(
                      x =  input$xcol,
                      y = input$ycol
                    ))+ geom_point()+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y ,  expand = TRUE) + geom_hline(yintercept = c(0))
    
    print(plot1)
  })
  
  output$plot1 <- renderPlot({
    plot1()
  })
  
  output$down <- downloadHandler(
    filename = function(){paste(input$data, '.pdf', sep = '')},
    content = function(file){
      ggsave(file, plot1())
    },
    contentType = "application/pdf"
  )
  
  ### plot histogram of data distribution ###
  output$distPlot <- renderPlot({
    B =  data()
    x    <- B$Mean_LGR
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    hist(x, breaks = bins, col = "#75AADB", border = "white",
         xlab = "Mean LGR", 
         ylab = "Number of strains",
         main = "Distribution of LGRs")
  })
  
  #############################################################################################################
  ### Script for calculating q-values from mixture modelk fitting and plotting data as bimodal distribution ###
  ########################### Script adapted from Howell et al., 2019, G3 #####################################
  #############################################################################################################
  
  MMdata <- reactive({
    
    req(input$file6)
    
    MM <- read.csv(input$file6$datapath,
                   header = T)
    
    select = dplyr::select
    MM2 = MM %>% select(ID, Mean_LGR)
    
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
      paste("mixture_model-", Sys.Date(), ".csv", sep="")
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
  
})
# Create Shiny object
#shinyApp(ui = ui, server = server)