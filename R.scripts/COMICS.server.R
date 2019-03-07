
#This script is the core script for shiny that first calls "COMICS.fluid.R" which creates a fluid user interface, then performs the server function. 

library(shiny)   # Loads shiny package.
library(mvtnorm) # Loads the Multivariate Normal and t Distribution package.
library(ICS)     # Loads the ICS package.
library(moments) # Loads the momment (Distribution moments) package.
library(ICSOutlier) # Loads a sepeerate ICS package that calls outliers from the ICS.
library(ggplot2) # Loads ggplot 2 for figure design.
library(reshape2) # Loads the reshape package for melting data.
library(shinythemes)

source("shiny.design.v5.JTN.R")   # Calls the "shiny.design.R" script to create/design the user interface 
# for more information, look at "shiny.design.R" annotation.

server <- function(input, output) {  # Calls the shiny server. From hearafter, all R code focues on "actions and functions" and do not pertain to the design of the app. 

 
  observe({    # Creates a "reactive" environment that can read reactive values and expressions. Importantly, the "observe" expression will automatically re-execute when dependencies are changed. 
      
    file1 =  input$file1   # Defines file1 as input$file1, which is what we told the UI above in the "shiny.design.R" script. File format (headers) = chrom, midpoint, test1, test2, test3, test4, testn. 
    if (is.null(file1)) { # Checks to see if the file has data in it
      return(NULL) # If are absent/in the wrong format, the program will terminate with an error.
    }  # Close if statement.
    data.1 = read.table(paste(file1$datapath), header = TRUE)  # Reads in the same file as above, however, we save it as a different object, "data.1"
    
    fig.test.title <- input$TestTitle # saves the test title as a new variable. This is a continuation from the "shiny.design.R" script.
  
    TestX <- input$TestOfInterest  # Creats a new variable from the input "TestOfInterest" in the UI.
    TestXF <- data.1[,2 + TestX]  # Creats a new variable that tells the server which selection test to graph.
    
    output$GenomeTest <- renderPlot({   # The following code uses ggplot2 to construct/render a scatter plot for a individual selection test (pre-ICS data). This plot is called "GenomeTest" in the UI.
      ggplot(data = data.1, aes(x = data.1$Midpoint, y = TestXF)) + geom_point(size=1.0, aes(colour = factor(Chrom))) + 
        theme(axis.text.x= element_text(size = rel(1.1), hjust = 1, vjust = 1)) + ylab(expression(test_statistic)) + 
        xlab(expression(Position(bp))) + guides(colour=FALSE) + ggtitle("Test of interest") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        theme(plot.title = element_text(hjust =0.5)) +
        theme(axis.title.x = element_text(size = 14)) +
        theme(axis.text.y = element_text(size = rel(1.3), color = "black")) +
        theme(axis.title.y = element_text(hjust = 0.5, size = 14)) +
        theme(plot.title = element_text(face = "bold", size = 20))
    })
    
  data.1.name <- colnames(data.1)  # Saves the columns names of data.1 into data.1.name. This is used to prep data.1 for metling with the reshape2 package.
  data.1.name <- data.1.name[3:length(data.1)] # Removes the first two column names in a file so that data.1.names will be a vector with length equal to the number of selection tests performed.
    
  data.1.melt <- melt(data.1, id.vars = c("Midpoint", "Chrom"), measure.vars = data.1.name) # data.1 is melted by test name and identified by the "Midpoint" and "Chrom" columns. 
    
    output$GenomeMelt <- renderPlot({   # The following code uses ggplot2 to construct/render a melted scatter plot for all selection test (pre-ICS data) individually. This plot is called "GenomeMelt" in the UI. Also not the "free_scale" flag accounting for different test statistics. 
      ggplot(data = data.1.melt, aes(x = Midpoint, y = value )) + geom_point(size=1.0, aes(colour = factor(Chrom))) + 
      theme(axis.text.x= element_text(size = rel(1.1), hjust = 1, vjust = 1)) + ylab(expression(test_statistic)) + 
      xlab(expression(Position(bp))) + guides(colour=FALSE) + ggtitle("Individual Selection Scans") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
      theme(plot.title = element_text(hjust =0.5)) +
      theme(axis.title.x = element_text(size = 14)) +
      theme(axis.text.y = element_text(size = rel(1.3), color = "black")) +
      theme(axis.title.y = element_text(hjust = 0.5, size = 14)) +
      theme(plot.title = element_text(face = "bold", size = 20)) + 
      facet_grid(variable ~ ., scales = "free")
    })
  
  output$summary.SingleTest <- renderPrint({
      summary(data.1[,-1:-2])
    })
    
          
}) # Closes the first "observe" environment.
  
  
  observe({     # Creates a "reactive" environment that can read reactive values and expressions. Importantly, the "observe" expression will automatically re-execute when dependencies are changed. 
    
      
   file1 =  input$file1   # Defines file1 as input$file1, which is what we told the UI above. Note, this is the same input file as data.1, used in the previous observe function.
    if (is.null(file1)) { # Checks to see if the file has data in it.
      return(NULL) # If are absent/in the wrong format, the program will terminate with an error.
    }
    dataX = read.table(paste(file1$datapath), header = TRUE)  # Reads in the same file as above, however, we save it as a different object, "dataX".
    X = dataX  # Makes a new object from "dataX", called "X".
    

    configuration = input$configuration   # Defines "configuration" as input$configuration, which is what we defined in the UI above.
    if (is.null(configuration)) {  # Checks to see if the file has data in it
      return(NULL) # If are absent/in the wrong format, the program will terminate with an error.
    }
    dataY = read.table(paste(configuration$datapath), header = TRUE)  # Reads in the same file as above, however, we save it as a different object, "dataY"
    chr.length = dataY  # Makes a new object from "dataY", called "chr.length". This file will be used later in the pipeline for figure porduction. This is file is the second of two input files, called the Genome Configuration File (GCF).
    

########Start ICS############
    X[,1:2] <- NULL  # Removes the first 2 columns in the file (prepares file for ICS).
    Z <- ics2(X)  # Performs ICS on object "X" and saves it as a new object, "Z".
    Z_comp <- comp.norm.test(Z) # Identifies invariant coordinates that are non normal using univariate normality tests and saves it as "Z_comp".
    Z_dist <- ics.distances(Z) # Computes the squared ICS distances and saves it as "Z_dist".
    Z_frame <- data.frame(Z_dist) # Makes the ICS distance a data frame, "Z_frame".
    #icsOutlierJB <- ics.outlier(Z, test = "jarque", level.dist = 0.05, level.test = 0.05, mDist = 100) # Identifies outliers based on the two scatter matrices using the Jarque Test. Creats a new variable called "icsOutliersJB".
    #icsOutlierAG <- ics.outlier(Z, test = "anscombe", level.dist = 0.05, level.test = 0.05, mDist = 100) # Identifies outliers based on the two scatter matrices using the Anscombe Test. Creates a new variable called "icsOutliersAG".
    Z2 <- cbind(dataX$Midpoint, Z_frame$Z_dist) # Combines the Midpoints, outliers, and ics distance into one object, called "Z2".
    Z2 <- data.frame(Z2)  # Makes Z2 a data frame.
    Z2 <- cbind(dataX$Chrom, Z2) # Adds a chromosome column to the Z2 data frame.
    colnames(Z2) <- c("Chr","Midpoint","ics.distance") # Changes column names in Z2.
    Log.ICS <- log10(Z2$ics.distance) # Logs the ICS distance values.
    Z2 <- cbind(Z2,Log.ICS) # Combines "Z2" object with the log ICS distances.
    
    colnames(Z2) <- c("Chr","Midpoint","ICS.distance","Log.ICS") # Changes column names in Z2

    n<-input$Cutoff # makes the statistical cutoff a variable. "input$Cutoff" was designed in the UI.
    cutoff.applied <- quantile(Z2$Log.ICS, probs = 1 - n/100)

    
########End ICS#############   
#######Start figures for ICS########
    
  output$value <- renderText({ genome = input$Chromosomes})  # Defines the number of chromosomes based off the GCF.
  genome <- input$Chromosomes # Defines the number of chromosomes based off the GCF.
  
  output$ggplotTitle <- renderText({fig.title = input$Title}) # Makes title for the ICS genome scan.
  fig.title <- input$Title # saves title as a new variable
  
  if(Z2[,1] == 1){   # The following if statement takes the midpoint value of each window and adds it to the corresponding chromosome length found in the GCF.
    Z2[,5] <- Z2[,2]
  } # End if statement.
  
  for(j in 2:genome){      # The following loop continues adds to the previous if statment and aligns the positions in order for each chromosome. Note that this is a nested loop and an if statement. One loop per chromosomes and another loop for the posistions in each chromosomes.
    for(i in 1:dim(Z2)[1]){
      if (Z2[i,1] == j) {
        Z2[i,5] <- Z2[i,2] + sum(chr.length[1:(j - 1),2])}
    }
  } # End nested loops/if statments
  
  
  colnames(Z2) <- c("Chr","Midpoint","ICS.Distance","Log.ICS","Position") # Changes column names in Z2.
  
Z2$Outlier <- 0
  
Z2$Outlier[ Z2$Log.ICS >= cutoff.applied ] <- 1

Z2$Color <- "grey20"


for(i in seq(1,dim(dataY)[1],2)){
  Z2$Color[Z2$Chr == i] <- "grey58"
}    

Z2$Color[ Z2$Log.ICS >= cutoff.applied ] <- "darkred"
Z2$Color <- factor(Z2$Color)

  
  
  output$ICS.hist <- renderPlot({   # Renders a genome-wide histogram of the log_ICS distances and applies a cutoff. 
    #ggplot(data = Z2, aes(x=Log.ICS)) + geom_histogram(binwidth = 0.05) + geom_vline(mapping = NULL, data = NULL, xintercept = log(attributes(Outliers)$ics.dist.cutoff), color = "red")
    ggplot(data = Z2, aes(x=Log.ICS)) + geom_histogram(binwidth = 0.05) + geom_vline(mapping = NULL, data = NULL, xintercept = cutoff.applied, color = "red")
  })
  
 
  Z2.chr <- matrix(data = NA,nrow = dim(Z2)[1], ncol = 1) # Makes an empty matrix called "Z2.chr". This is used to render the "input$Chromosomes" (single chromosome histograms) defined in the UI.
  Z2.chr <- data.frame(Z2.chr)  # Makes Z2.chr a data frame.
  
   
  chr.interest <- input$Chromosomes  # Creats a new variable, "chr.interest". Here, the "input$Chromosome" is saved as "chr.interest".
  for(i in 1:dim(Z2)[1]){   # The following loop and if satement goes through each line of the Z2 data frame and looks for a chromosome match that equals "chr.interest". If there is a match, the corresponding log_ICS distance is save in "Z2.chr".
    if(Z2[i,1] == chr.interest){
    Z2.chr[i,1] <- Z2[i,4]
    }  # If there is no match, R will do nothing.
  }  # Close loop and if statement.
  
   output$ICS.chromosome.hist <- renderPlot({   # Renders the "ICS.chromosome.hist" historgram for any chromosome chosen by the user with the desired cutoff. This figured is produced from the "Z2.chr" data frame and will contain many "NA" warnings. 
    ggplot(data = Z2.chr, aes(x= Z2.chr)) + geom_histogram(binwidth = 0.8) + 
      geom_vline(mapping = NULL, data = NULL, xintercept = cutoff.applied, color = "red")
  }) 
  
   
    
    
  
  output$GenomeScan <- renderPlot({   # Renders the genome scan for ICS, called "GenomeScan". Note the different variables used in this figure ("Z2", "fig.title", and "Outliers$ics.dist.cutoff"). This figure will be color coded by chromosome ("colour = factor(Chr)").
    ggplot(data = Z2, aes(x = Midpoint, y = Log.ICS, color = Color)) + geom_point(size=1.0) + guides(colour = FALSE) + scale_color_manual(values = c("darkred", "grey20", "grey58")) + 
      theme(axis.text.x= element_text(size = rel(1.1), hjust = 1, vjust = 1)) + ylab(expression(Log_ICS_Distance)) + 
      xlab(expression(Position(bp))) + guides(colour=FALSE) + ggtitle("ICS genome scan") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
      geom_hline(yintercept=cutoff.applied) + theme(plot.title = element_text(hjust =0.5)) +
      theme(axis.title.x = element_text(size = 14)) +
      theme(axis.text.y = element_text(size = rel(1.3), color = "black")) +
      theme(axis.title.y = element_text(hjust = 0.5, size = 14)) +
      theme(plot.title = element_text(face = "bold", size = 20))
  })
  
  
  output$summary.ICS <- renderPrint({
    summary(Z2$Log.ICS)
  })
  
  ############End ICS figures##############
  output$downloadData <- downloadHandler(  # Renders the download for a dataset and the type of file format.
    filename = function() {
      paste(input$dataset, ".csv")
    },
    content = function(file) {    #Uused for downloading the Z2 data file. Descibes what the output format will look like.
      write.csv(datasetInput(), file, row.names = FALSE)
    }) 
  
  output$downloadPlot <- downloadHandler(
    filename = "test.png",
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
      ggsave(file, plot = output$GenomeScan(), device = device)
   }
  )


  
  
  datasetInput <- reactive({ # Uses switch function to switch names with datasets. Here, we switch the "input$dataset" to "Z2", so that the download option will appear as "ICS.Distance", where "ICS.Distance" equals the "Z2" data frame.
    switch(input$dataset,
           "ICS Distance" = Z2)
    
  })
  
  })  # Close the "observe" environment.


}  # Close the shiny server  

shinyApp(ui = ui, server = server) # Calls the shiny app and puts the UI and the server together.



