\name{server}
\alias{server}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
server(input, output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{input} here~~
}
  \item{output}{
%%     ~~Describe \code{output} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (input, output) 
{
    observe({
        file1 = input$file1
        if (is.null(file1)) {
            return(NULL)
        }
        data.1 = read.table(paste(file1$datapath), header = TRUE)
        fig.test.title <- input$TestTitle
        TestX <- input$TestOfInterest
        TestXF <- data.1[, 2 + TestX]
        output$GenomeTest <- renderPlot({
            ggplot(data = data.1, aes(x = data.1$Midpoint, y = TestXF)) + 
                geom_point(size = 1, aes(colour = factor(Chrom))) + 
                theme(axis.text.x = element_text(size = rel(1.1), 
                  hjust = 1, vjust = 1)) + ylab(expression(test_statistic)) + 
                xlab(expression(Position(bp))) + guides(colour = FALSE) + 
                ggtitle("Test of interest") + theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), panel.background = element_blank()) + 
                theme(plot.title = element_text(hjust = 0.5)) + 
                theme(axis.title.x = element_text(size = 14)) + 
                theme(axis.text.y = element_text(size = rel(1.3), 
                  color = "black")) + theme(axis.title.y = element_text(hjust = 0.5, 
                size = 14)) + theme(plot.title = element_text(face = "bold", 
                size = 20))
        })
        data.1.name <- colnames(data.1)
        data.1.name <- data.1.name[3:length(data.1)]
        data.1.melt <- melt(data.1, id.vars = c("Midpoint", "Chrom"), 
            measure.vars = data.1.name)
        output$GenomeMelt <- renderPlot({
            ggplot(data = data.1.melt, aes(x = Midpoint, y = value)) + 
                geom_point(size = 1, aes(colour = factor(Chrom))) + 
                theme(axis.text.x = element_text(size = rel(1.1), 
                  hjust = 1, vjust = 1)) + ylab(expression(test_statistic)) + 
                xlab(expression(Position(bp))) + guides(colour = FALSE) + 
                ggtitle("Individual Selection Scans") + theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), panel.background = element_blank()) + 
                theme(plot.title = element_text(hjust = 0.5)) + 
                theme(axis.title.x = element_text(size = 14)) + 
                theme(axis.text.y = element_text(size = rel(1.3), 
                  color = "black")) + theme(axis.title.y = element_text(hjust = 0.5, 
                size = 14)) + theme(plot.title = element_text(face = "bold", 
                size = 20)) + facet_grid(variable ~ ., scales = "free")
        })
        output$summary.SingleTest <- renderPrint({
            summary(data.1[, -1:-2])
        })
    })
    observe({
        file1 = input$file1
        if (is.null(file1)) {
            return(NULL)
        }
        dataX = read.table(paste(file1$datapath), header = TRUE)
        X = dataX
        configuration = input$configuration
        if (is.null(configuration)) {
            return(NULL)
        }
        dataY = read.table(paste(configuration$datapath), header = TRUE)
        chr.length = dataY
        X[, 1:2] <- NULL
        Z <- ics2(X)
        Z_comp <- comp.norm.test(Z)
        Z_dist <- ics.distances(Z)
        Z_frame <- data.frame(Z_dist)
        Z2 <- cbind(dataX$Midpoint, Z_frame$Z_dist)
        Z2 <- data.frame(Z2)
        Z2 <- cbind(dataX$Chrom, Z2)
        colnames(Z2) <- c("Chr", "Midpoint", "ics.distance")
        Log.ICS <- log10(Z2$ics.distance)
        Z2 <- cbind(Z2, Log.ICS)
        colnames(Z2) <- c("Chr", "Midpoint", "ICS.distance", 
            "Log.ICS")
        n <- input$Cutoff
        cutoff.applied <- quantile(Z2$Log.ICS, probs = 1 - n/100)
        output$value <- renderText({
            genome = input$Chromosomes
        })
        genome <- input$Chromosomes
        output$ggplotTitle <- renderText({
            fig.title = input$Title
        })
        fig.title <- input$Title
        if (Z2[, 1] == 1) {
            Z2[, 5] <- Z2[, 2]
        }
        for (j in 2:genome) {
            for (i in 1:dim(Z2)[1]) {
                if (Z2[i, 1] == j) {
                  Z2[i, 5] <- Z2[i, 2] + sum(chr.length[1:(j - 
                    1), 2])
                }
            }
        }
        colnames(Z2) <- c("Chr", "Midpoint", "ICS.Distance", 
            "Log.ICS", "Position")
        Z2$Outlier <- 0
        Z2$Outlier[Z2$Log.ICS >= cutoff.applied] <- 1
        Z2$Color <- "grey20"
        for (i in seq(1, dim(dataY)[1], 2)) {
            Z2$Color[Z2$Chr == i] <- "grey58"
        }
        Z2$Color[Z2$Log.ICS >= cutoff.applied] <- "darkred"
        Z2$Color <- factor(Z2$Color)
        output$ICS.hist <- renderPlot({
            ggplot(data = Z2, aes(x = Log.ICS)) + geom_histogram(binwidth = 0.05) + 
                geom_vline(mapping = NULL, data = NULL, xintercept = cutoff.applied, 
                  color = "red")
        })
        Z2.chr <- matrix(data = NA, nrow = dim(Z2)[1], ncol = 1)
        Z2.chr <- data.frame(Z2.chr)
        chr.interest <- input$Chromosomes
        for (i in 1:dim(Z2)[1]) {
            if (Z2[i, 1] == chr.interest) {
                Z2.chr[i, 1] <- Z2[i, 4]
            }
        }
        output$ICS.chromosome.hist <- renderPlot({
            ggplot(data = Z2.chr, aes(x = Z2.chr)) + geom_histogram(binwidth = 0.8) + 
                geom_vline(mapping = NULL, data = NULL, xintercept = cutoff.applied, 
                  color = "red")
        })
        output$GenomeScan <- renderPlot({
            ggplot(data = Z2, aes(x = Midpoint, y = Log.ICS, 
                color = Color)) + geom_point(size = 1) + guides(colour = FALSE) + 
                scale_color_manual(values = c("darkred", "grey20", 
                  "grey58")) + theme(axis.text.x = element_text(size = rel(1.1), 
                hjust = 1, vjust = 1)) + ylab(expression(Log_ICS_Distance)) + 
                xlab(expression(Position(bp))) + guides(colour = FALSE) + 
                ggtitle("ICS genome scan") + theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), panel.background = element_blank()) + 
                geom_hline(yintercept = cutoff.applied) + theme(plot.title = element_text(hjust = 0.5)) + 
                theme(axis.title.x = element_text(size = 14)) + 
                theme(axis.text.y = element_text(size = rel(1.3), 
                  color = "black")) + theme(axis.title.y = element_text(hjust = 0.5, 
                size = 14)) + theme(plot.title = element_text(face = "bold", 
                size = 20))
        })
        output$summary.ICS <- renderPrint({
            summary(Z2$Log.ICS)
        })
        output$downloadData <- downloadHandler(filename = function() {
            paste(input$dataset, ".csv")
        }, content = function(file) {
            write.csv(datasetInput(), file, row.names = FALSE)
        })
        output$downloadPlot <- downloadHandler(filename = function() {
            paste(input$dataplot, ".png")
        }, content = function(plots) {
            ggsave("ICS.plot.png")
        })
        datasetInput <- reactive({
            switch(input$dataset, `ICS Distance` = Z2)
        })
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
