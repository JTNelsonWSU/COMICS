

####### This script is responsible for the C.O.M.I.C.S app design and all input/output types, including location on the interface.

ui <- fluidPage(      # Tells shiny that the webpage is not going to be static/fixed
  theme = shinytheme("superhero"),
  titlePanel("Thanks for using C.O.M.I.C.S"),  # Adds a main title to the webpage
  sidebarLayout(   # Edits the side bar of the webpage
    sidebarPanel(   # Edits the panel of the sidebar 
       fileInput(inputId = "file1", label = "Data file",   # Adds an inputfile. Name of the inputID is "file1". This file will contain the information to run the ICS.
                accept = c("text/csv", "text/tab-seperated-vlaues,text/plain", ".csv"), # Does not allow more than text or csv files.
                width = NULL, buttonLabel = "Browse",   # Sizes of the box and the label inside it
                placeholder = "No file selected"), # Default placeholder for the input box.
       fileInput(inputId = "configuration", label = "Genome Configuration File", #Adds an inputfile. Name of the inputID is "configuration". Used to construct ICS plot and histogram figures per chromosome.
                 accept = c("text/csv", "text/tab-seperated-values,text/plain", ".csv"), # Does not allow more than text or csv files
                 width = NULL, buttonLabel = "Browse", # Sizes of the box and the label inside it
                 placeholder = "No file selected"), # Default placeholder for the input box.
       numericInput("Cutoff", label = "Statistical Cutoff", value = 5, min = 0, max = 10), # Inputs a numerical option for the cutoff for outlier points. This will be used to determine ICS outliers (Values represent the top X%).
      selectInput("dataset", "Choose a dataset",  # For downloading the final ICS dataset, including putative outliers.
                  choices = c("ICS Distance", "ICS genome scan")),  # Dataset choices for downloading.   
       numericInput("TestOfInterest", label = "Test of interest", value = 1, min = 0, max = 100), # Inputs a numerical option that makes a plot for the single selection test selected. 
       # This input uses data from the original input file "file1" and does not include ICS results.
       numericInput("Chromosomes", label = "Chromosome of interest", value = 1, min = 0, max = "configuration"), # Inputs a numerical value that makes a histogram for any given chromosome selected.
       # This input uses data from the final ICS calculation.

      downloadButton("downloadData", "Download ICS output"), # Makes the download button.
      downloadButton("downloadPlot", "Download ICS figure") # Makes the download button.
      ),
    
    
    mainPanel(     # Calls the main panel.
      h4("Summary.ICS"),
      verbatimTextOutput("summary.ICS"),
      plotOutput("ICS.hist"), # The histogram for the log_ICS distaces. This figure includes all estimates genome-wide.   
      plotOutput("ICS.chromosome.hist"), # Histogram for the log_ICS distances for a selected chromosome.
      plotOutput("GenomeScan"), # The final plot of ICS data that maps the log_ICS distances per chromosomes.
      h3("Summary.SingleTest"),
      verbatimTextOutput("summary.SingleTest"),
      plotOutput("GenomeTest"), # Calls a plot that uses non-ICS data from "file1" to plot a genome scan of a given selection test ("Test of interest"). 
      plotOutput("GenomeMelt"), # A plot that has the melted data from "file" to compare all selection tests prior to running ICS.
      verbatimTextOutput("nText")
      )
  )
)

