# shiny-app C.O.M.I.C.S
Calling Outlier loci from Multi-dimensional data using Invariant Coordinate Selection
A pipeline for identifying outlier loci from multiple selection scans.
# Using COMICS
##Installing Dependencies

COMICS relies on several dependencies:
    install.packages(c("shiny", "mvtnorm", "ICS", "moments", "ICSOutlier", "ggplot2", "reshape", "shinythemes"), dependencies = TRUE)

Users may also need to install additional tools for a more efficient ues of COMICS. See:  https://www.rstudio.com/products/rpackages/devtools/

##Installing COMICS

There are two R-scripts that make up COMICS. The first script contains code that designs the user interface for COMICS, otherwise known as the fluid page. -> COMICS.design.R. The second script includes the input of the data and invokes the ICS and figure generation.

For COMICS istallation after dependencies are downloaded
    install.packages("devtools", dependencies = TRUE)
    devtools::install_github("JTNelsonWSU/shiny-app", build_vignettes = TRUE)
    libraray(shiny-app)
    shiny-app()
   
##Input Files

COMICS incorperates the use of two different input files, both of which are necessary for the complete use of the package. See the example input files for more information about the format of the two files.

input file 1 -> Scan_data.txt -> the main inout file that contains the test statistics for each selection scan for each position on a given chromosome/genome.

input file 2 -> Genome_configuration.txt -> a data table that contains the number of chromosomes and their putativley lengths.

##Additional Questions

If there are any issues that arise, please open a new issue.


