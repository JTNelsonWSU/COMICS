#' Launches COMICS

#' @export launchApp

#' @return shiny application object

#' @importFrom shiny runApp
# @import mvtnorm
# @import ICS
# @import moments
# @import ICSOutlier
# @import ggplot2
# @import reshape2



#'


# wrapper for shiny::shinyApp()
launchApp <- function() {
  shinyApp(ui = ui, server = server)
}
