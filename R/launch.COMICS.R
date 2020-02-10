#' Launches COMICS

#' @export launchApp

#' @return shiny application object

#' @import shiny
#' @import reshape2

#'


# wrapper for shiny::shinyApp()
launchApp <- function() {
  shinyApp(ui = ui, server = server)
}
