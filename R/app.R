#' Run the Shiny app map viewer for BBNJ inputs and
#'
#' @return Launches the app locally in your web browser. If \code{inst/app} directory
#'   exists in your current path (as when working from the R project of the bbnj
#'   repository) it will launch that version.
#' @export
run_app <- function(){
  library(shiny)

  if (dir.exists("inst/app/")){
    dir_app <- "inst/app/"
  } else {
    dir_app <- system.file(package="bbnj", "app")
  }
  shiny::runApp(dir_app)
}
