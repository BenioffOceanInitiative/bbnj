#' Helper function to run prioritizr::solve() and log output and resulting raster
#'
#' @param p \code{\link[prioritizr]{problem}} to feed into \code{\link[prioritizr]{solve}}
#' @param log path to output log from solving problem
#'
#' @return output from \code{\link[prioritizr]{solve}}
#' @export
#' @import prioritizr raster
#'
#' @examples
#solve_log <- function(p, log, tif){
solve_log <- function(p, pfx= deparse(substitute(p))){
  library(prioritizr)
  library(raster)
  library(glue)

  log <- glue("{pfx}_log.txt")
  tif <- glue("{pfx}_sol.tif")

  sink(log, append=T)
  s <- prioritizr::solve(p)
  sink()

  if ("RasterLayer" %in% class(s)){
    raster::writeRaster(s, tif, overwrite=T)
  }

  s
}
