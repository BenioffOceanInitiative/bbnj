#' Get Half-Degree Cell Empty Grid
#'
#' @return
#' @export
#'
#' @examples
#' @import raster
get_fishing_empty_grid <- function(){

  # usethis::use_package("raster")
  # devtools::loaded_packages()
  #devtools::load_all()

  # raster specifications for 0.5 degree global raster
  raster(
    xmn = -180, xmx = 180, ymn = -90, ymx = 90,
    resolution=0.5, crs=leaflet:::epsg4326)
}
