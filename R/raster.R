#' Get empty global grid
#'
#' @param res resolution, in decimal degrees (default 0.1)
#' @param crs coordinate resolutions system, ie projection (default geographic coordinate system 4326, ie WGS84)
#'
#' @return global raster with extent -180 to 180 longitude, -90 to 90 latitude and specified resolution
#' @export
#'
#' @examples
#'
#' # get half-degree global grid
#' get_grid(res=0.5)
#' @import raster leaflet
get_grid <- function(res = 0.1, crs=leaflet:::epsg4326, val=NA){
  # raster specifications for 0.5 degree global raster

  raster::raster(
    xmn = -180, xmx = 180, ymn = -90, ymx = 90,
    resolution=res, crs=leaflet:::epsg4326, vals=val)
}

rescale_grid <- function(r, res_new=0.1, fxn="ngb"){
  # TODO: resample(x, y, method="bilinear")

  res_old <- raster::res(r)[1]
  if (res_old == res_new) r
  if (res_old %% res_new != 0) stop(glue("res_old {res_old} not evenly divisible by res_new {res_new}"))
  if (res_new > res_old){
    r <- raster::aggregate(r, res_new/res_old)
  } else {
    r <- raster::disaggregate(r, res_old/res_new)
  }
}
