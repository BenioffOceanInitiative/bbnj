#' Get empty global grid
#'
#' @param res resolution, in decimal degrees (default 0.5)
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
get_grid <- function(res = 0.5, crs=leaflet:::epsg4326, val=NA){
  # raster specifications for 0.5 degree global raster

  raster::raster(
    xmn = -180, xmx = 180, ymn = -90, ymx = 90,
    resolution=res, crs=leaflet:::epsg4326, vals=val)
}

#' Quick leaflet map of raster
#'
#' @param r raster
#' @param r_title title of raster layer
#' @param na0 convert 0s to NA (default = FALSE)
#' @param method interpolation method for web mercator, "ngb" (nearest neighbor) or "bilinear" (default "ngb")
#' @param boundary_shp path to high seas boundary shapefile
#'
#' @return
#' @export
#'
#' @examples
qmap_r <- function(
  r, r_title, na0=T, method="ngb",
  boundary_shp="~/Gdrive Ecoquants/projects/bbnj/data/derived/boundary/high_seas_s05.shp"){
  # b2_shp="~/Google Drive/projects/bbnj/data/derived/high-seas_boundary/high_seas_final_v2.shp"

  # debug
  # boundary_shp="~/Google Drive/projects/bbnj/data/derived/high-seas_boundary/high_seas_final.shp"
  # r = r_sol; r_title = "Solution"
  # r = raster("~/Google Drive/projects/bbnj/data/derived/biodiversity/spp_bivalves_be_0.5dd.tif")
  # r_title = "bivalves"
  # na0=T; method="ngb"

  boundary <- read_sf(boundary_shp) %>%
    lwgeom::st_make_valid() # plot(boundary)

  if (na0) r[r==0] <- NA
  #plot(r)

  pal <- colorNumeric("Spectral", values(r), reverse=T, na.color = "transparent")

  is_1val <- ifelse(length(unique(r)) == 1, T, F)
  if (is_1val) pal <- colorNumeric("Reds", values(r), reverse=T, na.color = "transparent")

  r

  suppressWarnings({ # b/c addRasterImage() -> rgdal::rawTransform ... point(s) not finite
  m <- leaflet(
    options = leafletOptions(
      worldCopyJump = T)) %>%
    addProviderTiles(
      providers$Esri.OceanBasemap,
      options = providerTileOptions(opacity=0.8),
      group = "Esri.OceanBasemap") %>%
    addProviderTiles(
      providers$Stamen.TonerLite,
      options = providerTileOptions(opacity=0.8),
      group = "Stamen.TonerLite") %>%
    addRasterImage(
      r, colors = pal, opacity=0.8, group=r_title, method=method) %>%
    addPolygons(
      data=boundary, fill=F, weight=2, group="high seas") %>%
    addLayersControl(
      baseGroups = c("Esri.OceanBasemap", "Stamen.TonerLite"),
      overlayGroups = c(r_title, "high seas"))
  })

  if (is_1val){
    m <- m %>%
      addLegend(
        group = r_title, colors = pal(unique(r)),
        labels = unique(r),
        position = "bottomright",
        title = r_title)
  } else {
    m <- m %>%
      addLegend(
        group = r_title, pal = pal,
        position = "bottomright",
        values = values(r),
        title = r_title)
  }
  m
}


#' Gap-fill rasters with value based on function and mask
#'
#' @param r raster input
#' @param fxn function, default="min"
#' @param r_mask raster mask, default=r_pu_id
#'
#' @return raster after filling NAs with result of function and mask
#' @export
#'
#' @examples
gap_fill_raster <- function(r, fxn="min", r_mask=r_pu_id){
  r[is.na(r)] <- cellStats(r, fxn)
  mask(r, r_pu_id)
}

#' Rescale raster
#'
#' @param r input raster
#' @param rescale T/F to rescale values 0 to 1, default=T
#' @param log T/F to log transform, defaul=F
#' @param inverse T/F to invert by subtracting from 1
#' @param multiply_area T/F to multiply raster by 1
#'
#' @return raster
#' @export
#'
#' @examples
rescale_raster <- function(r, rescale = T, log = F, inverse = F,  multiply_area=F){
  if (multiply_area){
    r <- r * area(r)
  }

  v <- raster::values(r)

  if (log){
    if (rescale){
      v <- scales::rescale(v)
    }
    v[v == 0] <- 0.0001
    v <- log(v)
  }

  if (rescale){
    v <- scales::rescale(v)
  }

  if (inverse){
    v <- 1 - v
  }

  raster::values(r) <- v

  r

}

#' Rescale stack
#'
#' @param s stack
#' @param from_all rescale based on range from all rasters in stack
#'   (default=TRUE) or within each raster (FALSE)
#' @param rescale T/F to rescale values 0 to 1, default=T
#' @param log T/F to log transform, defaul=F
#' @param inverse T/F to invert by subtracting from 1
#' @param multiply_area T/F to multiply raster by 1
#'
#' @return raster
#' @export
#'
#' @examples
rescale_stack <- function(s, from_all = T, rescale = T, log = F, inverse = F,  multiply_area=F){

  s_maxs <- cellStats(s, "max")
  s_mins <- cellStats(s, "min")
  r_a <- area(s)

  for (i in 1:nlayers(s)){

    r <- raster(s, i)

    if (from_all){
      rng <- c(min(s_mins), max(s_maxs))
    } else {
      rng <- c(s_mins[i], s_maxs[i])
    }

    if (multiply_area){
      r <- r * r_a
    }

    v <- raster::values(r)

    if (log){
      if (rescale){
        v <- scales::rescale(v, from=rng)
      }
      v[v == 0] <- 0.0001
      v <- log(v)
    }

    if (rescale){
      if (log){
        v <- scales::rescale(v)
      } else {
        v <- scales::rescale(v, from=rng)
      }
    }

    if (inverse){
      v <- 1 - v
    }

    raster::values(r) <- v

    s[[i]] <- r
  }

  s
}

