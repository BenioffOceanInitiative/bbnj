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

  library(lwgeom)

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
gap_fill_raster <- function(r, fxn="min", r_mask=r_pu_id, debug=F){
  r[is.na(r)] <- cellStats(r, fxn)
  mask(r, r_mask)
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
    if (get_r_projection(s)$prj == "gcs"){
      r_a <- area(s)
    } else {
      # assume in meters
      cell_km2 <- prod(res(s)) / (1000 * 1000)
      r_a <- raster::setValues(s[[1]], cell_km2)
    }
    names(r_a) = "area_km2"

    r <- r * r_a
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

  # s = s_seamounts; from_all = T

  s_maxs <- cellStats(s, "max")
  s_mins <- cellStats(s, "min")

  for (i in 1:nlayers(s)){

    r <- raster(s, i)

    if (from_all){
      rng <- c(min(s_mins), max(s_maxs))
    } else {
      rng <- c(s_mins[i], s_maxs[i])
    }

    if (multiply_area){
      if (get_r_projection(s)$prj == "gcs"){
        r_a <- area(s)
      } else {
        # assume in meters
        cell_km2 <- prod(res(s)) / (1000 * 1000)
        r_a <- raster::setValues(s[[1]], cell_km2)
      }
      names(r_a) = "area_km2"

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

#' Convert raster tif to polygon shapefile in geographic coordinates for mapping
#' prioritizr solution in non-Mercator projection
#'
#' @param tif path to raster tif
#'
#' @return True if successfully writes shapefile to same path as tif (.tif ->
#'   _gcs.shp)
#' @export
#'
#' @examples
tif_to_shp_gcs <- function(tif){
  # devtools::load_all
  # tif = "/Users/bbest/github/bbnj/inst/app/www/scenarios/s00a.bio.30pct.gl.mer_sol.tif"
  library(rmapshaper)
  library(lwgeom)

  shp <- glue("{fs::path_ext_remove(tif)}_gcs.shp")

  r <- raster(tif)
  names(r) <- "value"
  #ply <- rasterToPolygons(r, dissolve=T) %>%
  ply <- rasterToPolygons(r) %>%
    st_as_sf() %>%
    filter(!is.na(value), value > 0) %>%
    #group_by(value) %>%
    #summarize() %>%
    ms_dissolve() %>%
    #st_make_valid() %>%
    st_transform(4326) %>%
    lwgeom::st_make_valid() #%>%
  # ply0 <- ply
  ply <- suppressMessages(suppressWarnings(
    st_buffer(ply, 0)))
  ply <- ply %>%
    st_cast("POLYGON", warn=F, do_split=T) %>%
    mutate(
      area_km2 = st_area(geometry) %>% units::set_units(km^2)) %>%
    arrange(desc(area_km2))
  #plot(ply)

  write_sf(ply, shp)
  shp
}
# for (tif in list.files(here("inst/app/www/scenarios"), "tif$", full.names = T)){
#   tif_to_shp_gcs(tif)
# }

get_global_bb <- function(crs=4326){
  # globe in geographic coordinates
  bb = st_sf(
    tibble(
      one = 1,
      geom = st_sfc(st_polygon(list(rbind(
        c(-180, -90),
        c( 180, -90),
        c( 180,  90),
        c(-180,  90),
        c(-180, -90)))), crs = 4326)))
  if (crs != 4326){
    # project
    #crs=54009
    #crs=3857
    #bb_gcs <- bb
    bb <- bb_gcs
    plot(bb)
    bb <- st_buffer(bb, -0.1) %>%
      st_transform(crs)
    plot(bb)
  }
  bb
}


#' Get dataset object in given projection-resolution
#'
#' @param dataset see datasets
#' @param prjres projection-resolution, per projections_tbl$prjres
#' @param debug defaults to False
#'
#' @return polygon (s), raster (r), or stack (s) per prefix of dataset in given projection-resolution (prjres)
#' @export
get_d_prjres <- function(dataset, prjres="", debug=F){
  library(dplyr)
  library(glue)
  library(stringr)
  # dataset = "p_eez_s05"

  dir_data <- system.file("data", package="bbnj")

  type  <- str_sub(dataset, end=1)
  name  <- str_sub(dataset, start=3)

  P <- projections_tbl %>%
    dplyr::filter(prjres == !!prjres)

  #browser()
  if (type == "p"){
    # polygon
    prj <- ifelse(P$prj == "gcs", "", glue("_{P$prj}"))
    path <- glue("{dir_data}/{name}{prj}.shp")
    if (!file.exists(path)) stop(glue("Missing path: {path}"))
    stopifnot(nrow(P) == 1)
    if (debug) message(glue("path: {path}"))
    p <- read_sf(path)
    st_crs(p) <- P$proj
    return(p)
  }
  if (type == "r"){
    # raster
    path <- glue("{dir_data}/{name}{prjres}.tif")
    if (!file.exists(path)) stop(glue("Missing path: {path}"))
    if (debug) message(glue("path: {path}"))
    r <- raster(path)
    crs(r) <- P$proj
    return(r)
  }
  if (type == "s"){
    # stack
    path <- glue("{dir_data}/{name}{prjres}")
    if (!dir.exists(path)) stop(glue("Missing path: {path}"))
    if (debug) message(glue("path: {path}"))
    s <- stack(list.files(path, ".*\\.tif$", full.names=T))
    crs(s) <- P$proj
    return(s)
  }
  # TODO: "s" for stack is folder w/ tifs inside
  #list.files(file.path(system.file(package="bbnj"),"data"))
  #devtools::load_all()

}

#' Get biogeographic stack, given taxonomic grouping, modeling period and projection-resolution
#'
#' @param dataset see datasets
#' @param prjres projection-resolution, per projections_tbl$prjres
#' @param debug defaults to False
#'
#' @return polygon (s), raster (r), or stack (s) per prefix of dataset in given projection-resolution (prjres)
#' @export
get_gmbi_grpsmdl_prjres <- function(grpsmdl="groups00", prjres="", debug=F){
  library(glue)
  library(stringr)
  # dataset = "p_eez_s05"
  # grpsmdl="groups00"; prjres=""; debug=F; dataset = "s_bio_gmbi"
  # get_gmbi_grpsmdl_prjres("groups05", "mol50km")
  # get_gmbi_grpsmdl_prjres("groups06_2100", "mol50km")

  dir_data <- system.file("data", package="bbnj")

  P <- projections_tbl %>%
    filter(prjres == !!prjres)

  # stack
  path <- glue("{dir_data}/bio_gmbi_{grpsmdl}{prjres}")
  if (!dir.exists(path)) stop(glue("Missing path: {path}"))
  stopifnot(nrow(P) == 1)
  if (debug) message(glue("path: {path}"))
  s <- stack(list.files(path, ".*\\.tif$", full.names=T))
  crs(s) <- P$proj
  return(s)
}

#' Project raster to projection-raster
#'
#' @param prjres projection-resolution, per projections_tbl$prjres
#' @param r raster
#' @param name dataset name without type prefix (p_|r_|s_)
#' @param method defaults to bilinear, or nearest neighbor (ngb)
#' @param debug defaults to False
#'
#' @return
#' @export
r_to_prjres <- function(prjres, r, name, method = c("bilinear", "ngb"), debug=F){
  # prjres = projections_tbl$prjres[2]; r = r_vgpm_0; name = "vgpm"; method = c("bilinear","ngb")

  P <- projections_tbl %>%
    filter(prjres == !!prjres)
  stopifnot(nrow(p) == 1)
  r_pr_tif <- glue("{dir_data}/{name}{prjres}.tif")
  message(basename(r_pr_tif))

  #devtools::load_all()
  r_pu_id_pr <- get_dataset_prj_res("r_pu_id", P$prj, P$res)

  r_pr <- suppressWarnings(
    projectRaster(
      from = r,
      to = r_pu_id_pr,
      method = method)) %>%
    mask(r_pu_id_pr)

  #plot(r_pu_id_pr, main = glue("r_pu_id{pr}"))
  #plot(r_pr, main = glue("r_vgpm{pr}"))
  #message(paste("  res:", res(r_pu_id_pr)))

  writeRaster(r_pr, r_pr_tif, overwrite=T)
}

#' Project stack to projection-raster
#'
#' @param prjres projection-resolution, per projections_tbl$prjres
#' @param s stack object
#' @param name dataset name without type prefix (p_|r_|s_)
#' @param method defaults to bilinear, or nearest neighbor (ngb)
#' @param debug defaults to False
#'
#' @return
#' @export
#'
#' @examples
s_to_prjres <- function(prjres, s, name, method = c("bilinear", "ngb"), debug=F){
  # prjres = projections_tbl$prjres[2]; r = r_vgpm_0; name = "vgpm"; method = c("bilinear","ngb")
  # s = s_bio_gmbi; name = "bio_gmbi"; method = c("bilinear","ngb")

  p <- projections_tbl %>%
    filter(prjres == !!prjres)
  stopifnot(nrow(p) == 1)

  dir_s <- glue("{name}{prjres}")
  if (debug) message(message(glue("dir_s: {dir_s}")))

  #devtools::load_all()
  r_pu_id_pr <- get_d_prjres("r_pu_id", prjres)

  s_pr <- suppressWarnings(
    projectRaster(
      from = s,
      to = r_pu_id_pr,
      method = method)) %>%
    mask(r_pu_id_pr)

  # write tifs
  map(names(s_pr), lyr_to_tif, s_pr, dir_s, debug)
}

#' Write layer in raster stack to tif
#'
#' @param lyr layer name
#' @param s stack
#' @param s_dir directory name
#' @param debug defaults to False
#'
#' @return
#' @export
#'
#' @examples
lyr_to_tif <- function(lyr, s, s_dir, debug=F){
  dir_tif <- here(glue("inst/data/{s_dir}"))
  if (!dir.exists(dir_tif)) dir.create(dir_tif)
  tif <- glue("{dir_tif}/{lyr}.tif")
  if (debug) message(glue("tif: {tif}"))
  r <- raster(s, lyr) #%>%
  #mask(r_pu_id)
  writeRaster(r, tif, overwrite=T)
}

#' Get projection info for tif of raster
#'
#' @param tif path to raster tif
#'
#' @return character string of projection_tbl$prjres
#' @export
get_tif_projection <- function(tif, debug=F){
  r <- raster(tif)
  get_r_projection(r)
}

#' Get projection info for raster
#'
#' @param r raster
#'
#' @return character string of projection_tbl$prjres
#' @export
get_r_projection <- function(r, debug=F){
  #tif  = list.files(dir_scenarios, "^s.*\\_sol.tif$", full.names=T)[1]

  # get projection with matching projection and closest resolution
  r_proj_str11 <- str_sub(as.character(crs(r)), end=11)
  r_res_num    <- res(r)[1]
  P <- projections_tbl %>%
    mutate(
      proj_str11  = str_sub(projections_tbl$proj, end=11),
      res_num_dif = abs(res_num - !!r_res_num)) %>%
    filter(
      proj_str11  == !!r_proj_str11,
      res_num_dif == min(res_num_dif))

  #pull(prjres)
  if (debug){
    message(glue(
      "tif: {basename(tif)}
       nrow(P): {nrow(P)}
       crs(r)[1-11]: {str_sub(as.character(crs(r)), end=11)}
       res(r)[1]: {res(r)[1]}"))
  }
  stopifnot(nrow(P) == 1)
  P
}

