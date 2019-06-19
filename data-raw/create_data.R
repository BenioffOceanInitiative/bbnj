# 1) use_data(). for lazy loading with library(bbnj). document datasets in R/data.R.
# 2) write tifs. for use with other software external to R, eg QGIS
library(tidyverse)
library(xml2)
library(purrr)
library(raster)
library(usethis)
library(sf)
library(glue)
library(here)
library(rmapshaper)
library(geojsonio)
library(bbnj)
library(rlang)
select = dplyr::select

library(gmbi) #devtools::install_github("marinebon/gmbi", force=T) #devtools::install_local("~/github/gmbi")
#remotes::install_local("~/github/gmbi", for)

devtools::load_all()

# paths ----
dir_data                 <- here("inst/data")
dir_gdata                <- "~/Gdrive Ecoquants/projects/bbnj/data" # on Ben Best's laptop
raw_eez_shp              <- glue("{dir_gdata}/raw/marineregions.org_boundaries/World_EEZ_v10_20180221/eez_v10.shp")
raw_eez_iho_shp          <- glue("{dir_gdata}/raw/marineregions.org_boundaries/Intersect_EEZ_IHO_v3_2018/EEZ_IHO_v3.shp")
raw_ppow_shp             <- glue("{dir_gdata}/raw/UNEP-WCMC/DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1/01_Data/WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp")
raw_fish_gfw_csv         <- glue("{dir_gdata}/raw/Sala et al 2018/half_degree_binned_results.csv")
raw_fish_saup_now_tif    <- glue("{dir_gdata}/raw/UBC-exploited-fish-projections/Current_MCP1.tif")
raw_fish_saup_future_tif <- glue("{dir_gdata}/raw/UBC-exploited-fish-projections/MCP2050_RCP85.tif")
raw_mine_claims_shp      <- glue("{dir_gdata}/raw/ISA_shapefiles/ISA_claim_areas_update_20181202.shp")
raw_phys_seamounts_kml   <- glue("{dir_gdata}/raw/Seamounts - Kim and Wessel 2011/KWSMTSv01.kml")
raw_phys_seamounts_txt   <- glue("{dir_gdata}/raw/Seamounts - Kim and Wessel 2011/KWSMTSv01.txt")
raw_phys_scapes_arcinfo  <- glue("{dir_gdata}/raw/Harris and Whiteway 2009/Global_Seascapes/class_11")
raw_phys_vents_csv       <- glue("{dir_gdata}/raw/Hydrothermal vents - Interridge Vent Database v3.4/vent_fields_all.csv")
abnj_shp                 <- glue("{dir_data}/abnj.shp")
abnj_s05_shp             <- glue("{dir_data}/abnj_s05.shp")
iho_shp                  <- glue("{dir_data}/iho.shp")
iho_s05_shp              <- glue("{dir_data}/iho_s05.shp")
ihor_shp                 <- glue("{dir_data}/ihor.shp")
ihor_s05_shp             <- glue("{dir_data}/ihor_s05.shp")
ihor_tif                 <- glue("{dir_data}/ihor.tif")
ppow_shp                 <- glue("{dir_data}/ppow.shp")
ppow_s05_shp             <- glue("{dir_data}/ppow_s05.shp")
eez_shp                  <- glue("{dir_data}/eez.shp")
eez_s05_shp              <- glue("{dir_data}/eez_s05.shp")
pu_id_tif                <- glue("{dir_data}/pu_id.tif")
mine_claims_shp          <- glue("{dir_data}/mine-claims.shp")
mine_claims_tif          <- glue("{dir_data}/mine-claims.tif")
phys_seamounts_tif       <- glue("{dir_data}/phys_seamounts.tif")
phys_vents_tif           <- glue("{dir_data}/phys_vents.tif")

# variables ----
redo_eez  = F
redo_abnj = F
redo_ihor = F
redo_lyrs = F

# helper functions ----
lyr_to_tif <- function(lyr, s, pfx){
  dir_tif <- here(glue("inst/data/{pfx}"))
  if (!dir.exists(dir_tif)) dir.create(dir_tif)
  tif <- glue("{dir_tif}/{lyr}.tif")
  r <- raster(s, lyr) %>%
    mask(r_pu_id)
  writeRaster(r, tif, overwrite=T)
}

# globe
bb = st_sf(
  tibble(
    one = 1,
    geom = st_sfc(st_polygon(list(rbind(
      c(-180, -90),
      c( 180, -90),
      c( 180,  90),
      c(-180,  90),
      c(-180, -90)))), crs = 4326)))

# p_eez ----
if (!file.exists(eez_shp) | redo_eez){
  p_eez <- read_sf(raw_eez_shp)
  #write_sf(p_eez, eez_shp)  # TOO BIG: 159.9 MB
  #use_data(p_eez, overwrite = TRUE) # TOO BIG: 121.8 MB

  p_eez_s05 <- p_eez %>%
    geojson_json() %>% # convert to geojson for faster ms_simplify; 69.6 sec
    ms_simplify(keep=0.05, keep_shapes=T)
  # TODO: crashing on BB's Mac 2019-03-27 so using prior version
  # p_eez_s05 <- read_sf(eez_s05_shp)

  use_data(p_eez_s05, overwrite = TRUE)
  write_sf(p_eez_s05, eez_s05_shp)
}

# p_abnj, p_iho ----
if (!file.exists(abnj_s05_shp) | redo_abnj){
  eez_iho <- read_sf(raw_eez_iho_shp)

  p_iho <- eez_iho %>%
    filter(is.na(EEZ)) %>%
    select(fid, MarRegion, MRGID, IHO_Sea, IHO_MRGID, Longitude, Latitude, Area_km2)

  abnj <- p_iho %>%
    mutate(name="Areas Beyond National Jurisdiction, i.e. High Seas") %>%
    group_by(name) %>%
    summarize()

  abnj_parts <- st_cast(abnj, "POLYGON") %>%
    mutate(
      area_km2 = st_area(geometry) %>%  units::set_units(km^2))
  #write_sf(abnj_parts, glue("{dir_gdata}/derived/boundary/abnj_parts.shp"))
  # manually edited out islands

  abnj_parts <- read_sf(glue("{dir_gdata}/derived/boundary/abnj_parts.shp"))
  p_abnj <- st_union(abnj_parts) %>%
    st_intersection(bb)
  write_sf(p_abnj, abnj_shp)
  p_abnj <- read_sf(abnj_shp)
  use_data(p_abnj, overwrite = TRUE)

  # TODO: clip iho to high seas
  #p_iho <- st_intersection(p_iho, p_abnj)
  write_sf(p_iho, iho_shp)
  use_data(p_iho, overwrite = TRUE)

  #p_iho <- read_sf(iho_shp)
  p_iho_s05 <- rmapshaper::ms_simplify(p_iho, keep = 0.05)
  use_data(p_iho_s05, overwrite = TRUE)
  write_sf(p_iho_s05, iho_s05_shp)

  p_abnj_s05 <- rmapshaper::ms_simplify(p_abnj, keep = 0.05)
  use_data(p_abnj_s05, overwrite = TRUE)
  write_sf(p_abnj_s05, abnj_s05_shp)

  p_abnj_s05 <- rmapshaper::ms_simplify(p_abnj, keep = 0.05)
  use_data(p_abnj_s05, overwrite = TRUE)
  write_sf(p_abnj_s05, abnj_s05_shp)
}

# p_ihor, s_ihor ----
# raster (and polygon) of IHO seas revised (ihor), ie 7 seas
if (!file.exists(ihor_s05_shp) | redo_ihor){

  p_iho <- read_sf(iho_shp)

  sz <- 2500000 # small (sm) < 2.5M km2 <= large (lg)
  sm <- filter(p_iho, Area_km2  < sz)
  lg <- filter(p_iho, Area_km2 >= sz)

  p_ihor <- rbind(
    lg,
    sm %>%
      select(-IHO_Sea) %>%
      left_join(
        st_join(
          st_centroid(sm) %>% select(fid),
          lg %>% select(IHO_Sea),
          join = st_nearest_feature) %>%
          st_drop_geometry(),
        by = "fid"))

  p_ihor <- p_ihor %>%
    group_by(IHO_Sea) %>%
    summarize(
      Area_km2 = sum(Area_km2)) %>%
    # Warning message:
    #   In st_is_longlat(x) :
    #     bounding box has potentially an invalid value range for longlat data
    st_crop(bb) %>%
    mutate(
      area_km2 = st_area(geometry) %>% units::set_units(km2)) %>%
    arrange(IHO_Sea) %>%
    tibble::rowid_to_column("seaid") %>%
    select(seaid, sea=IHO_Sea, area_km2)

  write_sf(p_ihor, ihor_shp)
  p_ihor <- read_sf(ihor_shp)
  use_data(p_ihor, overwrite = TRUE)

  p_ihor_s05 <- rmapshaper::ms_simplify(p_ihor, keep = 0.05)
  use_data(p_ihor_s05, overwrite = TRUE)
  write_sf(p_ihor_s05, ihor_s05_shp)

  r_ihor <- rasterize(p_ihor, r_pu_id, field="seaid") %>%
    mask(r_pu_id) # plot(r_ihor)
  # OLD: all in one raster
  # names(r_ihor) <- "seaid"
  # r_ihor <- ratify(r_ihor)
  # levels(r_ihor) <- p_ihor %>%
  #   st_drop_geometry() %>%
  #   select(ID=seaid, sea, area_km2) %>%
  #   as.data.frame()
  #factorValues(r_ihor, 1:7)
  # writeRaster(s_ihor, ihor_tif, overwrite = TRUE)
  # use_data(r_ihor, overwrite = TRUE)

  # NEW: seperate each sea into own raster to use as representative target
  seas = c("Arctic", "Indian", "N_Atlantic", "N_Pacific", "S_Atlantic", "S_Pacific", "Southern")
  s_ihor = stack()
  for (i in 1:length(seas)){ # i = 1
    s_ihor = stack(s_ihor, r_ihor == i)
    names(s_ihor)[i] = seas[i]
    # plot(raster(s_ihor, i), main = seas[i])
  }

  # write tifs
  map(names(s_ihor), lyr_to_tif, s_ihor, "bnd_ihor")

  # load into memory so not referencing local file and use_data() works
  s_ihor <- raster::readAll(s_ihor)

  # use_data()
  use_data(s_ihor, overwrite = TRUE)
}

# p_ppow ----
if (!file.exists(ppow_s05_shp) | redo_lyrs){
  ppow <- read_sf(raw_ppow_shp)
  abnj <- read_sf(abnj_shp)

  p_ppow <- st_intersection(ppow, abnj) %>%
    select(-FID, -one) %>%
    mutate(
      area_km2 = st_area(geometry) %>%  units::set_units(km^2) %>% as.integer())
  write_sf(p_ppow, ppow_shp)
  p_ppow <- read_sf(ppow_shp)
  use_data(p_ppow, overwrite = TRUE)

  p_ppow_s05 <- rmapshaper::ms_simplify(p_ppow, keep = 0.05)
  use_data(p_ppow_s05, overwrite = TRUE)
  write_sf(p_ppow_s05, ppow_s05_shp)

  # TODO: p_ppow: merge small parts into larger
  # TODO: r_ppow: assign numeric id for raster
}

# r_na ----
# for assigning fresh values later
r_na <- raster(
  xmn = -180, xmx = 180, ymn = -90, ymx = 90,
  resolution=0.5, crs=leaflet:::epsg4326, vals=NA)

# r_pu_id ----
if (!file.exists(pu_id_tif) | redo_lyrs){
  r_pu_id <- r_na
  values(r_pu_id) <- 1:ncell(r_na)

  r_abnj  <- rasterize(abnj, r_na)
  r_pu_id <- mask(r_pu_id, r_abnj)

  writeRaster(r_pu_id, pu_id_tif, overwrite=T)
  use_data(r_pu_id, overwrite = TRUE)
}

# s_bio_gmbi ----
if (!dir.exists("inst/data/bio_gmbi") | redo_lyrs){
  # get bio raster stack from gmbi
  data(gmbi_stack)

  # mask stack
  s_bio_gmbi <- gmbi_stack %>%
    mask(r_pu_id)

  # remove layers
  mins <- cellStats(s_bio_gmbi, "min")
  #length(names(s_bio_gmbi))
  s_bio_gmbi <- raster::subset(s_bio_gmbi, names(s_bio_gmbi)[mins!=Inf])

  # write tifs
  map(names(s_bio_gmbi), lyr_to_tif, s_bio_gmbi, "bio_gmbi")

  # load into memory so not referencing local file and use_data() works
  s_bio_gmbi <- raster::readAll(s_bio_gmbi)

  # use_data()
  use_data(s_bio_gmbi, overwrite = TRUE)
}

# s_bio_gmbi: seagrass? ----

# r_sg <- raster(s_bio_gmbi, "nspp_seagrasses")
# #xyFromCell(object, cell, spatial=FALSE, ...)
# ?raster
#   !is.na(r_sg)
# sum(!is.na(values(raster(s_bio_gmbi, "nspp_seagrasses"))))


# s_fish_gfw ----
if (!dir.exists("inst/data/fish_gfw") | T){

  # read csv, identify cell_id
  d_fish_gfw <- read_csv(raw_fish_gfw_csv) %>%
    mutate(
      cell_id = cellFromXY(r_pu_id, matrix(data = c(lon_bin_center, lat_bin_center), ncol=2)))

  # create raster stack
  if (exists("s_fish_gfw")) rm(s_fish_gfw)
  lyrs <- names(d_fish_gfw)[4:(ncol(d_fish_gfw) - 1)]
  for (lyr in lyrs){ # lyr = lyrs[2]
    r <- r_na
    names(r) <- lyr
    r[d_fish_gfw$cell_id] <- d_fish_gfw[[lyr]]

    #if (!exists("s_fish_gfw")){
    if (!env_has(global_env(), "s_fish_gfw")){
      s_fish_gfw <- r
    } else {
      s_fish_gfw <- stack(s_fish_gfw, r)
    }
  }

  # mask stack
  s_fish_gfw <- mask(s_fish_gfw, r_pu_id)

  # write tifs
  map(lyrs, lyr_to_tif, s_fish_gfw, "fish_gfw")

  # use_data()
  use_data(s_fish_gfw, overwrite = TRUE)
}

# s_fish_saup ----
if (!dir.exists("inst/data/fish_saup") | redo_lyrs){

  s_fish_saup  <- stack(c(raw_fish_saup_now_tif, raw_fish_saup_future_tif))
  lyrs <- names(s_fish_saup) <- c("mcp_2004", "mcp_2050")

  # mask stack
  s_fish_saup <- mask(s_fish_saup, r_pu_id)

  # write tifs
  map(lyrs, lyr_to_tif, s_fish_saup, "fish_saup")

  # use_data()
  use_data(s_fish_saup, overwrite = TRUE)
}

# s_phys_scapes ----
if (!dir.exists("inst/data/phys_scapes") | T){

  # load raster in 0.1 deg resolution
  r_scapes <- raster(raw_phys_scapes_arcinfo)

  # create raster stack of benthicscape classes
  # with each layer containing amount of class in 0.5 deg cell
  if (exists("s_scapes")) rm(s_scapes)
  r_a <- area(r_scapes)
  for (i in 1:11){ # i = 2
    r_i <- (r_scapes == i) * r_a
    r_ia <- aggregate(r_i, fact=5, fun=sum, expand=F, na.rm=T)
    #plot(r_ia, col=cols)

    #if (!exists("s_scapes")){
    if (!env_has(global_env(), "s_scapes")){
      s_scapes <- r_ia
    } else {
      s_scapes <- stack(s_scapes, r_ia)
    }
  }
  lyrs <- names(s_scapes) <- glue("class_{1:11}")

  # mask stack
  s_phys_scapes <- s_scapes <- mask(s_scapes, r_pu_id)

  # write tifs
  map(lyrs, lyr_to_tif, s_phys_scapes, "phys_scapes")

  # use_data()
  use_data(s_phys_scapes, overwrite = TRUE)
}

# r_phys_seamounts ----
if (!file.exists(phys_seamounts_tif) | redo_lyrs){


  x <- read_xml(raw_phys_seamounts_kml)

  xpaths <- list(
    name = "/d1:kml/d1:Document/d1:Folder/d1:Folder/d1:name",
    xyz  = "/d1:kml/d1:Document/d1:Folder/d1:Folder/d1:Placemark/d1:Point/d1:coordinates")
  pts <- tibble(
    names = xml_find_all(x, xpaths$name) %>% xml_text(),
    xyz   = xml_find_all(x, xpaths$xyz) %>% xml_text()) %>%
    separate(xyz, c("x", "y", "z"), sep=",", convert=T) %>%
    st_as_sf(coords = c("x", "y"), crs = 4326)

  readr::read_delim()
  x <- read_csv(raw_phys_seamounts_txt)

  pts <- st_as_sf(x, coords = c("x", "y"), crs = 4326)

  r_phys_seamounts <- rasterize(st_coordinates(pts), r_pu_id, fun='count', background=0) %>%
    mask(r_pu_id) # plot(r_phys_seamounts)
  names(r_phys_seamounts) <- "count"

  writeRaster(r_phys_seamounts, phys_seamounts_tif, overwrite = TRUE)
  use_data(r_phys_seamounts, overwrite = TRUE)
}

# r_phys_vents ----
if (!file.exists(phys_vents_tif) | redo_lyrs){
  pts_vents <- read_csv(raw_phys_vents_csv) %>%
    st_as_sf(
      coords = c("Longitude", "Latitude"), crs = 4326, remove=F)

  r_phys_vents <- rasterize(
    st_coordinates(pts_vents), r_pu_id, fun='count', background=0) %>%
    mask(r_pu_id) # plot(r_phys_seamounts)
  names(r_phys_vents) <- "count"

  writeRaster(r_phys_vents, phys_vents_tif, overwrite = TRUE)
  use_data(r_phys_vents, overwrite = TRUE)
}

# r_mine_claims ----
if (!file.exists(mine_claims_tif) | redo_lyrs){

  p_mine_claims <- read_sf(raw_mine_claims_shp) %>%
    # remove: APEI (Area of Particular Environmental Interest)
    filter(area_type != "apei")
  #table(p_mine_claims$area_type)
  # apei    claim reserved
  #    9     1330      174
  #table(p_mine_claims$Region)
  # Atlantic Ocean  Pacific Ocean
  #            144            440
  #table(p_mine_claims$Min_type)
  #   SMS  crust nodule
  #   700    725     79
  #table(p_mine_claims$holder)
  #   BGR (Germany)     CIIC (Cook Islands)
  #             102                       3
  #     CMC (China)           COMRA (China)
  #               8                     255
  # ...
  write_sf(p_mine_claims, mine_claims_shp)
  use_data(p_mine_claims, overwrite = TRUE)

  r_mine_claims <- rasterize(p_mine_claims, r_pu_id, field=1, background=0) %>%
    mask(r_pu_id) # plot(r_mine_claims)
  names(r_mine_claims) <- "present"
  writeRaster(r_mine_claims, mine_claims_tif, overwrite = TRUE)
  use_data(r_mine_claims, overwrite = TRUE)
}
