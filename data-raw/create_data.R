# 1) use_data(). for lazy loading with library(bbnj). document datasets in R/data.R.
# 2) write tifs. for use with other software external to R, eg QGIS
library(tidyverse)
library(xml2)
library(purrr)
library(raster)
library(fasterize)
library(rnaturalearth)
library(usethis)
library(sf)
library(glue)
library(here)
library(rmapshaper)
library(geojsonio)
library(rlang)
select = dplyr::select

library(gmbi) #devtools::install_github("marinebon/gmbi", force=T) #devtools::install_local("~/github/gmbi")
#remotes::install_local("~/github/gmbi", for)

#library(bbnj) #devtools::install_github("ecoquants/bbnj", force=T) #devtools::install_local("~/github/bbnj")
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
raw_vgpm_dir             <- glue("{dir_gdata}/raw/VGPM_oregonstate.edu")
countries_shp            <- glue("{dir_data}/countries.shp")
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
scapes_tif               <- sprintf("%s/class_11.tif", raw_phys_scapes_arcinfo %>% dirname() %>% dirname())

# variables ----
redo_eez  = F
redo_abnj = F
redo_ihor = F
redo_lyrs = F
redo_project_polygons = F
redo_project_pu_id_tifs = F

# helper functions ----

# globe
bb <- get_global_bb(crs=4326)

# projections ----

# r_gcs     <- get_grid(val=1)
# r_gcs_km2 <- area(r_gcs)
# c_gcs_km2 <- cellStats(r_gcs_km2, stat='mean')
# c_gcs_km  <- c_gcs_km2^0.5 # 44.28597 # 50

# list of projections, possibly with multiple resolutions
projections_lst <- list(
  gcs  =
    list(
      default = T,
      name = "Geographic",
      proj = leaflet:::epsg4326,
      epsg = 4326,
      res_num = list(
        "0.5d" = 0.5)),
  mer  =
    list(
      default = F,
      name = "Mercator",
      proj = leaflet:::epsg3857,
      epsg = 3857,
      res_num = list(
        #"56x16km" = c(55659.75, 156862.8))),
        "36km" = 36000)), # split difference of original non-rectilinear cell resolution
  mol =
    list(
      default = F,
      name = "Mollweide",
      proj = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
      epsg = 54009,
      res_num  = list(
        "50km"  =  50000,
        "10km"  =  10000,
        "100km" = 100000)))

# table of projections, flattening nested resolutions
projections_tbl <- map_df(
  projections_lst,
  function(x){
    as_tibble(x) %>% unnest(res_num, .id="res")},
  .id="prj") %>%
  mutate(
    prjres = ifelse(
      prj == "gcs" & res == "0.5d",
      "",
      glue("_{prj}{res}")))
use_data(projections_lst, overwrite = TRUE)
use_data(projections_tbl, overwrite = TRUE)

# p_eez ----
if (!file.exists(eez_s05_shp) | redo_eez){
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

  #cat(paste(p_ihor$sea, collapse='` = "",\n`'))
  sea_keys <- c(
    `Arctic Ocean`         = "Arctic",
    `Indian Ocean`         = "Indian",
    `North Atlantic Ocean` = "N_Atlantic",
    `North Pacific Ocean`  = "N_Pacific",
    `South Atlantic Ocean` = "S_Atlantic",
    `South Pacific Ocean`  = "S_Pacific",
    `Southern Ocean`       = "Southern")

  # make alternate polygon for regions
  # see section: p_*_[prj] for non-gcs projections
  # make rasters in projections
  for (prjres in projections_tbl$prjres){ # prjres = "_mol50km"
    dir_s   <- glue("bnd_ihor{prjres}")

    #projections_lst$mer$proj
    r_pu_id <- get_d_prjres("r_pu_id", prjres)
    p_ihor  <- get_d_prjres("p_ihor", prjres) #%>%

    r_ihor <- fasterize(p_ihor, r_pu_id, field="seaid") %>%
      mask(r_pu_id)

    s_ihor <- prioritizr::binary_stack(r_ihor)
    names(s_ihor) <- sea_keys

    map(names(s_ihor), lyr_to_tif, s_ihor, dir_s)

    if (prjres == ""){
      s_ihor <- raster::readAll(s_ihor)
      use_data(s_ihor, overwrite = TRUE)
    }
  }
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

# p_countries ----
p_countries <- rnaturalearth::ne_countries(
  scale = 110, type = "countries", returnclass = "sf") %>%
  mutate(
    pop_est = pop_est / 1000) %>%
  rename(pop_est_k = pop_est)

write_sf(p_countries, countries_shp)
p_countries <- read_sf(countries_shp)
use_data(p_countries, overwrite = TRUE)

# p_*_[prj] for non-gcs projections ----
shps_tbl <- tribble(
  ~pfx       , ~shp,
  "countries", countries_shp,
  "eez_s05"  , eez_s05_shp,
  "abnj"     , abnj_shp,
  "abnj_s05" , abnj_s05_shp,
  "ihor"     , ihor_shp,
  "ihor_s05" , ihor_s05_shp,
  "ppow"     , ppow_shp,
  "ppow_s05" , ppow_s05_shp)

# iterate over new projections (prj, res)
prjs <- projections_tbl %>%
  group_by(prj) %>%
  summarize(
    epsg = first(epsg))

if (redo_project_polygons){
  for (i in 2:length(projections_lst)){ # i=3
    prj  <- names(projections_lst)[i]
    epsg <- projections_lst[[i]]$epsg

    # iterate over shapefiles
    for (j in 1:nrow(shps_gcs)){ # j=1
      s <- shps_tbl[j,]

      o_shp <- glue("{dirname(s$shp)}/{s$pfx}_{prj}.shp")

      read_sf(s$shp, quiet=T) %>%
        st_transform(o, crs = epsg) %>%
        st_write(o_shp, delete_layer = T)
    }
  }
}

# r_na ----
# for assigning fresh values later
r_na <- raster(
  xmn = -180, xmx = 180, ymn = -90, ymx = 90,
  resolution=0.5, crs=projections_lst$gcs$proj, vals=NA)
#r_na_mer <- leaflet::projectRasterForLeaflet(r_na, "ngb")

# r_pu_id ----
if (!file.exists(pu_id_tif) | redo_lyrs){
  r_pu_id <- r_na
  values(r_pu_id) <- 1:ncell(r_na)

  r_abnj <- read_sf(abnj_shp) %>%
    fasterize(r_na) # rasterize() resulted in horizontal slivers

  r_pu_id <- mask(r_pu_id, r_abnj) # plot(r_pu_id)

  writeRaster(r_pu_id, pu_id_tif, overwrite=T)
  use_data(r_pu_id, overwrite = TRUE)
}

# pu_id_[prj][res].tif ----
if (redo_project_pu_id_tifs){

  for (prjres in projections_tbl$prjres){ # prjres = "" # _mer36km _mer36km
    P <- projections_tbl %>% filter(prjres == !!prjres)
    pu_id_pr_tif <- glue("{dir_data}/pu_id{prjres}.tif")
    message(glue("{prjres}: {pu_id_pr_tif}"))

    r_na_pr <- suppressWarnings(raster::projectRaster(
      r_na, raster::projectExtent(r_na, crs = sp::CRS(P$proj)),
      res = P$res_num))

    r_abnj_pr <- get_d_prjres("p_abnj", prjres) %>%
      fasterize(r_na_pr)     # rasterize() resulted in horizontal slivers
    crs(r_abnj_pr) <- P$proj # presume same projection for raster

    r_pu_id_pr <- r_na_pr %>%
      setValues(1:ncell(r_na_pr)) %>%
      mask(r_abnj_pr)
    #mapview::mapview(r_pu_id_pr)

    if (prjres == ""){
      r_pu_id <- r_pu_id_pr
      use_data(r_pu_id, overwrite=T)
    }

    writeRaster(r_pu_id_pr, pu_id_pr_tif, overwrite=T)
  }
}

# r_vpgm ----
if (!file.exists(vgpm_tif) | redo_lyrs){

  # Downloaded vgpm.v.[YYYY].xyz.tar from [Ocean Productivity: Online Data: Standard VGPM - 2160 by 4320 Monthly XYZ files from VIIRS Data](http://orca.science.oregonstate.edu/2160.by.4320.monthly.xyz.vgpm.v.chl.v.sst.php)
  # Feb 1 2013 - Jan 31 2019

  # untar and ungzip
  tars <- list.files(raw_vgpm_dir, ".tar", full.names=T)
  lapply(tars, untar, exdir = path.expand(raw_vgpm_dir))
  gzs <- list.files(raw_vgpm_dir, ".gz", full.names=T)
  lapply(gzs, R.utils::gunzip)

  # xyz to raster, in parallel
  source(here("data-raw/vgpm_func.R"))
  xyzs <- list.files(raw_vgpm_dir, ".xyz", full.names=T)
  library(furrr)
  plan(multiprocess)
  out <- future_map(xyzs, vgpm.raster)
  # vgpm.raster(file.path(raw_vgpm_dir, "vgpm.2019001.all.xyz"))

  tifs <- list.files(raw_vgpm_dir, ".*\\.tif$", full.names = T)
  r_vgpm_0 <- stack(tifs) %>%
    mean(na.rm=T) #%>%
  #aggregate(fact=6) %>%
  #mask(r_pu_id)
  #plot(log(r_vgpm_0))

  #for (prjres in projections_tbl$prjres){ # prjres = projections_tbl$prjres[1]
  for (prjres in projections_tbl$prjres){ # prjres = projections_tbl$prjres[4]

    vgpm_tif   <- glue("{dir_data}/vgpm{prjres}.tif")
    r_pu_id_pr <- get_d_prjres("r_pu_id", prjres)

    if (prjres == ""){
      # gcs 0.5d works out to exactly 6 cells
      r_vgpm <- r_vgpm_0 %>%
        aggregate(fact=6, fun=mean) %>%
        mask(r_pu_id_pr)
      #plot(log(r_vgpm))

      use_data(r_vgpm, overwrite=T)
    } else {
      r_vgpm <- projectRaster(
        from = r_vgpm_0,
        to = r_pu_id_pr,
        method = "bilinear") %>%
        mask(r_pu_id_pr)
      # TODO: check above since previously created int'l dateline gap when mapping in leaflet
      # preferred projectRaster():
          # prjres <- "_mer36km"
          # P <- projections_tbl %>% filter(prjres == !!prjres)
          # r_pu_id_pr <- get_d_prjres("r_pu_id", prjres)
          # s_features <- suppressWarnings(
          #   raster::projectRaster(
          #     s_features, raster::projectExtent(s_features, crs = sp::CRS(P$proj)),
          #     res = P$res_num)) %>%
          #   mask(r_pu_id_pr)
    }

    # write gcs raster for general use
    writeRaster(r_vgpm, vgpm_tif, overwrite=T)
  }
}

# s_bio_gmbi ----
if (!dir.exists("inst/data/bio_gmbi") | redo_lyrs){
  # get bio raster stack from gmbi
  data(gmbi_stack)
  data("r_pu_id")

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

  # create stacks in other projections
  walk(
    projections_tbl$prjres[-1],
    s_to_prjres, s_bio_gmbi, "bio_gmbi", "bilinear", debug=F)
}

# s_bio_gmbi: seagrass? ----

# r_sg <- raster(s_bio_gmbi, "nspp_seagrasses")
# #xyFromCell(object, cell, spatial=FALSE, ...)
# ?raster
#   !is.na(r_sg)
# sum(!is.na(values(raster(s_bio_gmbi, "nspp_seagrasses"))))


# s_fish_gfw ----
if (!dir.exists("inst/data/fish_gfw") | redo_lyrs){

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

  # load into memory so not referencing local file and use_data() works
  s_fish_gfw <- raster::readAll(s_fish_gfw)

  # use_data()
  use_data(s_fish_gfw, overwrite = TRUE)

  # create stacks in other projections
  walk(
    projections_tbl$prjres[-1],
    s_to_prjres, s_fish_gfw, "fish_gfw", "bilinear", debug=F)
}

# s_fish_saup ----
if (!dir.exists("inst/data/fish_saup") | redo_lyrs){

  s_fish_saup  <- stack(c(raw_fish_saup_now_tif, raw_fish_saup_future_tif))
  lyrs <- names(s_fish_saup) <- c("mcp_2004", "mcp_2050")

  # mask stack
  s_fish_saup <- mask(s_fish_saup, r_pu_id)

  # write tifs
  map(lyrs, lyr_to_tif, s_fish_saup, "fish_saup")

  # load into memory so not referencing local file and use_data() works
  s_fish_saup <- raster::readAll(s_fish_saup)

  # use_data()
  use_data(s_fish_saup, overwrite = TRUE)

  # create stacks in other projections
  walk(
    projections_tbl$prjres[-1],
    s_to_prjres, s_fish_saup, "fish_saup", "bilinear", debug=F)
}

# s_phys_scapes ----
if (!dir.exists("inst/data/phys_scapes") | redo_lyrs){

  # load raster in 0.1 deg resolution
  r_scapes <- raster(raw_phys_scapes_arcinfo)
  writeRaster(r_scapes, scapes_tif)
  r_scapes <- raster(scapes_tif)

  for (prjres in projections_tbl$prjres){ # prjres = ""
    #P <- projections_tbl %>% filter(prjres == !!prjres)
    r_pu_id_pr <- get_d_prjres("r_pu_id", prjres)
    r_scapes_pr <- r_scapes %>%
      projectRaster(r_pu_id_pr, method = "ngb") %>%
      mask(r_pu_id_pr)

    s_phys_scapes        <- prioritizr::binary_stack(r_scapes_pr)
    names(s_phys_scapes) <- glue("class_{1:11}")

    # write tifs
    map(names(s_phys_scapes), lyr_to_tif, s_phys_scapes, glue("phys_scapes{prjres}"))

    if (prjres == ""){
      # load into memory so not referencing local file and use_data() works
      s_phys_scapes <- raster::readAll(s_phys_scapes)

      # use_data()
      use_data(s_phys_scapes, overwrite = TRUE)
    }
  }

  # # create stacks in other projections
  # walk(
  #   projections_tbl$prjres[-1],
  #   s_to_prjres, s_phys_scapes, "phys_scapes", "bilinear", debug=F)
}

# s_phys_seamounts ----
if (!file.exists(phys_seamounts_tif) | redo_lyrs){

  # OLD: read in XML since read_sf() was taking days
  # x <- read_xml(raw_phys_seamounts_kml)
  #
  # xpaths <- list(
  #   name = "/d1:kml/d1:Document/d1:Folder/d1:Folder/d1:name",
  #   xyz  = "/d1:kml/d1:Document/d1:Folder/d1:Folder/d1:Placemark/d1:Point/d1:coordinates")
  # pts <- tibble(
  #   names = xml_find_all(x, xpaths$name) %>% xml_text(),
  #   xyz   = xml_find_all(x, xpaths$xyz) %>% xml_text()) %>%
  #   separate(xyz, c("x", "y", "z"), sep=",", convert=T) %>%
  #   st_as_sf(coords = c("x", "y"), crs = 4326)

  # NEW TODO: easier to read CSV
  # Global Seamount Database from Kim, S.-S., and P. Wessel (2011), Geophys. J. Int., in revision.
  # This file contains 24,646 potential seamounts from the entire ocean basins.
  # The columns are:
  # Col 1: Longitude (-180/+180). Center of each seamount (in degrees)
  # Col 2: Latitude (-90/+90). Center of each seamount (in degrees)
  # Col 3: Estimated azimuth of the basal ellipse (in degree)
  # Col 4: Estimated major axis of the basal ellipse of each seamount (in km)
  # Col 5: Estimated minor axis of the basal ellipse of each seamount (in km)
  # Col 6: Seamount height obtained from the prediced bathymetry TOPO V12 (in m)
  # Col 7: Maximum amplitude of the Free-Air Gravity Anomaly (in mGal)
  # Col 8: Maximum amplitude of the Vertical Gravity Gradient Anomaly (in Eotvos)
  # Col 9: Regional depth of each seamount (in m)
  # Col 10: Age of underlying seafloor from the AGE 3.2 grid (in Myr)
  # Col 11: ID for each seamount (plate_###)
  d <- read_tsv(
    raw_phys_seamounts_txt, skip=18,
    col_names = c("lon", "lat", "azimuth_deg", "axis_major_km", "axis_minor_km", "height_m", "max_faa_mgal", "max_vgg_eotvos", "depth_m", "crustage_myr", "id"))

  brks <- c(0, 200, 800, Inf)
  lbls <- c("lteq200m","gt200lteq800m","gt800m")

  pts <- d %>%
    # remove odd lines like "> AF 3888 seamounts" (original line 18),"> AN 4837 seamounts" (imported line 3889)
    slice(setdiff(1:nrow(d), problems(d) %>% pull(row))) %>%
    mutate(
      lon = as.numeric(lon),
      lat = as.numeric(lat),
      summit_depth_m     = (-1 * depth_m) - height_m,
      summit_depth_m_brk = cut(summit_depth_m, brks)) %>%
    rowid_to_column() %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove=F) %>%
    filter(summit_depth_m > 0) # presumably islands: remove 158 pts of 24,643

  # range(pts$depth_m)          #  -7699.2   -90.9
  # range(pts$height_m)         #    100.0  6593.2
  # range(pts$summit_depth_m)   #  -2783.7  5962.3
  # sum(pts$summit_depth_m < 0) # 158 pts of 24,643

  # # scatter plot
  # p <- ggplot(data = pts %>% st_drop_geometry(), aes(x = depth_m, y = summit_depth_m, label = rowid)) +
  #   geom_point()
  # plotly::ggplotly(p) # slow b/c many pts
  #
  # # before: filter(summit_depth_m < 0)
  # pts_q <- pts %>%
  #   #st_drop_geometry() %>%
  #   filter(summit_depth_m < 0) %>%
  #   arrange(summit_depth_m) # %>% View()
  #
  # # map all islands and zoom to specific ones
  # library(mapview)
  # names(leaflet::providers) %>% sort()
  # mapviewOptions(basemaps = c("Esri.OceanBasemap", "Esri.WorldTopoMap"))
  #
  # pt <- pts_q %>%
  #   filter(rowid == 22358)
  # mapview(pts_q, zcol = "summit_depth_m", legend=T) %>%
  #   leafem::addMouseCoordinates() %>%
  #   leaflet::flyTo(pt$lon, pt$lat, zoom = 10)

  r_lvl <- function(brk, prjres){ # brk="(800,Inf]"
    r_pu_id_pr <- get_d_prjres("r_pu_id", prjres)
    epsg <- projections_tbl %>% filter(prjres==!!prjres) %>% pull(epsg)
    p <- pts %>%
      filter(summit_depth_m_brk == brk) %>%
      st_transform(epsg) %>%
      st_coordinates()
    r <- rasterize(p, r_pu_id_pr, fun='count', background=0) %>%
      mask(r_pu_id_pr)
    #plot(r)
    r
  }

  lvls <- levels(pts$summit_depth_m_brk)

  for (prjres in projections_tbl$prjres){ # prjres = projections_tbl$prjres[1]
    #for (prjres in projections_tbl$prjres[-1]){ # prjres = projections_tbl$prjres[1]

    s_phys_seamounts        <- stack(sapply(lvls, r_lvl, prjres))
    names(s_phys_seamounts) <- lbls

    # write tifs
    map(lbls, lyr_to_tif, s_phys_seamounts, glue("phys_seamounts{prjres}"))

    if (prjres == ""){
      # load into memory so not referencing local file and use_data() works
      s_phys_seamounts <- raster::readAll(s_phys_seamounts)

      # use_data()
      use_data(s_phys_seamounts, overwrite = TRUE)
    }
  }
}

# r_phys_vents ----
if (!file.exists(phys_vents_tif) | redo_lyrs){
  pts_vents <- read_csv(raw_phys_vents_csv) %>%
    st_as_sf(
      coords = c("Longitude", "Latitude"), crs = 4326, remove=F)

  for (prjres in projections_tbl$prjres){ # prjres = projections_tbl$prjres[1]

    phys_vents_tif <- glue("{dir_data}/phys_vents{prjres}.tif")
    r_pu_id_pr     <- get_d_prjres("r_pu_id", prjres)

    r_phys_vents <- rasterize(
      st_coordinates(pts_vents), r_pu_id_pr, fun='count', background=0) %>%
      mask(r_pu_id_pr) # plot(r_phys_seamounts)
    names(r_phys_vents) <- "count"

    writeRaster(r_phys_vents, phys_vents_tif, overwrite = TRUE)

    if (prjres == ""){
      use_data(r_phys_vents, overwrite = TRUE)
    }
  }
}

# r_mine_claims ----
if (!file.exists(mine_claims_tif) | redo_lyrs){

  p_mine_claims <- read_sf(raw_mine_claims_shp) %>%
    # remove: APEI (Area of Particular Environmental Interest)
    filter(area_type != "apei") %>%
    lwgeom::st_make_valid()
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

  p_mine_claims <- read_sf(mine_claims_shp) %>%
    lwgeom::st_make_valid()
  for (prjres in projections_tbl$prjres){ # prjres = projections_tbl$prjres[1]
    #for (prjres in projections_tbl$prjres[-1]){ # prjres = projections_tbl$prjres[2]

    mine_claims_tif  <- glue("{dir_data}/mine-claims{prjres}.tif")
    r_pu_id_pr       <- get_d_prjres("r_pu_id", prjres)

    r_mine_claims <- rasterize(p_mine_claims, r_pu_id_pr, field=1, background=0) %>%
      mask(r_pu_id_pr) # plot(r_mine_claims)
    names(r_mine_claims) <- "present"
    writeRaster(r_mine_claims, mine_claims_tif, overwrite = TRUE)

    if (prjres == ""){
      use_data(r_mine_claims, overwrite = TRUE)
    }
  }

}
