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
select = dplyr::select

library(gmbi) #devtools::install_github("marinebon/gmbi") #devtools::install_local("~/github/gmbi")
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
raw_phys_scapes_arcinfo  <- glue("{dir_gdata}/raw/Harris and Whiteway 2009/Global_Seascapes/class_11")
raw_phys_vents_csv       <- glue("{dir_gdata}/raw/Hydrothermal vents - Interridge Vent Database v3.4/vent_fields_all.csv")
abnj_shp                 <- glue("{dir_data}/abnj.shp")
abnj_s05_shp             <- glue("{dir_data}/abnj_s05.shp")
iho_shp                  <- glue("{dir_data}/iho.shp")
iho_s05_shp              <- glue("{dir_data}/iho_s05.shp")
ppow_shp                 <- glue("{dir_data}/ppow.shp")
ppow_s05_shp             <- glue("{dir_data}/ppow_s05.shp")
eez_shp                  <- glue("{dir_data}/eez.shp")
eez_s05_shp              <- glue("{dir_data}/eez_s05.shp")
pu_id_tif                <- glue("{dir_data}/pu_id.tif")
mine_claims_shp          <- glue("{dir_data}/mine-claims.shp")
mine_claims_tif          <- glue("{dir_data}/mine-claims.tif")

# helper functions ----
lyr_to_tif <- function(lyr, s, pfx){
  dir_tif <- here(glue("data-raw/{pfx}"))
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

# p_abnj ----
if (!file.exists(abnj_s05_shp)){
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
  #write_sf(p_iho, iho_shp)
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

if (!file.exists(abnj_ppow_s05_shp)){
  ppow <- read_sf(raw_ppow_shp)
  abnj <- read_sf(abnj_shp)

  abnj_ppow <- st_intersection(ppow, abnj)
  abnj_ppow <- abnj_ppow %>%
    select(-FID, -one) %>%
    mutate(
      area_km2 = st_area(geometry) %>%  units::set_units(km^2) %>% as.integer())
  #p_abnj_ppow <- abnj_ppow
  write_sf(p_abnj_ppow, abnj_ppow_shp)
  p_abnj_ppow <- read_sf(abnj_ppow_shp)
  use_data(p_abnj_ppow, overwrite = TRUE)

  p_abnj_ppow_s05 <- rmapshaper::ms_simplify(p_abnj_ppow, keep = 0.05)
  use_data(p_abnj_ppow_s05, overwrite = TRUE)
  write_sf(p_abnj_ppow_s05, abnj_ppow_s05_shp)
}

# r_na ----
# for assigning fresh values later
r_na <- raster(
  xmn = -180, xmx = 180, ymn = -90, ymx = 90,
  resolution=0.5, crs=leaflet:::epsg4326, vals=NA)

# r_pu_id ----
if (!file.exists(pu_id_tif)){
  r_abnj <- rasterize(abnj, r_na)
  r_pu_id <- r_na
  values(r_pu_id) <- 1:ncell(r_abnj)
  r_pu_id <- mask(r_pu_id, r_abnj)

  writeRaster(r_pu_id, pu_id_tif, overwrite=T)
  use_data(r_pu_id * 1, overwrite = TRUE)
}

# testing ----
# library(tidyverse)
# library(raster)
# raster("inst/data/bio_gmbi/nspp_all.tif") %>%
#   plot()

# library(sf)
# load("data/r_phys_seamounts.rda")
# r_phys_seamounts
# plot(r_mine_claim)

# use_data()
# load("data/p_abnj_ppow_s05.rda")
# p_abnj_ppow_s05
# plot(p_abnj_ppow_s05['ECOREGION'])
# class(p_abnj_ppow_s05)

# p_eez ----
if (!file.exists(eez_shp)){
  p_eez <- read_sf(raw_eez_shp)
  write_sf(p_eez, here("data-raw/eez.shp"))
  use_data(p_eez, overwrite = TRUE)

  p_eez_s05 <- rmapshaper::ms_simplify(p_abnj, keep = 0.05)
  use_data(p_eez_s05, overwrite = TRUE)
  write_sf(p_eez_s05, eez_s05_shp)
}

# p_highseas ----
p_highseas <- read_sf(highseas_shp)
write_sf(p_highseas, here("data-raw/highseas.shp"))
use_data(p_highseas, overwrite = TRUE)

# s_bio ----

# get bio raster stack from gmbi
data(gmbi_indicators)

# mask stack
s_bio_gmbi <- gmbi_indicators %>%
  mask(r_pu_id)

# write tifs
lyrs <- names(s_bio_gmbi)
map(lyrs, lyr_to_tif, s_bio_gmbi, "bio_gmbi")

# use_data()
use_data(s_bio_gmbi, overwrite = TRUE)

# s_fish_gfw ----

# read csv, identify cell_id
d_fish_gfw <- read_csv(fish_gfw_csv) %>%
  mutate(
    cell_id = cellFromXY(r_pu_id, matrix(data = c(lon_bin_center, lat_bin_center), ncol=2)))

# create raster stack
if (exists("s_fish_gfw")) rm(s_fish_gfw)
lyrs <- names(d_fish_gfw)[4:(ncol(d_fish_gfw) - 1)]
for (lyr in lyrs){ # lyr = lyrs[2]
  r <- r_na
  names(r) <- lyr
  r[d_fish_gfw$cell_id] <- d_fish_gfw[[lyr]]

  if (!exists("s_fish_gfw")){
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

# s_fish_saup ----
fish_saup_now_tif           <- glue("{dir_gdata}/raw/UBC-exploited-fish-projections/Current_MCP1.tif")
fish_saup_future_tif        <- glue("{dir_gdata}/raw/UBC-exploited-fish-projections/MCP2050_RCP85.tif")

s_fish_saup  <- stack(c(fish_saup_now_tif, fish_saup_future_tif))
lyrs <- names(s_fish_saup) <- c("mcp_2004", "mcp_2050")

# mask stack
s_fish_saup <- mask(s_fish_saup, r_pu_id)

# write tifs
map(lyrs, lyr_to_tif, s_fish_saup, "fish_saup")

# use_data()
use_data(s_fish_saup, overwrite = TRUE)


# r_phys_seamounts ----
x <- read_xml(phys_seamounts_kml)

xpaths <- list(
  name = "/d1:kml/d1:Document/d1:Folder/d1:Folder/d1:name",
  xyz  = "/d1:kml/d1:Document/d1:Folder/d1:Folder/d1:Placemark/d1:Point/d1:coordinates")
pts <- tibble(
  names = xml_find_all(x, xpaths$name) %>% xml_text(),
  xyz   = xml_find_all(x, xpaths$xyz) %>% xml_text()) %>%
  separate(xyz, c("x", "y", "z"), sep=",", convert=T) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

r_phys_seamounts <- rasterize(st_coordinates(pts), r_pu_id, fun='count', background=0) %>%
  mask(r_pu_id) # plot(r_phys_seamounts)
names(r_phys_seamounts) <- "count"

writeRaster(r_phys_seamounts, here("data-raw/phys_seamount_count.tif"), overwrite = TRUE)
use_data(r_phys_seamounts, overwrite = TRUE)

# s_phys_scapes ----

# load raster in 0.1 deg resolution
r_scapes <- raster(phys_scapes_arcinfo)

# create raster stack of benthicscape classes
# with each layer containing amount of class in 0.5 deg cell
if (exists("s_scapes")) rm(s_scapes)
r_a <- area(r_scapes)
for (i in 1:11){ # i = 2
  r_i <- (r_scapes == i) * r_a
  r_ia <- aggregate(r_i, fact=5, fun=sum, expand=F, na.rm=T)
  #plot(r_ia, col=cols)

  if (!exists("s_scapes")){
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

# r_phys_vents ----
pts_vents <- read_csv(phys_vents_csv) %>%
  st_as_sf(
    coords = c("Longitude", "Latitude"), crs = 4326, remove=F)

r_phys_vents <- rasterize(
  st_coordinates(pts_vents), r_pu_id, fun='count', background=0) %>%
  mask(r_pu_id) # plot(r_phys_seamounts)
names(r_phys_vents) <- "count"

writeRaster(r_phys_vents, here("data-raw/phys_vents_count.tif"), overwrite = TRUE)
use_data(r_phys_vents, overwrite = TRUE)

# r_mine_claims ----
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
