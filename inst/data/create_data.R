# 1) use_data(). for lazy loading with library(bbnj). document datasets in R/data.R.
# 2) write tifs. for use with other software external to R, eg QGIS

library(xml2)
library(purrr)
library(raster)
library(usethis)
library(sf)
library(glue)
library(here)
library(rmapshaper)

library(gmbi) #devtools::install_github("marinebon/gmbi") #devtools::install_local("~/github/gmbi")
devtools::load_all()

# paths ----
dir_gdata           <- "~/Gdrive Ecoquants/projects/bbnj/data" # on Ben Best's laptop
cell_res            <- 0.5
eez_shp             <- "/Users/bbest/github/sdg14-shiny/info-gl/data/eez_s005005.shp"
highseas_shp        <- glue("{dir_gdata}/derived/boundary/high_seas.shp")
pu_id_tif           <- glue("{dir_gdata}/derived/boundary/high_seas_cellid_{cell_res}dd.tif")
fish_gfw_csv        <- glue("{dir_gdata}/raw/Sala et al 2018/half_degree_binned_results.csv")
phys_seamounts_kml  <- glue("{dir_gdata}/raw/Seamounts - Kim and Wessel 2011/KWSMTSv01.kml")
phys_scapes_arcinfo <- glue("{dir_gdata}/raw/Harris and Whiteway 2009/Global_Seascapes/class_11")
phys_vents_csv      <- glue("{dir_gdata}/raw/Hydrothermal vents - Interridge Vent Database v3.4/vent_fields_all.csv")
fish_ubc_now_tif    <- glue("{dir_gdata}/raw/UBC-exploited-fish-projections/Current_MCP1.tif")
fish_ubc_future_tif <- glue("{dir_gdata}/raw/UBC-exploited-fish-projections/MCP2050_RCP85.tif")
mine_claims_shp     <- glue("{dir_gdata}/raw/ISA_shapefiles/ISA_claim_areas_update_20181202.shp")

# helper functions ----
lyr_to_tif <- function(lyr, s, pfx){
  dir_tif <- here(glue("data-raw/{pfx}"))
  if (!dir.exists(dir_tif)) dir.create(dir_tif)
  tif <- glue("{dir_tif}/{lyr}.tif")
  r <- raster(s, lyr) %>%
    mask(r_pu_id)
  writeRaster(r, tif, overwrite=T)
}

# r_pu_id ----
r_pu_id <- raster(pu_id_tif) # plot(r_pu_id)
writeRaster(r_pu_id, here("data-raw/pu_id.tif"), overwrite = TRUE)
use_data(r_pu_id, overwrite = TRUE)

# r_na ----
# for assigning fresh values later
r_na <- r_pu_id
values(r_na) <- NA

# p_eez ----
p_eez <- read_sf(eez_shp)
write_sf(p_eez, here("data-raw/eez.shp"))
use_data(p_eez, overwrite = TRUE)

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

# s_fish_ubc ----
fish_ubc_now_tif           <- glue("{dir_gdata}/raw/UBC-exploited-fish-projections/Current_MCP1.tif")
fish_ubc_future_tif        <- glue("{dir_gdata}/raw/UBC-exploited-fish-projections/MCP2050_RCP85.tif")

s_fish_ubc  <- stack(c(fish_ubc_now_tif, fish_ubc_future_tif))
lyrs <- names(s_fish_ubc) <- c("mcp_2004", "mcp_2050")

# mask stack
s_fish_ubc <- mask(s_fish_ubc, r_pu_id)

# write tifs
map(lyrs, lyr_to_tif, s_fish_ubc, "fish_ubc")

# use_data()
use_data(s_fish_ubc, overwrite = TRUE)


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

# r_mine_claims
p_mine_claims <- read_sf(mine_claims_shp)
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
write_sf(p_mine_claims, here("data-raw/mine_claims.shp"))
use_data(p_mine_claims, overwrite = TRUE)

r_mine_claim <- rasterize(p_mine_claims, r_pu_id, field=1, background=0) %>%
  mask(r_pu_id) # plot(r_mine_claims)
names(r_mine_claim) <- "present"
writeRaster(r_mine_claim, here("data-raw/mine_claim_present.tif"), overwrite = TRUE)
use_data(r_mine_claim, overwrite = TRUE)

