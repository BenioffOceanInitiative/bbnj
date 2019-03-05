# 1) use_data(). for lazy loading with library(bbnj). document datasets in R/data.R.
# 2) write tifs. for use with other software external to R, eg QGIS

#devtools::install_github("marinebon/gmbi")
devtools::load_all()
library(raster)
library(usethis)

# paths ----
dir_gdata      <- "~/Gdrive Ecoquants/projects/bbnj/data"
cell_res       <- 0.5  # n=  113,401
p_boundary_shp <- glue("{dir_gdata}/derived/boundary/high_seas.shp")
pu_id_tif      <- glue("{dir_gdata}/derived/boundary/high_seas_cellid_{cell_res}dd.tif")
fish_gfw_csv   <- glue("{dir_gdata}/raw/Sala et al 2018/half_degree_binned_results.csv")

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
r_pu_id <- raster(pu_id_tif)
writeRaster(r_pu_id, here("data-raw/pu_id.tif"), overwrite = TRUE)
use_data(r_pu_id, overwrite = TRUE)

# s_bio ----

# get bio raster stack from gmbi
library(gmbi)
data(gmbi_indicators)

# mask stack
s_bio <- gmbi_indicators %>%
  mask(r_pu_id)

# write tifs
lyrs <- names(s_bio)
map(lyrs, lyr_to_tif, s_bio, "bio")

# use_data()
use_data(s_bio, overwrite = TRUE)

# s_fish_gfw ----

# read csv, identify cell_id
d_fish_gfw <- read_csv(fish_gfw_csv) %>%
  mutate(
    cell_id = cellFromXY(r_pu_id, matrix(data = c(lon_bin_center, lat_bin_center), ncol=2)))

# create raster stack
if (exists("s_fish_gfw")) rm(s_fish_gfw)
lyrs <- names(d_fish_gfw)[4:(ncol(d_fish_gfw) - 1)]
for (lyr in lyrs){ # lyr = lyrs[2]
  r <- r_pu_id
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




