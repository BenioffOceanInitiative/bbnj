library(sf)
library(rmapshaper)

# Caroline prepared this
dir_b   <- "~/Google Drive/projects/bbnj/data/derived/boundary"
b_shp   <- file.path(dir_b, "high_seas.shp")
b_s_shp <- file.path(dir_b, "high_seas_s05.shp")
#b_s_shp <- file.path(dir_b, "high_seas_final.shp")

# TODO: redo high_seas_boundary.shp

if (!file.exists(b_s_shp)){
  b <- read_sf(b_shp)
  b_s <- ms_simplify(b, keep = 0.05)
  write_sf(b_s, b_s_shp)
}
