# libraries ----
if (!require("pacman")) install.packages("pacman")
library(pacman)   # pacman::p_load() installs package if not present vs just library()
p_load(tidyverse) # ggplot2,tibble,tidyr,readr,purrr,dplyr,stringr,forcats
p_load(glue, here, fs)                              # text, paths
p_load(rgdal, ncdf4, raster, leaflet)               # spatial
p_load(sf, rmapshaper, lwgeom)                       # vector
p_load(fasterize) # raster to vector; p_load_gh("ecohealthalliance/fasterize")
p_load(knitr, rmarkdown, htmltools, DT, htmltools)  # reporting
p_load(RColorBrewer)                                # graphics
#p_update()
select <- dplyr::select # vs raster::select

# paths ----
dir_data_gdrive <- "/Users/bbest/Google Drive/projects/Pew BBNJ/data"
elev_nc         <- file.path(dir_data_gdrive, "raw/gebco.net_depth/GEBCO_2014_2D.nc")
eez_land_shp    <- file.path(dir_data_gdrive, "raw/marineregions.org_boundaries/EEZ_land_union_v2_201410/EEZ_land_v2_201410.shp")
eez_shp         <- file.path(dir_data_gdrive, "raw/marineregions.org_boundaries/World_EEZ_v10_20180221/eez_v10.shp")
gfw_csv         <- file.path(dir_data_gdrive, "raw/Sala et al 2018/half_degree_binned_results.csv")

