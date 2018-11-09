source(here::here("inst/scripts/setup.R"))

eez_sf      <- sf::read_sf(eez_shp)
eez_land_sf <- sf::read_sf(eez_land_shp)
depth_r     <- raster::raster(elev_nc, layer = "elevation") * -1
depth_r[depth_r < 0] <- NA

depth_r <- depth
# coral triangle
depth_r_ct <- raster::crop(depth_r, extent(97, 130, -2, 13))


arcgrd <- file.path(dir_data_gdrive, "data/raw/Harris and Whiteway 2009/Global_Seascapes/class_11")
seascapes <- raster(arcgrd)

plot(depth_r_ct)
plot()

eez_land_sf <- eez_land_sf %>%
  mutate(one = 1)
eez_land_r <- fasterize(eez_land_1, depth, field=one)

# QGIS: menu Vector > Geoprocessing Tools > Dissolve
# QGIS: menu Vector > Geoprocessing Tools > Dissolve
eez_land_d <- ms_dissolve(eez_land)
eez_s      <- ms_simplify(eez, keep = 0.05)
eez_sd     <- ms_dissolve(eez_s)

fasterize(eez_land)
# Error in CPL_geos_union(st_geometry(x), by_feature) :
#   Evaluation error: TopologyException: Input geom 1 is invalid: Ring Self-intersection at or near point -70.931989984449672 -54.755801231741259 at -70.931989984449672 -54.755801231741259.

eez_v <- eez %>%
  lwgeom::st_make_valid()

st_union(x)

eez_land_v <- eez_land %>%
  lwgeom::st_make_valid()

eez_d <- eez_v %>%
  mutate(
    one = 1) %>%
  group_by(one) %>%
  summarise() %>%
  st_cast() %>%
  mapview(zcol = "m")
