source(here::here("inst/scripts/setup.R"))

# World Vector Shoreline
# https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/

# data
eez_sf  <- sf::read_sf(eez_shp)

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
plot(bb)

#Error in st_crs(x) == st_crs(y) :
#  Evaluation error: TopologyException: Input geom 1 is invalid: Ring Self-intersection at or near point -70.931989984449672 -54.755801231741259 at -70.931989984449672 -54.755801231741259.

sf_erase <- function(x, y){
  st_difference(x, st_union(st_combine(y))) } # 215.873

# x <- sf_erase(bb, eez_sf)
eez_b <- st_buffer(eez_sf, 0) # 396.303 sec

system.time(eez_bs   <- ms_simplify(eez_b, keep = 0.05)) # 569 sec
system.time(eez_bse  <- ms_erase(bb, eez_bs, remove_slivers = T)) # 11.5 sec
system.time(eez_bsee <- ms_explode(eez_bse)) # 8.9 sec

eez_bsee <- eez_bsee %>%
  mutate(area_km2 = st_area(geometry)/(1000*1000)) # 0.56 sec
# 195019807.766397536

write_sf(eez_bsee, file.path(dir_data_gdrive, "derived/eez_bsee.shp"))

eez_bsee <- read_sf(file.path(dir_data_gdrive, "derived/eez_bsee.shp"))
eez_bse <- read_sf(file.path(dir_data_gdrive, "derived/eez_bse.shp"))


library(rnaturalearth)

land110 <- ne_download(
  scale = 110, type = "land", category = c("physical"),
  returnclass = "sf",
  destdir = file.path(dir_data_gdrive, "raw/rnaturalearth"))
land_l <- ne_download(
  scale = "large", type = "land", category = c("physical"),
  returnclass = "sf",
  destdir = file.path(dir_data_gdrive, "raw/rnaturalearth"))
land10 <- ne_download(
  scale = 10, type = "land", category = c("physical"),
  returnclass = "sf",
  destdir = file.path(dir_data_gdrive, "raw/rnaturalearth"))
islands_l <- ne_download(
  scale = "large", type = "minor_islands", category = c("physical"),
  returnclass = "sf",
  destdir = file.path(dir_data_gdrive, "raw/rnaturalearth"))

islands_l <- islands_l %>%
  mutate(
    geometry = st_cast(geometry, "MULTIPOLYGON"))

names(land_l)
names(islands_l)
lands <- rbind(
  land_l %>% select(-scalerank),
  islands_l %>% select(-scalerank))

#system.time(lands <- st_union(land_l, islands_l)) #, by_feature = T)
system.time(lands_b100 <- st_buffer(lands, 1/60*100)) # 144.4 sec
#system.time(land110_b100 <- st_buffer(land110, 1/60*100)) # 0.8 sec
#system.time(land10_b100  <- st_buffer(land10, 1/60*100)) # 169.8 sec

#land110_b100 <- st_crop(land110_b100, bb)
lands_b100 <- st_crop(lands_b100, bb)

lands_b100 <- ms_dissolve(lands_b100)

# land10_b100 <- land10_b100 %>%
#   mutate(geometry = st_cast(geometry, "MULTIPOLYGON"))
# st_geometry_type(land110_b100$geometry) %>% table()

# land110_b100 <- land110_b100 %>%
#   mutate(geometry = st_make_valid(geometry)) # lwgeom::st_make_valid()

lands_b100_u <- st_sf(
  tibble(
    one = 1,
    geometry = st_union(st_buffer(lands_b100, 0.0)), crs = 4326))

#system.time(eez_bsel  <- ms_erase(eez_bse, land110_b100)) # 11.5 sec
#system.time(eez_bsel  <- ms_erase(eez_bse, land10_b100)) # 18.4 sec
#system.time(eez_bsel  <- ms_erase(eez_bse, lands_b100)) # 22.4 sec
system.time(eez_bsel  <- ms_erase(eez_bse, lands_b100_u)) # 6.1 sec

#land110_b100 <- st_crop(land110_b100, bb)


plot(eez_bsel)
lands_b100$one = 1
plot(lands_b100)
plot(lands_b100_u)

system.time(eez_bsele <- ms_explode(eez_bsel)) # 8.9 sec


bb_l <- st_difference(bb, st_combine(land110_b100))

eez_bsee$geometry

# https://www.r-spatial.org/r/2017/03/19/invalid.html
x <- land110_b100$geometry
bb <- bb %>% st_cast("MULTIPOLYGON")
x <- bb$geom
any(is.na(st_dimension(x)))
any(is.na(st_is_valid(x)))
any(na.omit(st_is_valid(x)) == FALSE)
st_is_valid(x, reason = TRUE)

st_cast()

system.time(eez_bsd <- ms_dissolve(eez_bs))
system.time(eez_bse <- ms_erase(bb, eez_bsd, remove_slivers = T))


system.time(eez_bcu <- eez_b %>% st_combine() %>% st_union())
system.time(x <- st_difference(bb, eez_b)) # 215 sec



system.time(x <- st_difference(bb, eez_b)) # 215 sec


# Dissolve in QGIS

write_sf(x, file.path(dir_data_gdrive, "derived/x.shp"))

system.time(x <- st_difference(bb, eez_b)) # 215 sec
system.time(
  x <- x %>%
    mutate(area_km2 = st_area()))
plot(x['GeoName'])
write_sf(x, here("data/tmp_x.shp"))  # qgis
earth = land %>%
  filter(
    st_intersects(land, bb_ply, sparse=F)[,1]) # plot(earth)

globe <- st_bbox() %>% st_as_sfc()
st_bbox(obj, ..., crs = NA_crs_)

raster::extent(-180, 180,-90,90) %>% sf::st_as_sf()


st_as_sfc(x, ..., crs = NA_crs_)


# old ...
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
