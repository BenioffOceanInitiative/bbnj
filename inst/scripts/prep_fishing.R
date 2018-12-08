source(here::here("inst/scripts/setup.R"))


#sp_id_to_tif <- function(sp_id){
  # , cells=cells, r_na=r_na, dir_spp=dir_spp

  #devtools::load_all()
  #r <- get_fishing_empty_grid()


dir_gfw_out <- file.path(dir_data_gdrive, "derived/fishing")

  r_na <- raster(
    xmn = -180, xmx = 180, ymn = -90, ymx = 90,
    resolution=0.5, crs=leaflet:::epsg4326)

  gfw <- read_csv(gfw_csv)

  # get cells
  cells <- gfw %>%
    mutate(
      cell = cellFromXY(r_na, matrix(data = c(lon_bin_center, lat_bin_center), ncol=2)))
  cat(paste(names(cells), collapse='", "'))
  # "year", "lat_bin_center", "lon_bin_center", "fishing_KWH", "mean_costs", "revenue", "mean_scaled_profits", "mean_scaled_profits_with_subsidies", "scaled_profits_low_labor_cost", "cell"

  for (nm in names(cells)[4:(ncol(cells)-1)]){
    tif <- glue("gfw_{nm}.tif")
    cat(glue("{tif}\n\n"))

    # assign to raster
    r <- r_na
    r[cells$cell] <- cells[[nm]]
    r <- raster::trim(r) # plot(r)

    # write out
    writeRaster(r, file.path(dir_gfw_out, tif), overwrite=T)
  }
#}
