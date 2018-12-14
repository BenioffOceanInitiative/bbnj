source(here::here("inst/scripts/setup.R"))
devtools::load_all()

res_out <- 0.5
dir_gfw_out <- file.path(dir_data_gdrive, "derived/fishing")

gfw   <- read_csv(gfw_csv)
res_gfw <- 0.5
r_gfw <- get_grid(res=res_gfw)

# get cells
cells <- gfw %>%
  mutate(
    cell = cellFromXY(r_gfw, matrix(data = c(lon_bin_center, lat_bin_center), ncol=2)))
#cat(paste(names(cells), collapse='", "'))
# "year", "lat_bin_center", "lon_bin_center", "fishing_KWH", "mean_costs", "revenue", "mean_scaled_profits", "mean_scaled_profits_with_subsidies", "scaled_profits_low_labor_cost", "cell"

for (nm in names(cells)[4:(ncol(cells)-1)]){ # nm = names(cells)[4] # nm = names(cells)[ncol(cells)]
  tif <- glue("gfw_{nm}_{res_out}dd.tif")
  cat(glue("{tif}\n\n"))

  # assign to raster
  r <- r_gfw
  r[cells$cell] <- cells[[nm]]
  #r <- raster::trim(r) # skip so all aligned
  #plot(r)
  #plot(r_new)
  #devtools::load_all()
  r_new <- rescale_grid(r, res_out)

  # write out
  writeRaster(r, file.path(dir_gfw_out, tif), overwrite=T)
}

