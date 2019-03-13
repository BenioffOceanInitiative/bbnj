p_feature_names <- names(lyrs) %>%
  str_subset("^Features.rescaled.*")
p_features <- subset(lyrs, p_feature_names)
names(p_features) <- str_replace(p_feature_names, "Features.rescaled_", "")

#devtools::install_local(force=T)
dir_dat          <- system.file("data-raw", package="bbnj")
r_pu_id          <- raster(file.path(dir_dat, "pu_id.tif"))
s_bio_gmbi       <- stack(list.files(file.path(dir_dat, "bio_gmbi"), "\\.tif$", full.names=T))
s_fish_gfw       <- stack(list.files(file.path(dir_dat, "fish_gfw"), "\\.tif$", full.names=T))
s_fish_ubc       <- stack(list.files(file.path(dir_dat, "fish_ubc"), "\\.tif$", full.names=T))
r_phys_seamounts <- raster(file.path(dir_dat, "phys_seamounts.tif"))
r_phys_vents     <- raster(file.path(dir_dat, "phys_vents.tif"))
s_phys_scapes    <- stack(list.files(file.path(dir_dat, "phys_scapes"), "\\.tif$", full.names=T))
r_mine_claims    <- raster(file.path(dir_dat, "mine_claims.tif"))

p_features <- stack(
  raster(s_bio_gmbi, "nspp_all") %>%
    rescale_raster(multiply_area=T),
  raster(s_bio_gmbi, "rls_all") %>%
    rescale_raster(multiply_area=T),
  raster(s_fish_gfw, "mean_scaled_profits_with_subsidies") %>%
    gap_fill_raster() %>%
    rescale_raster(inverse=T),
  rescale_stack(s_fish_ubc, inverse=T),
  rescale_raster(r_phys_seamounts),
  rescale_raster(r_phys_vents),
  rescale_stack(s_phys_scapes))
names(p_features) <- c(
  "bio_nspp",
  "bio_rls",
  "fish_profit.subs",
  "fish_mcp.2004",
  "fish_mcp.2050",
  "phys_seamounts",
  "phys_vents",
  glue("phys_scape.{1:11}"))

r_pu_area <- area(r_pu_id) %>%  # in km2
  mask(r_pu_id)
A <- cellStats(r_pu_area, "sum")
r_pu_areas <- r_pu_area / A

# sudo apt-get install coinor-libsymphony-dev coinor-libcgl-dev autotools-dev
# install.packages("Rsymphony")
#install.packages("gurobi")

p <- problem(r_pu_areas, p_features) %>%
  add_max_utility_objective(budget = 0.1) #%>%  # 10% of total high seas area
#add_gurobi_solver()

p_sol <- solve_log(p, redo=F)
# p01_sol_Rsymphony <- p_sol
#writeRaster(p_sol, "inst/app/p01_sol_Rsymphony.tif")
plot(p_sol)

p01_sol <- raster("vignettes/p01_sol.tif")
plot(p01_sol)

plot(p_sol)
plot(p01_sol, col = c("grey90", "darkgreen"), main = "p01 solution")

# area of solution
cellStats(r_pu_areas * p01_sol, "sum")

# calculate how well features are represented in the solution
feature_representation(p01, p01_sol)
