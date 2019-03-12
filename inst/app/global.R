library(raster)
library(tidyverse)
library(stringr)
select = dplyr::select
library(sf)
library(leaflet)
library(mapview)
library(here)
library(glue)
library(shiny)
library(shinydashboard)
#library(plotly)
library(bbnj) # devtools::load_all() # devtools::install_github("ecoquants/bbnj")

lyrs_grd <- here("inst/app/lyrs.grd")

# 0. eez
data(p_eez)

if (!file.exists(lyrs_grd)){

  # 1. Features, Original
  features_original <- stack(
    raster(s_bio_gmbi, "nspp_all"),
    raster(s_bio_gmbi, "rls_all"),
    raster(s_fish_gfw, "mean_scaled_profits_with_subsidies"),
    s_fish_ubc,
    r_phys_seamounts,
    r_phys_vents,
    r_mine_claim,
    s_phys_scapes)
  names(features_original) <- c(
    "bio_nspp",
    "bio_rls",
    "fish_profit.subs",
    "fish_mcp.2004",
    "fish_mcp.2050",
    "phys_seamounts",
    "phys_vents",
    "mine_claim",
    glue("phys_scape.{1:11}"))
  names(features_original) <- glue("Features.original_{names(features_original)}")

  # 2. Features, Rescaled
  features_rescaled <- stack(
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
    r_mine_claim,
    rescale_stack(s_phys_scapes))
  names(features_rescaled) <- c(
    "bio_nspp",
    "bio_rls",
    "fish_profit.subs",
    "fish_mcp.2004",
    "fish_mcp.2050",
    "phys_seamounts",
    "phys_vents",
    "mine_claim",
    glue("phys_scape.{1:11}"))
  names(features_rescaled) <- glue("Features.rescaled_{names(features_rescaled)}")

  lyrs <- stack(
    features_original,
    features_rescaled)

  lyrs <- projectRasterForLeaflet(lyrs, "ngb")

  writeRaster(lyrs, lyrs_grd)
} else {
  lyrs <- stack(lyrs_grd)
}

lyr_choices <- list(
  `Features, original` =
    setNames(
      names(lyrs) %>%
        str_subset("^Features.original.*"),
      names(lyrs) %>%
        str_subset("^Features.original.*") %>%
        str_replace("Features.original_", "")),
  `Features, rescaled` =
    setNames(
      names(lyrs) %>%
        str_subset("^Features.rescaled.*"),
      names(lyrs) %>%
        str_subset("^Features.rescaled.*") %>%
        str_replace("Features.rescaled_", "")))


# 2. Features, 10%
# 3. Features, 30%
# 4. Scenarios
