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
library(bbnj)
# devtools::load_all()
# devtools::install_github("ecoquants/bbnj", force=T)
# devtools::install_local(force=T)
#library(prioritizr)

# 0. eez
data(p_eez_s05)
data(p_abnj_s05)
data(p_ppow_s05)

lyrs_rda <- file.path(system.file(package="bbnj"), "app/lyrs.rda")

if (!file.exists(lyrs_rda)){

  # 1. Features, Original
  features_original <- stack(
    raster(s_bio_gmbi, "nspp_all"),
    raster(s_bio_gmbi, "rls_all"),
    raster(s_fish_gfw, "mean_scaled_profits_with_subsidies"),
    s_fish_ubc,
    r_phys_seamounts,
    r_phys_vents,
    r_mine_claims,
    s_phys_scapes)
  names(features_original) <- c(
    "bio_nspp",
    "bio_rls",
    "fish_profit.subs",
    "fish_mcp.2004",
    "fish_mcp.2050",
    "phys_seamounts",
    "phys_vents",
    "mine_claims",
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
    r_mine_claims,
    rescale_stack(s_phys_scapes))
  names(features_rescaled) <- c(
    "bio_nspp",
    "bio_rls",
    "fish_profit.subs",
    "fish_mcp.2004",
    "fish_mcp.2050",
    "phys_seamounts",
    "phys_vents",
    "mine_claims",
    glue("phys_scape.{1:11}"))
  names(features_rescaled) <- glue("Features.rescaled_{names(features_rescaled)}")

  # 3. Taxa, nspp
  lyrs_nspp <- names(s_bio_gmbi) %>% str_subset("^nspp_(?!all)")

  taxa_nspp <- subset(s_bio_gmbi, lyrs_nspp)
  names(taxa_nspp) <- str_replace(lyrs_nspp, "nspp_", "")
  names(taxa_nspp) <- glue("Taxa.nspp_{names(taxa_nspp)}")

  # 4. Taxa, rls
  lyrs_rls <- names(s_bio_gmbi) %>% str_subset("^rls_(?!all)")

  taxa_rls <- subset(s_bio_gmbi, lyrs_rls)
  names(taxa_rls) <- str_replace(lyrs_rls, "rls_", "")
  names(taxa_rls) <- glue("Taxa.rls_{names(taxa_rls)}")

  lyrs_geo <- stack(
    features_original,
    features_rescaled,
    taxa_nspp,
    taxa_rls)

  lyrs_mer <- projectRasterForLeaflet(lyrs_geo, "ngb")

  lyrs_geo <- readAll(lyrs_geo)
  lyrs_mer <- readAll(lyrs_mer)

  #writeRaster(lyrs, lyrs_grd)
  save(lyrs_geo, lyrs_mer, file=lyrs_rda, compress="bzip2")
} else {
  #lyrs <- stack(lyrs_grd)
  load(lyrs_rda)
  #lyrs_geo <- readAll(lyrs_geo)
  #lyrs_mer <- readAll(lyrs_mer)
}

lyr_choices <- list(
  `Features, original` =
    setNames(
      names(lyrs_mer) %>%
        str_subset("^Features.original.*"),
      names(lyrs_mer) %>%
        str_subset("^Features.original.*") %>%
        str_replace("Features.original_", "")),
  `Features, rescaled` =
    setNames(
      names(lyrs_mer) %>%
        str_subset("^Features.rescaled.*"),
      names(lyrs_mer) %>%
        str_subset("^Features.rescaled.*") %>%
        str_replace("Features.rescaled_", "")),
  `Taxa, nspp` =
    setNames(
      names(lyrs_mer) %>%
        str_subset("^Taxa.nspp.*"),
      names(lyrs_mer) %>%
        str_subset("^Taxa.nspp.*") %>%
        str_replace("Taxa.nspp_", ""))
  ,
  `Taxa, rls` =
    setNames(
      names(lyrs_mer) %>%
        str_subset("^Taxa.rls.*"),
      names(lyrs_mer) %>%
        str_subset("^Taxa.rls.*") %>%
        str_replace("Taxa.rls_", "")))

# 2. Features, 10%
# 3. Features, 30%
# 4. Scenarios
