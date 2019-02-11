library(shiny)
library(tidyverse)
library(leaflet)
library(glue)
library(shinydashboard)
library(plotly)
library(sf)
library(raster)

# paths
dir_gdata <- "~/Gdrive Ecoquants/projects/bbnj/data"
abnj_shp  <- file.path(dir_gdata, "derived/boundary/high_seas_s05.shp")
dir_gfw   <- file.path(dir_gdata, "derived/fishing")

abnj <- read_sf(abnj_shp) %>%
  st_buffer(dist=0)

tbl_gfw <- tibble(
  tif = list.files(dir_gfw, ".*\\.tif$", full.names = T),
  name = str_replace(tif, "(.*)/gfw_(.*).tif$", "\\2"))


# https://github.com/SFG-UCSB/The-economics-of-fishing-the-high-seas
raster("/Users/bbest/Gdrive Ecoquants/projects/bbnj/data/derived/fishing/gfw_revenue_0.1dd.tif")
raster("/Users/bbest/Gdrive Ecoquants/projects/bbnj/data/derived/fishing/gfw_revenue_0.5dd.tif")
