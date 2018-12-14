library(shiny)
library(tidyverse)
library(leaflet)
library(glue)
library(shinydashboard)
library(plotly)
library(sf)
library(raster)

# paths
dir_gdata <- "~/Google Drive/projects/Pew BBNJ/data"
abnj_shp  <- file.path(dir_gdata, "derived/Caroline - high seas layer/high_seas_final.shp")
dir_gfw   <- file.path(dir_gdata, "derived/fishing")

abnj <- read_sf(abnj_shp) %>%
  st_buffer(dist=0)

tbl_gfw <- tibble(
  tif = list.files(dir_gfw, ".*\\.tif$", full.names = T),
  name = str_replace(tif, "(.*)/gfw_(.*).tif$", "\\2"))


