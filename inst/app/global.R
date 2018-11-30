library(shiny)
library(tidyverse)
library(leaflet)
library(glue)
library(shinydashboard)
library(plotly)
library(sf)

# paths
dir_gdata <- "~/Google Drive"
abnj_shp <- file.path(dir_gdata, "projects/Pew BBNJ/data/derived/Caroline - high seas layer/high_seas_final.shp")

abnj <- read_sf(abnj_shp) %>%
  st_buffer(dist=0)

