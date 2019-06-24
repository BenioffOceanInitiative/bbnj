library(rgdal)
library(raster)
library(tidyverse)
library(stringr)
library(sf)
library(leaflet)
library(mapview)
library(leafem)
library(here)
library(glue)
library(shiny)
library(htmltools)
library(shinydashboard)
library(bbnj)
# devtools::load_all() # devtools::install_local(force=T) # devtools::install_github("ecoquants/bbnj")
library(bsplus)
select = dplyr::select

# rsconnect::showLogs() # show log on shinyapps.io to debug

# setwd if working from github repo vs shinyapps.io; reset: setwd(here())
if (!"global.R" %in% list.files(getwd())) setwd(here("inst/app"))

# paths
input_features_csv <- "data/input_features.csv"
id_gcs2mer_csv     <- "data/id_gcs2mer.csv"
na_mer_tif         <- "data/na_mer.tif"
dir_scenarios      <- "www/scenarios"

# input features ----
features        <- read_csv(input_features_csv) %>%
  filter(active)
features_label <- features %>%
  mutate(
    label = map(r_labels, function(x) eval(parse(text=x), envir=globalenv()))) %>%
  select(type, label) %>%
  unnest(label) # View(lyrs_labels)
features_stack    <- map(features$r_data, function(x) eval(parse(text=x))) %>% stack()
names(features_stack) <- features_label$label

# output scenarios ----
scenarios_label <- tibble(
  type = "output_scenario",
  tif  = list.files(dir_scenarios, "^s.*\\_sol.tif$", full.names=T)) %>%
  mutate(
    label     = map_chr(basename(tif), function(x) str_replace(x, "^(s.*)\\_sol.tif$", "\\1")))
scenarios_stack        <- stack(scenarios_label$tif)
names(scenarios_stack) <- scenarios_label$label

# layers: features + scenarios ----
lyrs_gcs <- stack(features_stack, scenarios_stack)

# project raster stack from gcs to mer for leaflet ----
if (any(!file.exists(na_mer_tif, id_gcs2mer_csv))){
  r_pu_id_mer      <- projectRasterForLeaflet(r_pu_id, "ngb")
  r_na_mer         <- r_pu_id_mer
  values(r_pu_id_mer) <-  NA
  writeRaster(r_pu_id_mer, na_mer_tif)

  id <- tibble(
    mer = 1:ncell(r_pu_id_mer),
    gcs = values(r_pu_id_mer)) %>%
    filter(!is.na(gcs))
  write_csv(id, id_gcs2mer_csv)
} else {
  id       <- read_csv(id_gcs2mer_csv)
  r_na_mer <- raster(na_mer_tif)
}

lyrs_mer <- stack(lapply(1:nlayers(lyrs_gcs), function(x) r_na_mer))
for (i in 1:nlayers(lyrs_gcs)){ # i = 1
  lyrs_mer[[i]][id$mer] <- lyrs_gcs[[i]][id$gcs]
}
names(lyrs_mer) <- names(lyrs_gcs)

# setup ui layer choices ----
lyr_choices <- bind_rows(
    features_label,
    scenarios_label) %>%
  select(type, label) %>%
  arrange(type, label)

lyrs_mismatch <- setdiff(names(lyrs_mer), lyr_choices$label)
stopifnot(length(lyrs_mismatch) == 0)
