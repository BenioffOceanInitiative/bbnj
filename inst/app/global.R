library(rgdal)
library(raster)
library(readr)
library(purrr)
library(dplyr)
library(tidyr)
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
#devtools::load_all() # devtools::install_local(force=T) # devtools::install_github("ecoquants/bbnj")
library(bsplus)
select = dplyr::select

# rsconnect::showLogs() # show log on shinyapps.io to debug

# setwd if working from github repo vs shinyapps.io; reset: setwd(here())
if (!"global.R" %in% list.files(getwd())) setwd(here("inst/app"))

# variables
#prjres            <- "" # see projections_tbl$prjres
redo_features_rds <- T
debug <- T

# paths
features_csv <- "data/input_features.csv"
features_rds <- "data/features.rds"
#feature_layers_csv  <- "data/input_feature_layers.csv"
#id_gcs2mer_csv   <- "data/id_gcs2mer.csv"
#na_mer_tif       <- "data/na_mer.tif"
#id_gcs2mol_csv   <- "data/id_gcs2mol.csv"
#na_mol_tif       <- "data/na_mol.tif"
dir_scenarios <- "www/scenarios"

# input features ----
features <- read_csv(features_csv) %>%
  filter(active) %>%
  mutate(
    #prjres = ifelse(is.na(prjres), "", prjres))
    prjres = "_mer36km")
feature_labels <- features %>%
  mutate(
    label = map(r_labels, function(x) eval(parse(text=x), envir=globalenv()))) %>%
  select(type, label, prjres) %>%
  unnest(label)

# redo features_rds if features_csv modified more recently
mod_features_rds <- features_csv %>%
  fs::file_info() %>%
  pull(modification_time) >
  features_rds %>%
  fs::file_info() %>%
  pull(modification_time)
#redo_features_rds <- T
if (redo_features_rds | mod_features_rds){
  #for (prjres in unique(features$prjres)){ # prjres = "_mer36km"
  # TODO: update this to reflect that only geographic projection mappable

  # raster(s_fish_gfw, "mean_scaled_profits_with_subsidies") # "fish_profit.subs"
  r_pu_id <- get_d_prjres("r_pu_id", "")
  s_features <- features %>%
    #filter(prjres == !!prjres) %>%
    mutate(
      raster = map(r_data, function(x) eval(parse(text=x)))) %>%
    pull(raster) %>%
    stack()

  s_lbls <- feature_labels %>%
    #filter(prjres == !!prjres) %>%
    pull(label)

  nlyrs <- nlayers(s_features)
  nlbls <- length(s_lbls)
  stopifnot(nlyrs == nlbls)

  prjres <- "_mer36km"
  P <- projections_tbl %>% filter(prjres == !!prjres)
  r_pu_id_pr <- get_d_prjres("r_pu_id", prjres)
  s_features <- suppressWarnings(
    raster::projectRaster(
      s_features, raster::projectExtent(s_features, crs = sp::CRS(P$proj)),
      res = P$res_num)) %>%
    mask(r_pu_id_pr)

  names(s_features) <- s_lbls
  s_features <- raster::readAll(s_features)
  saveRDS(s_features, features_rds)
  #}
}

# read in features
# for (prjres in unique(features$prjres)){ # prjres = ""
#   assign(
#     glue("s_features{prjres}"),
#     readRDS(glue("data/features{prjres}.rds")))
# }
s_features <-  readRDS(features_rds)

# %>%
  #write_csv(feature_layers_csv) # features
 # View(lyrs_labels)
  #features_stack    <- map(features$r_data, function(x) eval(parse(text=x))) %>% stack()
  #names(features_stack) <- features_label$label
#}

# output scenarios ----
scenarios <- tibble(
  type = "output_scenario",
  tif  = list.files(dir_scenarios, "^s.*\\_sol.tif$", full.names=T)) %>%
  mutate(
    label   = map_chr(basename(tif), function(x) str_replace(x, "^(s.*)\\_sol.tif$", "\\1")),
    prjres  = map_chr(tif, function(x) get_tif_projection(x, debug=F)$prjres),
    gcs_shp = map_chr(tif, function(x) glue("{fs::path_ext_remove(x)}_gcs.shp")),
    has_shp = map_lgl(gcs_shp, file.exists))
# TODO: show scenario overview table; get total & percent area

stopifnot(all(scenarios$has_shp))
# scenarios <- scenarios %>%
#   mutate(
#     has_shp = map_lgl(gcs_shp, file.exists),
#     tif_to_shp_gcs())
# scenarios %>%
#   filter(!has_shp) %>%
#   pull(tif) %>%
#   walk(tif_to_shp_gcs)

# scenarios_stack        <- stack(scenarios_label$tif)
# names(scenarios_stack) <- scenarios_label$label
# for (prjres in unique(scenarios$prjres)){ # prjres = "_mol50km"
#   v_pr <- glue("s_scenarios{prjres}")
#   d_pr <- filter(scenarios, prjres == !!prjres)
#
#   s_pr        <- pull(d_pr, tif) %>% stack()
#   names(s_pr) <- pull(d_pr, label)
#   assign(v_pr, s_pr)
# }

# layers: features + scenarios ----
#prjres_lyrs <- unique(c(features$prjres, scenarios$prjres))
prjres_lyrs <- unique(scenarios$prjres)
# for (prjres in prjres_lyrs){ # prjres = ""
#   v_pr <- glue("s_layers{prjres}")
#   s_pr <- stack(
#     get(glue("s_features{prjres}")),
#     get(glue("s_scenarios{prjres}")))
#   assign(v_pr, s_pr)
# }

# setup ui choices ----
lyr_choices <- bind_rows(
    feature_labels,
    scenarios %>% select(type, label, prjres)) %>%
  select(type, label, prjres) %>%
  arrange(type, label, prjres)
# View(lyr_choices)

prjres_choices <- projections_tbl %>%
  mutate(
    lbl = glue("{name} {res}")) %>%
  filter(prjres %in% prjres_lyrs) %>%
  select(lbl, prjres)
