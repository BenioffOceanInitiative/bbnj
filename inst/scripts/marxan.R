library(pacman)

p_load(
  "doParallel",
  "foreach",
  "foreign",
  "gplots",
  "Hmisc",
  "iptools",
  "labdsv",
  "leaflet",
  "maptools",
  "PBSmapping",
  "png",
  "rgdal",
  "rgeos",
  "rhandsontable",
  "rjson",
  "shiny",
  "shinyBS",
  "sp",
  "sqldf",
  "vegan",
  "xtable")

#Launching the user interfaces from R
library(marxanui)       # Load the R package

launch_app("import")    # Launch the import app to import your own data

launch_app("marxan")    # Launch the marxan app to run Marxan

launch_app("mxptest")   # Launch the parameter testing app to do BLM, SPF calibration, and target sensitivity testing

launch_app("marzone")   # Launch the marzone app to run MarZone

launch_app("manage")    # Launch the manage app to manage your datasets
