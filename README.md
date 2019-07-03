# bbnj
Biodiversity conservation for areas beyond national jurisdiction, funded by Pew with Benioff for UN

This repository is for storing code, including:

- R package of functions
- Shiny app for interactive display

## Install R package

From R console (such as in RStudio IDE):

```r
# install devtools if you need it
# install.packages("devtools")

# install package
devtools::install_github("ecoquants/bbnj")

# load library
library(bbnj)

# get help
help(package="bbnj")
```

## Run Shiny app

```r
run_app()
```

## Scenario Building Notes

The scenarios are located under the package's `inst/app/www/scenarios`, which you can access via the local copy of the repository (for editing) or the installed package (for viewing only; ).



## Developer Notes

### R package creation

```r
usethis::create_package(path = "../bbnj")
devtools::check()
```

### Include packages

```r
usethis::use_package("raster")
usethis::use_package("sf")
usethis::use_package("leaflet")
```

### Vignettes

```r
usethis::use_vignette("prep_climate_ubc")
```

### Documentation

```r
pkgdown::build_site()
```
