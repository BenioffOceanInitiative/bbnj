# bbnj

Biodiversity conservation for areas beyond national jurisdiction, funded by Pew with Benioff for UN

This repository is for storing code, including:

- [bbnj](https://benioffoceaninitiative.github.io/bbnj): R package of functions to import, analyze and visualize the conservation prioritization process
  - NOTE (2019-11-26): The current documentation, particularly the vignette [Explore prioritizr](https://benioffoceaninitiative.github.io/bbnj/articles/prioritizr_explore.html), needs to be updated with the latest functions using `pkgdown::build_site()`.
- [app](http://bbnj.ecoquants.com): Shiny app for interactive display of input layers and output solutions
  - NOTE (2019-11-26): This mapping application will be updated to use just one geographic projection (Mollweide) instead of three (geographic, web Mercator and Mollweide) and easily allow the selection of any input data layer and overlay with any output scenario.

Please also see this repository:

- [bbnj-scripts](https://github.com/BenioffOceanInitiative/bbnj-scripts): scripts to run various analyses, such as overlays between data layers and generation of manuscript figures.

## Install R package

From R console (e.g. RStudio IDE):

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
bbnj::run_app()
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
