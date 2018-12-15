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
# install shiny if you need it
# install.packages("shiny")

# run shiny app
shiny::runGitHub("ecoquants/bbnj", subdir="inst/app")
```

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

### Documentation

```r
pkgdown::build_site()
```
