# bbnj
Biodiversity conservation for areas beyond national jurisdiction, funded by Pew with Benioff for UN

This repository is for storing code, including:

- R package of functions
- Shiny app for interactive display

## R package creation

```r
usethis::create_package(path = "../bbnj")
devtools::check()
```

### Include packages

```r
usethis::use_package("raster")
usethis::use_package("leaflet")
```

