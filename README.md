[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3554536.svg)](https://doi.org/10.5281/zenodo.3554536)

# bbnj

Biodiversity conservation for areas beyond national jurisdiction, funded by Pew with Benioff for UN

This repository is for storing code, including:

- [bbnj](https://benioffoceaninitiative.github.io/bbnj): R package of functions to import, analyze and visualize the conservation prioritization process
  - NOTE (2019-11-26): The current documentation, particularly the vignette [Explore prioritizr](https://benioffoceaninitiative.github.io/bbnj/articles/prioritizr_explore.html), needs to be updated with the latest functions using `pkgdown::build_site()`.
- [app](http://bbnj.ecoquants.com): Shiny app for interactive display of input layers and output solutions
  - NOTE (2019-11-26): This mapping application will be updated to use just one geographic projection (Mollweide) instead of three (geographic, web Mercator and Mollweide) and easily allow the selection of any input data layer and overlay with any output scenario.

Please also see this repository:

- [bbnj-scripts](https://github.com/BenioffOceanInitiative/bbnj-scripts): scripts to run various analyses, such as overlays between data layers and generation of manuscript figures.

## Peer-Reviewed Article

The outputs of this analysis were summarized here:

- Visalli, M. E., Best, B. D., Cabral, R. B., Cheung, W. W. L., Clark, N. A., Garilao, C., Kaschner, K., Kesner-Reyes, K., Lam, V. W. Y., Maxwell, S. M., Mayorga, J., Moeller, H. V., Morgan, L., Crespo, G. O., Pinsky, M. L., White, T. D., & McCauley, D. J. (2020). **Data-driven approach for highlighting priority areas for protection in marine areas beyond national jurisdiction**. _Marine Policy_, 122, 103927. https://doi.org/10.1016/j.marpol.2020.103927

Here are code and raw data outputs related to the figures:

### Figure 1

![](https://raw.githubusercontent.com/BenioffOceanInitiative/bbnj/master/inst/app/www/scenarios/s04a.biofish.alltime.mol50km_sol_levelplot_manual-edits.png)
- [s04a.biofish.alltime.mol50km.Rmd](https://github.com/BenioffOceanInitiative/bbnj/blob/master/inst/app/www/scenarios/s04a.biofish.alltime.mol50km.Rmd)\
  code to generate solution
- [s04a.biofish.alltime.mol50km\_sol.tif](https://github.com/BenioffOceanInitiative/bbnj/blob/master/inst/app/www/scenarios/s04a.biofish.alltime.mol50km_sol.tif)\
  raw raster data of solution

### Figure 2

![](https://raw.githubusercontent.com/BenioffOceanInitiative/bbnj-scripts/master/fig2/fig2-panels.png)\
Fig. 2. Summary of percent targets conserved by the solution highlighting priority areas to be considered for protection in marine areas beyond national jurisdiction (ABNJ) (Fig. 1). A minimum target of 30% of each conservation feature target was required in this implementation of the planning algorithm (dotted vertical line). Conservation features are described in Methods and Table S2 and are grouped here as: A. summary biological features reporting on ocean productivity, contemporary and future (i.e. climate change 2100 forecast) aggregated species richness, and the Red List sums of endangered species; B. benthic physical features; C. contemporary species richness aggregated by major taxonomic group; and D. future species richness (i.e. climate change 2100 forecast) aggregated by major taxonomic group. (For interpretation of the references to colour in this figure legend, the reader is referred to the Web version of this article.) 
- [figure2.Rmd](https://github.com/BenioffOceanInitiative/bbnj-scripts/blob/master/figure2.Rmd)\
    code to generate figure

### Figure 3

![](https://raw.githubusercontent.com/BenioffOceanInitiative/bbnj/master/inst/app/www/scenarios/diffs/s04a.biofish.alltime.mol50km%20-%20s02a.bio.alltime.mol50km.png)
Fig. 3. Influence of fishing effort as a cost in the prioritization process for protection of marine areas beyond national jurisdiction (ABNJ). This map compares model outputs that used fishing effort as the cost layer to a distinct model output that used uniform planning unit area as the cost layer and did not consider f ishing effort. Orange indicates regions that were dropped from the MPA solution when the model was parameterized to avoid areas of highest fishing effort. Dark blue indicates areas that were added as the model sought to identify new high priority regions for protection to counterbalance the loss of valued high fishing effort areas. Yellow indicates areas that were part of the planning solution regardless of whether fishing effort was included as the cost layer. (For interpretation of the references to colour in this figure legend, the reader is referred to the Web version of this article.) 

- [scenario_diffs.Rmd](https://github.com/BenioffOceanInitiative/bbnj-scripts/blob/master/scenario_diffs.Rmd)\
  source code to generate the figure
  * [`scenarios_diff()`](https://benioffoceaninitiative.github.io/bbnj/reference/scenarios_diff.html) ([source code](https://github.com/BenioffOceanInitiative/bbnj/blob/HEAD/R/prioritizr.R#L519-L572))\
  * [s04a.biofish.alltime.mol50km - s02a.bio.alltime.mol50km.tif](https://github.com/BenioffOceanInitiative/bbnj/blob/master/inst/app/www/scenarios/diffs/s04a.biofish.alltime.mol50km - s02a.bio.alltime.mol50km.tif)\
  raw raster data of difference


### Figure 4

![](https://raw.githubusercontent.com/BenioffOceanInitiative/bbnj/master/inst/app/www/scenarios/diffs/s02b.bio.future.mol50km%20-%20s01a.bio.now.mol50km.png)
Fig. 4. Influence of climate-biodiversity forecast data as additional targets in the prioritization process for protection of marine areas beyond national jurisdiction (ABNJ). This map compares conservation planning model outputs without versus with forecasted future (2100) biodiversity as an additional set of conservation feature target layers; both include contemporary biodiversity. Orange indicates regions that were dropped from the planning solution when the model was parameterized to include both contemporary and forecasted future biodiversity data. Dark blue indicates areas that were added as the model sought to identify new high priority regions for protection to offset loss of orange areas. Yellow indicates areas that were part of the planning solution regardless of whether forecasted future biodiversity data were included. (For interpretation of the references to colour in this figure legend, the reader is referred to the Web version of this article.)

- [scenario_diffs.Rmd](https://github.com/BenioffOceanInitiative/bbnj-scripts/blob/master/scenario_diffs.Rmd)\
  source code to generate the figure
  * [`scenarios_diff()`](https://benioffoceaninitiative.github.io/bbnj/reference/scenarios_diff.html) ([source code](https://github.com/BenioffOceanInitiative/bbnj/blob/HEAD/R/prioritizr.R#L519-L572))\
    function to calculate differences between scenarios
  * [s04a.biofish.alltime.mol50km - s03a.biofish.now.mol50km.tif](https://github.com/BenioffOceanInitiative/bbnj/blob/master/inst/app/www/scenarios/diffs/s04a.biofish.alltime.mol50km - s03a.biofish.now.mol50km.tif)\
  raw raster data of difference

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

## Acknowledgements

### Funding

This work has been funded by [Pew Charitable Trusts: Protecting Ocean Life on the High Seas](https://www.pewtrusts.org/en/projects/protecting-ocean-life-on-the-high-seas) with support by the [Benioff Ocean Initiative | UC Santa Barbara](https://boi.ucsb.edu/).

### Algorithm

Huge thanks to @jeffreyhanson for the excellent [prioritizr](https://prioritizr.net/) conservation prioritization R package that did the heavy lifting.

### Data

TODO: Many datasets to acknowledge here...

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
