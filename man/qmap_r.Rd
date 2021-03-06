% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial.R
\name{qmap_r}
\alias{qmap_r}
\title{Quick leaflet map of raster}
\usage{
qmap_r(r, r_title, na0 = T, method = "ngb",
  boundary_shp = "~/Gdrive Ecoquants/projects/bbnj/data/derived/boundary/high_seas_s05.shp")
}
\arguments{
\item{r}{raster}

\item{r_title}{title of raster layer}

\item{na0}{convert 0s to NA (default = FALSE)}

\item{method}{interpolation method for web mercator, "ngb" (nearest neighbor) or "bilinear" (default "ngb")}

\item{boundary_shp}{path to high seas boundary shapefile}
}
\description{
Quick leaflet map of raster
}
