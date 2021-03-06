% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial.R
\name{gap_fill_raster}
\alias{gap_fill_raster}
\title{Gap-fill rasters with value based on function and mask}
\usage{
gap_fill_raster(r, fxn = "min", r_mask = r_pu_id)
}
\arguments{
\item{r}{raster input}

\item{fxn}{function, default="min"}

\item{r_mask}{raster mask, default=r_pu_id}
}
\value{
raster after filling NAs with result of function and mask
}
\description{
Gap-fill rasters with value based on function and mask
}
