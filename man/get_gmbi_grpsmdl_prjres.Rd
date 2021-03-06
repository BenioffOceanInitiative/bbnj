% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial.R
\name{get_gmbi_grpsmdl_prjres}
\alias{get_gmbi_grpsmdl_prjres}
\title{Get biogeographic stack, given taxonomic grouping, modeling period and projection-resolution}
\usage{
get_gmbi_grpsmdl_prjres(grpsmdl = "groups00", prjres = "", debug = F)
}
\arguments{
\item{prjres}{projection-resolution, per projections_tbl$prjres}

\item{debug}{defaults to False}

\item{dataset}{see datasets}
}
\value{
polygon (s), raster (r), or stack (s) per prefix of dataset in given projection-resolution (prjres)
}
\description{
Get biogeographic stack, given taxonomic grouping, modeling period and projection-resolution
}
