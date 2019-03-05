#' Raster stack of global marine biodiversity indicators from AquaMaps and IUCN
#'
#' Dataset containing global half-degree rasters clipped to the high seas of
#' species richness (\code{nspp_*}) and Red List Sum (\code{rls_*})  package
#' \href{https://github.com/marinebon/bbnj}{marinebon/bbnj}.
#'
#' @format A \code{\link[raster]{stack}} with indicator (\code{nspp_*} or \code{rls_*}) and taxonomic group or all taxa (\code{*_all}).
#' @source \url{https://github.com/marinebon/bbnj}
#' @source \url{https://aquamaps.org}
#' @source \url{https://iucnredlist.org}
"s_bio"

#' Raster of planning unit id for high seas
#'
#' Global half-degree raster of high seas with unique id.
#'
#' @format A \code{\link[raster]{raster}} with integer of id unique to cell.
"r_pu_id"

#' Raster stack from Global Fishing Watch analysis of high seas (Sala et al, 2018)
#'
#' Global half-degree raster of high seas with unique id. Year of analysis is 2016.
#'
#' @format A \code{\link[raster]{raster}} with integer of id unique to cell.
"s_fish_gfw"
