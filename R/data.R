#' Raster of planning unit id for high seas
#'
#' Global half-degree raster of high seas with unique id.
#'
#' @format A \code{\link[raster]{raster}} with integer of id unique to cell.
"r_pu_id"

#' Polygon of high seas area
#'
#' Global area [-180,180,-90,90] with land and Exclusive Economic Zones (EEZs
#' frommMarineRegions.org) clipped out.
#'
#' @format A \code{\link[sf]{sf}} object of a polygon.
"p_highseas"

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
"s_bio_gmbi"

#' Raster stack from Global Fishing Watch analysis of high seas (Sala et al, 2018)
#'
#' Global half-degree raster of high seas. Year of analysis is 2016.
#'
#' @format A \code{\link[raster]{stack}} with layers of results from analyses.
"s_fish_gfw"

#' Raster stack from UBC (Cheung, Lam et al; in draft) of maximum catch
#' potential
#'
#' Global half-degree raster of high seas. Maximum catch potential (MCP;
#' landings in metric tons) of more than 1,000 fish and invertebrates species
#' for \code{mcp_2004} (average 1995 to 2014) and  \code{mcp_2050} (average 2041 to
#' 2060 under 'business as usual' climate change scenario GFDL 8.5).
#'
#' @format A \code{\link[raster]{stack}} with layers of results from analyses.
"s_fish_ubc"

#' Raster for count of seamounts (Kim & Wessel, 2011)
#'
#' Global half-degree raster of high seas.
#'
#' @format A \code{\link[raster]{raster}} count of seamount features within cell.
"r_phys_seamounts"

#' Raster stack for area (km2) of 1 thru 11 classes of benthic seascapes (Harris
#' & Whiteway, 2009)
#'
#' Global half-degree raster of high seas. Original 0.1 deg cells were
#' multiplied by area and summed to 0.5 deg cells.
#'
#' @format A \code{\link[raster]{raster}} integer class id of benthic seascape.
"s_phys_scapes"
