#' Raster of planning unit id for high seas
#'
#' Global half-degree raster of high seas with unique id.
#'
#' @format A \code{\link[raster]{raster}} with integer of id unique to cell.
"r_pu_id"

#' Polygon of areas beyond national jurisdiction, ie high seas area
#'
#' Global area [-180,180,-90,90] with land and Exclusive Economic Zones (EEZs)
#' clipped out using the "Intersect_EEZ_IHO_v3_2018" product from MarineRegions.org.
#'
#' @format A \code{\link[sf]{sf}} object of a polygon.
"p_abnj"

#' Polygon of Pelagic Provinces of the World (Spalding et al, 2012) for high
#' seas
#'
#' See http://data.unep-wcmc.org/datasets/38. Clipped to high seas. TODO: merge
#' small areas with larger ones.
#'
#' @format A \code{\link[sf]{sf}} object of a polygon.
"p_abnj_ppow"

#' Simplified polygon of Pelagic Provinces of the World (Spalding et al, 2012) for high
#' seas
#'
#' See http://data.unep-wcmc.org/datasets/38. Clipped to high seas. TODO: merge
#' small areas with larger ones. Simplified to 5% of original vertices for fast visualization.
#'
#' @format A \code{\link[sf]{sf}} object of a polygon.
"p_abnj_ppow_s005"

#' Simplified polygon of areas beyond national jurisdiction, ie high seas area
#'
#' Global area [-180,180,-90,90] with land and Exclusive Economic Zones (EEZs)
#' clipped out using the "Intersect_EEZ_IHO_v3_2018" product from MarineRegions.org.
#' Simplified to 5% of original vertices for fast visualization.
#'
#' @format A \code{\link[sf]{sf}} object of a polygon.
"p_abnj_s05"

#' Polygons of exclusive economic zones (EEZs)
#'
#' Downloaded World_EEZ_v10_20180221/eez_v10.shp from MarineRegions.org. For
#' quick mapping applied simplification using
#' \code{\link[rmapshaper]{ms_simplify}(keep = 0.5)} twice.
#'
#' @format A \code{\link[sf]{sf}} object of a polygon.
"p_eez"

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

#' Raster stack from Global Fishing Watch analysis of high seas (Sala et al,
#' 2018)
#'
#' Global half-degree raster of high seas. Year of analysis is 2016. See also
#' \url{https://github.com/SFG-UCSB/The-economics-of-fishing-the-high-seas}.
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

#' Raster for count of hydrothermal vents (Interridge Vent Database v3.4)
#'
#' Global half-degree raster of high seas.
#'
#' @format A \code{\link[raster]{raster}} count of hydrothermal vents features within cell.
"r_phys_vents"

#' Raster for presence of mining claim (Interridge Vent Database v3.4)
#'
#' Global half-degree raster of presence in high seas, per ISA_claim_areas_update_20181202.
#'
#' @format A \code{\link[raster]{raster}} prence (1 or 0) of mine claim in cell.
"r_mine_claim"

#' Polygons of mine claims (ISA)
#'
#' Per ISA_claim_areas_update_20181202.
#'
#' @format A \code{\link[sf]{sf}} object of mining claim polygons.
"p_mine_claims"

#' Raster stack for area (km2) of 1 thru 11 classes of benthic seascapes (Harris
#' & Whiteway, 2009)
#'
#' Global half-degree raster of high seas. Original 0.1 deg cells were
#' multiplied by area and summed to 0.5 deg cells.
#'
#' @format A \code{\link[raster]{stack}} in which each layer represents area for given integer class id of benthic seascape corresponding to :
#' \enumerate{
#'   \item Upper Bathyal, shallow shelf, low DO, very high PP. thick sediment very warm
#'   \item Lower Bathyal, deep shelf (submerged), marginal plateaus, very high DO, high PP, thick sediment, warm
#'   \item Lower Bathyal, other ridges and plateaus, marginal plateaus, marginal seas with hilly bottom, steep, very low DO
#'   \item Lower Bathyal, continental slope see high PP, very thick sediment, warm
#'   \item Lower Bathyal, island arcs, steep, high DO
#'   \item Lower Bathyal (Abyssal-Hadal), deep water trenches, island arcs, trenches controlled by fracture zones, volcanic ridges and plateaus, very steep
#'   \item Abyssal, volcanic ridges and highs, central rift zone, ridge flanks, microcontinents, cold
#'   \item Abyssal, flat sedimented plains of marginal seas, central rift zone, ridge flanks, marginal seas with hilly bottoms, flat, low DO, thin sediment
#'   \item Abyssal (Hadal), trenches controlled by fracture zones, deep water trenches, large arched uplifted structures, low PP thin sediment, cold
#'   \item Abyssal, plains with slightly undulating seafloor, flat abyssal plains, continental rise, very flat, high DO, low PP, very cold
#'   \item Abyssal, hilly plains, large (arched) uplifted structures, flat abyssal plains, flat, very low PP, very thin sediment
#' }
"s_phys_scapes"
