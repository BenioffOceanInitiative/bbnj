#' Helper function to run prioritizr::solve() and log output and resulting raster
#'
#' @param p \code{\link[prioritizr]{problem}} to feed into \code{\link[prioritizr]{solve}}
#' @param log path to output log from solving problem
#'
#' @return output from \code{\link[prioritizr]{solve}}
#' @export
#' @import prioritizr raster
#'
#' @examples
#solve_log <- function(p, log, tif){
solve_log <- function(p, pfx= deparse(substitute(p)), redo=F){
  library(prioritizr)
  library(raster)
  select = dplyr::select
  library(glue)

  log <- glue("{pfx}_log.txt")
  tif <- glue("{pfx}_sol.tif")

  if (all(file.exists(log,tif)) & !redo){
    s <- raster::raster(tif)
    message(glue("Files found ({basename(tif)}) & redo=F, returning previously computed raster solution"))
    return(s)
  }

  sink(log)
  s <- prioritizr::solve(p)
  sink()

  if ("RasterLayer" %in% class(s)){
    raster::writeRaster(s, tif, overwrite=T)
  }

  s
}

#' Produce diagnostics table for setting relative targets of features
#'
#' @param pu planning unit raster
#' @param features features stack
#' @param budget in percent of area
#'
#' @return table with rel_all and rel_each per feature
#' @export
#'
#' @examples
problem_diagnostics <- function(pu, features, budget, pfx= deparse(substitute(features)), redo=F){
  # pu = r_pu_areas; features = p01_features; budget = 0.1

  library(prioritizr)
  library(raster)
  select = dplyr::select
  library(glue)

  log <- glue("{pfx}_diagnostics_log.txt")
  csv <- glue("{pfx}_diagnostics.csv")

  if (all(file.exists(log,csv)) & !redo){
    message(glue("Files found ({basename(csv)}) & redo=F, returning previously computed diagnostic table"))
    return(read_csv(csv))
  }

  sink(log)

  # get feature representation for maximizing all features given budget
  p <- problem(pu, features) %>%
    add_max_utility_objective(budget = budget) %>%
    add_gurobi_solver()
  p_sol <- solve(p)
  d <- feature_representation(p, p_sol) %>%
    dplyr::select(feature, rel_all=relative_held) %>%
    dplyr::mutate(rel_each = NA)
  # d0 <- d

  for (i in 1:nlayers(features)){ # i=5
    message(glue("{i}: {names(features[[i]])}"))

    feature_i <- raster(features, i)
    #options(error = utils::recover)
    #options(error = NULL)
    p <- problem(pu, feature_i) %>%
      add_max_utility_objective(budget = budget) %>%
      add_gurobi_solver()

    p_sol <- solve(p)
    d$rel_each[[i]] = feature_representation(p, p_sol) %>%
      dplyr::pull(relative_held)
  }

  sink()

  write_csv(d, csv)
  d
}
