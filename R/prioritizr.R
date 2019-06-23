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
solve_log <- function(p, pfx=deparse(substitute(p)), redo=F){
  library(prioritizr)
  library(raster)
  select = dplyr::select
  library(glue)
  library(readr)

  log <- glue("{pfx}_log.txt")
  tif <- glue("{pfx}_sol.tif")
  rep <- glue("{pfx}_rep.csv")

  if (all(file.exists(log,tif,rep)) & !redo){
    s <- raster::raster(tif)
    attr(s, "feature_representation") <- read_csv(rep)

    message(glue("Files found ({basename(tif)}) & redo=F, returning previously computed raster solution"))
    return(s)
  }
  if (!file.exists(log)) message(glue("File NOT found 'log': {log}"))
  if (!file.exists(tif)) message(glue("File NOT found 'tif': {tif}"))
  if (!file.exists(rep)) message(glue("File NOT found 'log': {rep}"))

  sink(log)
  s <- prioritizr::solve(p)
  sink()

  if ("RasterLayer" %in% class(s)){
    raster::writeRaster(s, tif, overwrite=T)

    d <- feature_representation(p, s)
    write_csv(d, rep)
    attr(s, "feature_representation") <- d
  }

  s
}

#' Display table of target representaiton in conservation problem solution
#'
#' @param csv comma-seperated value output from \code{\link{solve_log}}
#'
#' @return \code{\link[formattable]{formattable}} in form of
#'   \code{\link[DT]{datatable}} with horizontal bars indicating percent held
#'   per target.
#' @export
#'
#' @examples
tbl_target_representation <- function(csv = glue("{pfx}_rep.csv")){
  # csv = here("inst/scenarios/s01a.bio.10.gl.now_rep.csv")
  library(readr)
  library(formattable)
  library(dplyr)
  library(DT)

  d <- readr::read_csv(csv)

  stopifnot(c("feature", "relative_held", "absolute_held") %in% names(d))

  d %>%
    dplyr::mutate(
      relative_held = ifelse(is.nan(relative_held), 0, relative_held),
      relative_held = formattable::percent(relative_held),
      absolute_held = formattable::accounting(absolute_held)) %>%
    formattable(list(
      relative_held = formattable::normalize_bar())) %>%
    formattable::as.datatable(
      options = list(
        pageLength = 1000,
        columnDefs = list(list(className = 'dt-right', targets = c(2,3)))))
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
