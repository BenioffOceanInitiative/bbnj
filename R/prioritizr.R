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
solve_log <- function(p, pfx=deparse(substitute(p)), redo=F, debug=F, ...){
  library(prioritizr)
  library(raster)
  select = dplyr::select
  library(glue)
  library(readr)
  library(lwgeom)

  log <- glue("{pfx}_log.txt")
  tif <- glue("{pfx}_sol.tif")
  rep <- glue("{pfx}_rep.csv")

  if (all(file.exists(log,tif,rep)) & !redo){
    r_sol <- raster::raster(tif)
    attr(r_sol, "feature_representation") <- read_csv(rep)

    if (debug) message(glue("Files found ({basename(tif)}) & redo=F, returning previously computed raster solution"))
    return(tif)
  }
  if (debug){
    if (!file.exists(log)) message(glue("File NOT found 'log': {log}"))
    if (!file.exists(tif)) message(glue("File NOT found 'tif': {tif}"))
    if (!file.exists(rep)) message(glue("File NOT found 'log': {rep}"))
  }

  sink(log)
  r_sol <- prioritizr::solve(p, ...)
  sink()

  if ("RasterLayer" %in% class(r_sol)){
    raster::writeRaster(r_sol, tif, overwrite=T)

    #P <- get_tif_projection(tif)
    #if (P$prj != "mer")
    attr(r_sol, "shapefile_gcs") <- tif_to_shp_gcs(tif)

    d <- feature_representation(p, r_sol)

    S <- get_tif_area_stats(tif)

    d_a <- tibble(
      feature = "_area_km2",
      absolute_held = S$solution_km2,
      relative_held = S$pct_solution)
    d <- bind_rows(d_a, d)

    write_csv(d, rep)
    attr(r_sol, "feature_representation") <- d
  }

  tif
}

#' get area statistics from tif of solution
#'
#' @param tif
#'
#' @return named list of area stats for high seas (\code{!is.na(value)}) and solution (\code{value=1})
#' @export
get_tif_area_stats <- function(tif){
  P <- get_tif_projection(tif)
  r_sol <- raster(tif)

  if (P$prj == "gcs"){
    #tif <- "/Users/bbest/github/bbnj/inst/app/www/scenarios/s00b.bio.30pct.gl.gcs0.5d_sol.tif"
    #r_sol <- raster(tif) # plot(r_sol==1); plot(r_sol)
    r_hs_a <- area(r_sol) %>%
      mask(r_sol, maskvalue=NA)
    r_sol_a <- r_hs_a %>%
      mask(r_sol, maskvalue=0)
    a_hs  <- sum(values(r_hs_a), na.rm=T)
    a_sol <- sum(values(r_sol_a), na.rm=T)
  } else {
    # tif <- "/Users/bbest/github/bbnj/inst/app/www/scenarios/s00a.bio.30pct.gl.mol100km_sol.tif"
    # r_sol <- raster(tif) # plot(r_sol==1); plot(r_sol)
    a_cell <- prod(res(r_sol))/(1000*1000)
    a_hs   <- sum(values(!is.na(r_sol))) * a_cell
    a_sol  <- sum(values(r_sol == 1), na.rm=T) * a_cell
  }

  A <- list(
    highseas_km2 = a_hs,
    solution_km2 = a_sol,
    pct_solution = a_sol/a_hs)
  A$summary = glue(
      "Solution : {str_pad(format(A$solution_km2, big.mark=','), nchar(format(A$highseas_km2, big.mark=',')))} km2 ({round(A$pct_solution*100, 2)}%)
       High seas: {format(A$highseas_km2, big.mark=',')} km2")


  A
}

#' get raster solution of tif
#'
#' @param tif
#'
#' @return raster solution of tif with attributes appended
#' @export
get_tif_solution <- function(tif){
  pfx <- str_replace(tif, "_sol.tif", "")
  log <- glue("{pfx}_log.txt")
  rep <- glue("{pfx}_rep.csv")
  shp <- glue("{fs::path_ext_remove(tif)}_gcs.shp")

  r_sol <- raster(tif)
  attr(r_sol, "feature_representation") <- read_csv(rep)
  attr(r_sol, "log")                    <- readLines(log)
  attr(r_sol, "area_stats")             <- get_tif_area_stats(tif)
  attr(r_sol, "projection")             <- get_tif_projection(tif)
  attr(r_sol, "shapefile_gcs")          <- shp

  r_sol
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

#' Report on conservation solution
#'
#' @param tif path to raster tif of solution
#'
#' @return for use in reports to spit out area stats, map and table of resulting scenario
#' @export
report_solution <- function(tif, fig_ht_in=2, redo=F){
  # tif <- "~/github/bbnj/inst/app/www/scenarios/s00a.bio.30pct.gl.gcs0.5d_sol.tif"

  # tif; redo=redo_map; fig_ht_in=2

  pfx     <- str_replace(tif, "_sol.tif", "")
  sol_png <- glue("{pfx}_sol_map.png")
  #sol_pdf <- glue("{fs::path_ext_remove(tif)}_map.pdf")
  r_sol   <- get_tif_solution(tif)
  P       <- get_tif_projection(tif)
  A       <- get_tif_area_stats(tif)

  # area ----
  cat(A$summary)
  #(markdown::renderMarkdown(text = "Hello World!"))

  # plot ----
  if (!file.exists(sol_png) | redo){
    # countries <- rnaturalearth::ne_countries(returnclass = "sf") %>%
    #   st_transform(P$epsg)
    # graticules <- st_graticule(countries)

    map_r2png(r_sol, sol_png)

    # map_sol <- function(){
    #   #op <- par(mar = rep(0, 4))
    #   #op <- par(bg=NA,mar=c(0,0,0,0),oma=c(0,0,0,0))
    #   plot(r_sol, legend=F, axes=F, box=F)
    #   plot(st_geometry(countries), add=T, col=gray(0.8), border=gray(0.7), lwd=0.5)
    #   plot(st_geometry(graticules), add=T, col=gray(0.6), lwd=0.5)
    #   #par(op)
    # }
    # png(sol_png, width=480*4, height = 480*4, res=300, type="cairo", units='px')
    # map_sol()
    # dev.off()
    # sol_png %>%
    #   magick::image_read() %>% magick::image_trim() %>%
    #   magick::image_write(sol_png)
  }
  img <- magick::image_read(sol_png) %>%
    grid::grid.raster()

  # table ----
  tbl_target_representation(glue("{pfx}_rep.csv"))
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

  readr::write_csv(d, csv)
  d
}

#' Calculate difference of two raster solutions
#'
#' @param r1 1st raster solution [0,1]
#' @param r2 2nd raster solution [0,1]
#'
#' @return r2 - r1
#' @export
#'
#' @examples
r_diff <- function(r1, r2){
  rng1 <- c(cellStats(r1, "min"), cellStats(r1, "max"))
  rng2 <- c(cellStats(r2, "min"), cellStats(r2, "max"))
  stopifnot(identical(rng1, c(0,1)))
  stopifnot(identical(rng2, c(0,1)))

  r_d  <- r2*10 - r1
  # plot(r_d, main="r_d")
  # table(values(r_d), useNA = "ifany")

  d_subs <- tribble(
    ~old,   ~new,
    0, -9999,
    9,     0,
    -1,    -1,
    10,     1)

  r_ds <- subs(r_d, d_subs, "old", "new")
  # table(values(r_ds), useNA = "ifany")
  # plot(r_ds, main="r_ds")
  r_dsm <- mask(r_ds, r_ds, maskvalue=-9999)
  # table(values(r_dsm), useNA = "ifany")
  # plot(r_dsm, main="r_dsm")

  r_dsm
}

#' map raster to png
#'
#' @param r
#' @param png
#'
#' @return
#' @export
#'
#' @examples
map_r2png <- function(r, png){
  library(RColorBrewer)
  library(sf)
  library(raster)
  library(magick)
  # library(tidyverse)
  # library(glue)
  # library(here)
  # scenarios_diff:  r <- rd   ; png <- png
  # report_solution: r <- r_sol; png <- sol_png

  col_hs <- ggplot2::alpha("lightblue", 0.3)

  # projection
  P <- get_r_projection(r)

  # overlays
  countries  <- rnaturalearth::ne_countries(returnclass = "sf") %>%
    st_transform(P$epsg)
  graticules <- st_graticule(
    lon = seq(-180, 180, by = 20),
    lat = seq(-80, 80, by = 20)) %>%
    st_transform(P$epsg)

  rng <- c(cellStats(r, "min"), cellStats(r, "max"))

  lgnd_inset4png <- 0.21

  plot_countries_graticules <- function(){
    plot(st_geometry(countries), add=T, col=gray(0.8), border=gray(0.7), lwd=0.5)
    plot(st_geometry(graticules), add=T, col=gray(0.6), lwd=0.5)
  }

  # begin plot
  png(png, width=480*4, height = 480*4, res=300, type="cairo", units='px')

  # plot difference
  if (identical(rng, c(-1,1))){
    #cols <- brewer.pal(5, "Set1")[1:3] # red, blue, green
    # override with: http://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=3
    cols <- c("#fc8d59", "#ffffbf", "#91bfdb")
    pal <- colorRampPalette(cols)(3)

    r_hs <- get_d_prjres("r_pu_id", P$prjres)
    r_hs[!is.na(r_pu_id)] = 1

    par(xpd = F)
    plot(r_hs, legend=F, axes=F, box=F, col=col_hs)

    plot(r, legend=F, axes=F, box=F, col=pal, add=T)

    plot_countries_graticules()

    par(xpd = T)
    legend(
      "bottom",
      legend = c("Loss", "Same", "Gain", "ABNJ"),
      fill = c(cols, col_hs),
      horiz = TRUE,
      cex = 0.6,
      inset = lgnd_inset4png)
  }
  # plot solution
  if (identical(rng, c(0, 1))){
    # cols <- terrain.colors(2) # gray, green
    # override gray
    col_sol <- ggplot2::alpha("darkgreen",0.7)
    cols <- c(col_hs, col_sol) # lightblue, darkgreen
    #pal <- colorRampPalette(cols)(2)

    par(xpd = F)
    plot(r, legend=F, axes=F, box=F, col=cols)

    plot_countries_graticules()

    par(xpd = T)
    legend(
      "bottom",
      legend = c("Areas Beyond National Jurisdiction", "Planning Analysis Solution"),
      fill = c(cols),
      horiz = TRUE,
      cex = 0.5,
      inset = lgnd_inset4png)
      #inset = -.0001)
  }

  # end plot
  dev.off()
  #browseURL(png)

  png %>%
    magick::image_read() %>%
    magick::image_trim() %>%
    magick::image_write(png)
  #browseURL(png)
}

#' calculate scenario difference
#'
#' @param scenarios a named list with scenarios
#' @param dir_scenarios directory to find scenarios
#' @param dir_diffs directory to output png difference
#'
#' @return
#' @export
#'
#' @examples
scenarios_diff <- function(scenarios, dir_scenarios, dir_diffs){
  # scenarios = scenarios[c("s1", "s2a")]
  # scenarios = scenarios[c("s2a", "s4")]
  stopifnot(length(scenarios) == 2)

  d <- tibble(
    s          = names(scenarios),
    scen       = unlist(scenarios),
    tif        = file.path(dir_scenarios, glue("{scen}_sol.tif")),
    tif_exists = file.exists(tif))
  stopifnot(all(d$tif_exists))

  png <- glue("{dir_diffs}/{scenarios[2]} - {scenarios[1]}.png")

  r1 <- raster(d$tif[1])
  r2 <- raster(d$tif[2])

  rd <- r_diff(r1, r2)

  a_km2 <- prod(res(rd)) / (1000 * 1000)

  table(values(rd))

  tbl <- tibble(
    value = na.omit(values(rd))) %>%
    group_by(value) %>%
    summarize(
      ncells = n()) %>%
    mutate(
      label = recode(
        value,
        `-1` = "loss",
        `0`  = "same",
        `1`  = "gain"),
      area_km2 = ncells * a_km2,
      pct = area_km2 / sum(area_km2)) %>%
    dplyr::select(value, label, ncells, area_km2, pct)

  map_r2png(rd, png)

  list(
    tbl = tbl,
    png = png)
}

