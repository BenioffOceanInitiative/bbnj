library(tidyverse)
library(raster)
library(here)
library(glue)
library(fs)
devtools::load_all()

dir_scenarios <- here("inst/app/www/scenarios")
rmds <- list.files(dir_scenarios, ".*\\.Rmd$")

#for (rmd in rmds){ # rmd = rmds[1]
for (rmd in rmds[c(-32)]){ # rmd = rmds[1]

  # 32/43: s03a.biofish.rls.effort.reclass.blm.grp01.30pct.gl.now.mol50km.Rmd -- 2019-08-07 08:42:44
  # label: solution
  # Quitting from lines 46-116 (s03a.biofish.rls.effort.reclass.blm.grp01.30pct.gl.now.mol50km.Rmd)
  # Error in .local(a, b, ...) :
  #   problem failed presolve checks. For more information see ?presolve_check

  pfx <- file.path(dir_scenarios, path_ext_remove(rmd))
  tif <- glue("{pfx}_sol.tif")
  rep <- glue("{pfx}_rep.csv")
  Rmd <- glue("{pfx}.Rmd")

  message(glue("{which(rmd == rmds)}/{length(rmds)}: {rmd} -- {Sys.time()}"))

  if (any(!file.exists(c(tif,rep)))){
    message(glue("  !file.exists(c(tif,rep)): {paste(!file.exists(c(tif,rep)), collapse=',')}"))
    next()
  }

  S <- get_tif_area_stats(tif)

  d <- bind_rows(
    tibble(
      feature = "_area_km2",
      absolute_held = S$solution_km2,
      relative_held = S$pct_solution),
    read_csv(rep) %>%
      filter(feature != "_area_km2"))

  write_csv(d, rep)

  rmarkdown::render(Rmd)
}

# b/c of merge conflicts
scenarios_redo <- c(
 "s01a.bio.now.mol50km.Rmd",
 "s02b.bio.future.mol50km.Rmd",
 "s03a.biofish.now.mol50km.Rmd",
 "s04a.biofish.alltime.mol50km.Rmd")
dir_scenarios <- here::here("inst/app/www/scenarios")
purrr::walk(scenarios_redo, function(rmd){
  path_rmd <- file.path(dir_scenarios, rmd)
  message(path_rmd)
  rmarkdown::render(path_rmd)})

# processing file: s01a.bio.now.mol50km.Rmd
# Error: 's05a.biofish.mcp3.rls.grp01.30pct.gl.alltime.mol50km_rep.csv' does not exist in current working directory ('/Users/bbest/github/bbnj/inst/app/www/scenarios').

