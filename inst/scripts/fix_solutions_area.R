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
