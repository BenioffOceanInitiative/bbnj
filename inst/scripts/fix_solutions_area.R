library(tidyverse)
library(raster)
library(here)
library(glue)
library(fs)

dir_scenarios <- here("inst/app/www/scenarios")
rmds <- list.files(dir_scenarios, ".*\\.Rmd$")

for (rmd in rmds){ # rmd = rmds[1]
  pfx <- file.path(dir_scenarios, path_ext_remove(rmd))
  tif <- glue("{pfx}_sol.tif")
  rep <- glue("{pfx}_rep.csv")
  Rmd <- glue("{pfx}.Rmd")

  message(glue("{which(rmd == rmds)}/{length(rmds)}: {rmd} -- {Sys.time()}"))

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
