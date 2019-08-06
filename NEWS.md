# bbnj 0.7.2

* in process fixing `_area_km2` in tbl of solutions
* + tbl gain/loss/same for `scenarios_diff()`
* fixed `sol_gcs.shp` w/ `ms_dissolve()`
* suggestion for relabel seamounts layers [#21](https://github.com/ecoquants/bbnj/issues/21)
* + `scenario_overlays.Rmd` w/ indiv EBSAs & ISAs in inst/scripts
* + candidate scenarios for manuscript

# bbnj 0.7.1

* + RLS for groups00* #19

# bbnj 0.7.0

* updated mining claims (`r_min_claims`,`p_min_claims`) to 20190719 (from 20181202) [#14](https://github.com/ecoquants/bbnj/issues/14)
* more scenarios
* +jellyfish, -squids, ∆cephalopods groups01-03 #3
* new treemap #15
* include both Red List Index (`rli`) and Red List Sum (`rls`) for `s_bio_gmbi_*`
* `get_gmbi_grpsmdl_prjres()` to retrieve `s_bio_gmbi_*` by taxonomic groups0# | groups0#_2100 for projection-resolution
* new functions for mapping differences between scenarios: `scenarios_diff_png()`, `r_diff()`, 
`map_r2png()`

# bbnj 0.6.1

* added `s_fish_saup_v2` w/out Chilean jack mackerel

# bbnj 0.6.0

* added `r_phys_scapes_hetero`

# bbnj 0.5.4

* fixed `get_tif_projection()` & prepending `_area_km2`

# bbnj 0.5.3

* fixed input features in app with `raster::readAll()`

# bbnj 0.5.2

* rendering 300dpi maps in scenario report

# bbnj 0.5.1 fixed `pu_id` projections

# bbnj 0.5.0

* now working with different projections

# bbnj 0.4.1

# bbnj 0.4

* Fixed modal window to resize and include full report
* Move scenarios into `inst/app/www/scenarios`
* Load layers on fly, so new scenario shows up with app restart

# bbnj 0.3

* Added `r_vgpm`

# bbnj 0.2

* Added `r_phys_vents` dataset, `problem_diagnostics()`

# bbnj 0.1

* Added a `NEWS.md` file to track changes to the package.
