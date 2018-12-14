library(tidyverse)
library(fs)
library(glue)
library(scales)
#library(bbnj)
devtools::load_all()

dir_in  <- "/Users/bbest/github/sdg14-shiny/am-ranges/data/groups"
dir_out <- "/Users/bbest/Google Drive/projects/Pew BBNJ/data/derived/biodiversity"
res_am  <- 0.5
res_out <- 0.1
r_cellid_tif <- file.path(dir_data, glue("Caroline - high seas layer/high_seas_final_cellid_{res_out}dd.tif")) # see marxan.Rmd for creation

# *_be.tif: [b]inary (threshold = 0.4) * [e]xtinction risk raster [*.tif]
# extinction risk: NA=0.2, DD=0.2, LC=0.2, NT=0.4, VU=0.6, EN=0.8, CR=1, EX=1.2
# TODO biggest groups w/ db approach: coastal fishes (n=11,409), Crustaceans (n=3,058), Gastropods (n=2,852)

tifs_in  <- list.files(dir_in, ".*_be.tif$")
#tifs_cp  <- glue("{path_ext_remove(tifs_in)}_{res_am}dd.tif")
#file.copy(file.path(dir_in, tifs_in), file.path(dir_out, tifs_cp))
tifs_out <- glue("{path_ext_remove(tifs_in)}_{res_out}dd.tif")

r_cellid <- raster(r_cellid_tif) # plot(r_cellid)
r_cellid

for (i in 1:length(tifs_in)){ # i=3
  r <- raster(file.path(dir_in, tifs_in[i]))
  r <- rescale_grid(r, res_out) # plot(r)

  r  <- mask(r, r_cellid) # plot(r)
  values(r) <- rescale(values(r)) # plot(r)

  writeRaster(r, file.path(dir_out, tifs_out[i]), overwrite=T)
}
