library(tidyverse)
library(raster)
library(fs)
library(glue)
library(here)

dir_ubc <- "~/Gdrive Ecoquants/projects/bbnj/data/raw/UBC-exploited-fish-projections"

# - Current_MCP1.tif: Maximum Catch Potential  (MCP) of more than 1,000 marine
# fish and invertebrates species in the high seas at the current status (20
# years average from 1995 to 2014) under GFDL 8.5 scenario
#   - under climate scenario?
#   - units in metric tons?
#   - list of species?
#   - individual species rasters to calculate profitability based on current market prices?
#
# - MCP2050_RCP85.tif: The projected Maximum Catch Potential of more than 1,000
# marine fish and invertebrates species in the high seas in the future (20 years
# average from 2041to 2060) under GFDL 8.5 scenario
#   - Why take average 2041to 2060 vs just 2050? A: to be similar to Current
#
# - hsmcpcap851.tif: The change in MCP in the 2050s relative to the current
# status under GFDL 8.5 scenario (I cap the change greater than 2 at 2).
#   - Why cap at change > 2?
#
# - MCP2050_RCP851.tif: ?

tifs <- list.files(dir_ubc, "*.tif$", full.names = T)
csvs <- list.files(dir_ubc, "*.csv$", recursive = T, full.names = T)

basename(tifs)
cat(paste(basename(tifs), collapse="\n"))

library(RColorBrewer)
pal <- colorRampPalette(brewer.pal(11, "Spectral"))
cols <- rev(pal(255))

for (tif in tifs){
  r <- raster(tif)
  title <- basename(tif) %>% path_ext_remove()
  plot(r, col = cols, main=title)
  #plot(log(r), col = cols, main=glue("log({title})"))
  hist(r)
}

r_now    <- raster(file.path(dir_ubc, "Current_MCP1.tif"))
r_future <- raster(file.path(dir_ubc, "MCP2050_RCP85.tif"))
r_dif0   <- raster(file.path(dir_ubc, "hsmcpcap851.tif"))
r_dif1   <- r_future - r_now
r_dif2   <- r_dif1
r_dif2[r_dif1 >  2] <-  2
r_dif2[r_dif1 < -1] <- -1

r_dif0
r_dif2
r_dif1

plot(r_dif0, col = cols, main="dif0")
plot(r_dif2, col = cols, main="dif2")
plot(r_dif1, col = cols, main="dif1")
hist(r_dif0, main="dif0")
hist(r_dif2, main="dif2")
hist(r_dif1, main="dif1")
