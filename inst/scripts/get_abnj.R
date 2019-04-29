library(sf)
library(here)

# full resolution
download.file("https://github.com/ecoquants/bbnj/raw/master/data/p_abnj.rda", "p_abnj.rda")
load("p_abnj.rda") # 7.5 MB
plot(p_abnj["one"])

# simplified to 5% of vertices
download.file("https://github.com/ecoquants/bbnj/raw/master/data/p_abnj_s05.rda", "p_abnj_s05.rda")
load("p_abnj_s05.rda") # 0.3 MB
plot(p_abnj_s05["one"])
