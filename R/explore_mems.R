library(tidyverse)
library(adespatial)
library(ade4)
library(sf)

# food web data
p <- st_read("data/marine.shp") %>% 
  transmute(foodweb.name = fdwb_nm,
            geographic.location = ggrphc_) %>% 
  left_join(read_csv("~/Proj/macro-net/env-net-dist/data/network_metrics.csv")) %>% 
  filter_all(all_vars(!is.na(.)))
net <- p %>% st_drop_geometry()
# location data
coords <- st_coordinates(p) %>% 
  as_tibble() %>% 
  transmute(lon = X, lat = Y)
d <- p %>% 
  st_drop_geometry() %>% 
  bind_cols(coords)
# MEMs ------------
mems <- read_rds("analysis-output/mems.rds")
mems$best.method$summary
mems_sel <- mems$best.method$MEM.select
plot(seq_len(length(mems_sel)),
     attr(mems_sel, "values"),
     type = "h",
     frame = FALSE,
     axes = FALSE,
     xlab = "",
     ylab = "Value")
axis(side = 1,
     at = seq_len(length(mems_sel)),
     labels = colnames(mems_sel),
     las = 2)
axis(side = 2,
     at = seq(floor(min(attr(mems_sel, "values"))),
              ceiling(max(attr(mems_sel, "values"))),
              by = 0.5),
     labels = seq(floor(min(attr(mems_sel, "values"))),
                  ceiling(max(attr(mems_sel, "values"))),
                  by = 0.5))
# plot MEMs in geographic space
s.value(d[, c("lon", "lat")],
        mems_sel[, colnames(mems_sel) == "MEM33"],
        adda = FALSE,
        grid = FALSE,
        sub = "MEM33")
