library(tidyverse)
library(sf)

marine <- read_csv("~/Proj/gateway/gateway.csv", guess_max = 200000) %>% 
  filter(ecosystem.type == "marine") %>%
  select(foodweb.name, geographic.location, longitude, latitude) %>% 
  distinct_all()

p <- marine %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326")

st_write(p, "data/marine.shp", delete_layer = TRUE)

p <- st_read("data/marine.shp")
p1 <- p[1:ceiling(nrow(p) / 4), ]
p2 <- p[ceiling(nrow(p) / 4 + 1):ceiling(2 * nrow(p) / 4), ]
p3 <- p[ceiling(2 * nrow(p) / 4 + 1):ceiling(3 * nrow(p) / 4), ]
p4 <- p[ceiling(3 * nrow(p) / 4 + 1):nrow(p), ]
st_write(p1, "data/to_GEE_1.shp", delete_layer = TRUE)
st_write(p2, "data/to_GEE_2.shp", delete_layer = TRUE)
st_write(p3, "data/to_GEE_3.shp", delete_layer = TRUE)
st_write(p4, "data/to_GEE_4.shp", delete_layer = TRUE)
