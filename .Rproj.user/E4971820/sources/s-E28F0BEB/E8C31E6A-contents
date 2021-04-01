library(tidyverse)
library(sf)
library(ncdf4)
library(raster)
library(automap) #for kriging interpolation

# foodwebs coordinates
p <- st_read("data/marine.shp")
coords <- st_coordinates(p)
# isothermal layer database XBT
xbt <- nc_open("/home/eb97ziwi/Documents/databases/NOAA_thermoclines/1.1/data/0-data/WODXBTtempNCEI.nc")
# isothermal layer depth from xbt ---------
xbt <- tibble(lon = ncvar_get(xbt, "Longitude"), #see Chu & Fan 2019
              lat = ncvar_get(xbt, "Latitude"),
              iso_depth = ncvar_get(xbt, "IsoLayerDepth"),
              year = ncvar_get(xbt, "year")) %>% 
  mutate(lon = ifelse(lon > 180, -(360 - lon), lon)) %>% 
  group_by(lon, lat) %>% 
  summarize(avg_iso_depth = mean(iso_depth, na.rm = TRUE)) %>% 
  ungroup() %>% 
  st_as_sf(coords = c("lon", "lat"), crs = "EPSG:4326")
# rasterize ctd & xbt --------
r_template <- raster(xbt, nrows = 1000, ncols = 1000)
r_xbt <- rasterize(xbt, r_template, field = "avg_iso_depth", background = NA)
r_xbt <- aggregate(r_xbt, 4, fun = mean)
plot(r_xbt, col = topo.colors(100))
land <- rnaturalearth::ne_countries() %>%
  rasterize(r_template, field = 1)
# interpolate
v <- getValues(r_xbt)
xy <- xyFromCell(r_xbt, 1:ncell(r_xbt))
tps <- fields::Tps(xy, v) # thin plate spline regression
r_ctd_tps <- interpolate(r_ctd, tps)
# comparison
tibble(original = v,
       interpolated = getValues(r_ctd_tps)) %>% 
  filter(!is.na(original)) %>% 
  ggplot() +
  aes(original, interpolated) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.2) +
  geom_smooth() +
  theme_bw()
par(mfrow = c(1, 2))
plot(r_ctd)
plot(mask(r_ctd_tps, land, inverse = TRUE))
# foodweb values --------
extract(r_ctd_tps, p)
