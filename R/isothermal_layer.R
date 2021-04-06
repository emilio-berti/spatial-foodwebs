library(tidyverse)
library(sf)
library(ncdf4)
library(raster)
library(mgcv)

rasterOptions(tmpdir = "/work/eb97ziwi/tmp") #Rstudio server temporary folder
#' @param x vector (observed values)
#' @param r raster (derived values)
#' @param sf spatial object (food web locations)
check_values <- function(x, r, sf) {
  if (!max(x) == maxValue(r)) {
    message("Values got modified: ")
    message("Max value of XBT: ", max(x))
    message("Max value of XBT raster: ", maxValue(r))
  }
  if (!identical(crs(r), crs(sf))) {
    message("Reprojecting sf object to raster CRS")
    sf <- st_transform(sf, crs(r))
  }
  if (any(is.na(raster::extract(r, sf)))) {
    message("Missing extracted values: you need to interpolate")
  }
}
# non long-lat crs
EPSG <- rgdal::make_EPSG()
EPSG <- EPSG[!grepl("longlat", EPSG$prj4), ]
EPSG <- EPSG[grepl("world", EPSG$note, ignore.case = TRUE), ]
# load data ------------
# foodwebs coordinates
p <- st_read("data/marine.shp")
coords <- st_coordinates(p)
# isothermal layer database XBT
xbt <- nc_open("data/WODXBTtempNCEI.nc")
# isothermal layer depth from xbt ---------
xbt <- tibble(lon = ncvar_get(xbt, "Longitude"), #see Chu & Fan 2019
              lat = ncvar_get(xbt, "Latitude"),
              iso_depth = ncvar_get(xbt, "IsoLayerDepth"),
              year = ncvar_get(xbt, "year"),
              month = ncvar_get(xbt, "month")) %>% 
  mutate(lon = ifelse(lon > 180, -(360 - lon), lon)) %>% 
  group_by(lon, lat) %>% 
  summarise(avg_iso_depth = mean(iso_depth, na.rm = TRUE)) %>% 
  ungroup() %>% 
  st_as_sf(coords = c("lon", "lat"), crs = "EPSG:4326")
# project to world mercator in meters
xbt_mercator <- xbt %>% st_transform(crs = EPSG$prj4[1])
# rasterize xbt --------
r_template <- raster(xbt, nrows = 500, ncols = 500)
r_xbt <- rasterize(xbt, r_template, field = "avg_iso_depth", background = NA)
plot(r_xbt, col = topo.colors(100))
check_values(xbt$avg_iso_depth, r_xbt, p)
# global interpolation ----------
coords <- xyFromCell(r_xbt, 1:ncell(r_xbt))
values <- getValues(r_xbt)
d <- data.frame(X = coords[, 1], Y = coords[, 2], val = values)
# weights <- sqrt(weights + abs(min(weights)) * 1.1)
m <- gam(val ~ s(X, Y, k = 1000),
         data = d %>% 
           filter(!is.na(val)) %>% 
           distinct_all())
#land surface
land <- rnaturalearth::ne_countries() %>% 
  rasterize(r_xbt, field = 1, background = NA)
r_xbt[values(!is.na(land))] <- 0
# interpolation
preds <- d %>% 
  transmute(X, Y, Z = as.vector(predict(m, d))) %>%
  rasterFromXYZ(crs = crs(r_xbt))
preds[preds < 0] <- 0
p <- p %>% 
  mutate(isolayer_depth = extract(preds, p))
# plot
preds[values(!is.na(land))] <- NA
cairo_ps("figures/isolayer_depth.eps", width = 6, height = 4)
pal <- RColorBrewer::brewer.pal(9, "BuPu")
pal <- colorRampPalette(pal)(10)
rasterVis::levelplot(preds,
                     at = seq(0, maxValue(preds), length.out = 10),
                     maxpixels = 10^6,
                     margin = FALSE,
                     par.setting = rasterVis::rasterTheme(
                       region = pal,
                       panel.background = list(col = "grey20")
                     ))
dev.off()
r <- stack(preds, land)
names(r) <- c("iso.layer.depth", "land.mask")
writeRaster(r, "data/isothermal_layer_depth.tif", overwrite = TRUE)

# comparison original vs interpolated values -----
tibble(Original = d %>% 
         distinct_all() %>% 
         filter(!is.na(val)) %>% 
         pull(val),
       Interpolated = m$fitted.values) %>%
  filter(!is.na(Original)) %>% 
  ggplot() +
  aes(Original, Interpolated) +
  geom_bin2d() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_smooth() +
  scale_fill_gradient(low = "steelblue", high = "tomato") +
  scale_x_sqrt() +
  scale_y_sqrt() +
  theme_bw()
ggsave("figures/isolayer_depth_biplot.eps", width = 6, height = 4)
