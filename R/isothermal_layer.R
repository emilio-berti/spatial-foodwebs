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

# local interpolation ----------
# each point location is interpolated using values within a buffer with radius =
# 500 km at 10 km resolution
vals <- data.frame(original = raster::extract(r_xbt, p))
interpolated <- list()
for (x in seq_len(nrow(vals))) {
  geom <- p[x, ] %>% 
    st_transform(crs = EPSG$prj4[1]) %>% 
    st_buffer(500000)
  lon_range <- as.vector(st_bbox(geom))[c(1, 3)]
  lat_range <- as.vector(st_bbox(geom))[c(2, 4)]
  r_tmp <- raster(res = 10000,
                  xmn = lon_range[1], xmx = lon_range[2],
                  ymn = lat_range[1], ymx = lat_range[2],
                  crs = crs(geom))
  r_xbt <- xbt_mercator %>% 
    mutate(X = st_coordinates(.)[, 1],
           Y = st_coordinates(.)[, 2]) %>% 
    filter(X >= lon_range[1], X <= lon_range[2],
           Y >= lat_range[1], Y <= lat_range[2]) %>%
    rasterize(r_tmp, field = "avg_iso_depth", background = NA)
  plot(r_xbt)
  # gam interpolation
  coords <- xyFromCell(r_xbt, 1:ncell(r_xbt))
  values <- getValues(r_xbt)
  d <- data.frame(X = coords[, 1], Y = coords[, 2], val = values)
  k_max <- d %>% 
    filter(!is.na(val)) %>% 
    distinct_all() %>% 
    nrow()
  if (k_max > 500) {
    k_max <- 500
  }
  # weights <- sqrt(weights + abs(min(weights)) * 1.1)
  m <- gam(val ~ s(X, Y, k = k_max),
           data = d %>% 
             filter(!is.na(val)) %>% 
             distinct_all())
  land <- rnaturalearth::ne_countries() %>% #land surface
    rasterize(r_xbt, field = 1, background = NA)
  r_xbt[values(!is.na(land))] <- 0
  preds <- d %>% 
    transmute(X, Y, Z = as.vector(predict(m, d))) %>%
    rasterFromXYZ(crs = crs(r_xbt))
  preds[values(!is.na(land))] <- 0
  preds[1] <- maxValue(r_xbt)
  # cairo_pdf(paste0("analysis-output/interpolations/point_", x, ".pdf"),
  #           width = 6, height = 4)
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 4))
  #rasters
  plot(r_xbt, col = c("grey80", rev(topo.colors(100))))
  plot(preds, col = c("grey80", rev(topo.colors(100))))
  plot(d$val, predict(m, d), #scatter plot
       xlab = "Original value",
       ylab = "Interpolated",
       frame = FALSE,
       pch = 20,
       col = adjustcolor("grey20", alpha.f = 0.1))
  abline(0, 1, lt = 2)
  fit <- loess(predict(m, d) ~ d$val, span = 1)
  lines(sort(fit$x), fit$fitted[order(fit$x)], col = "tomato", lw = 2)
  dev.off()
  #results to list
  ans <- data.frame(
    p.row = x,
    xbt.value = raster::extract(r_xbt, 
                                p[x, ] %>% 
                                  st_transform(crs = EPSG$prj4[1]) %>% 
                                  st_buffer(2000)) %>% 
      unlist() %>% 
      mean(na.rm = TRUE),
    avg.value = raster::extract(preds,
                                p[x, ] %>% 
                                  st_transform(crs = EPSG$prj4[1]) %>% 
                                  st_buffer(2000)
    ) %>% 
      unlist() %>% 
      mean(na.rm = TRUE)
  )
  interpolated <- append(interpolated, list(ans))
}

interpolated %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  ggplot() +
  aes(xbt.value, avg.value) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_smooth() +
  scale_x_sqrt() +
  scale_y_sqrt() +
  theme_bw()

# interpolation on whole dataset ----------
d <- xbt_mercator %>% #get unique values and coordinates
  mutate(X = st_coordinates(.)[, 1],
         Y = st_coordinates(.)[, 2]) %>% 
  st_drop_geometry() %>% 
  filter(!is.na(avg_iso_depth)) %>% 
  distinct_all()
r_tmpl <- raster(crs = EPSG$prj4[1], res = 10000, #raster template
                 xmn = min(d$X), xmx = max(d$X),
                 ymn = min(d$Y), ymx = max(d$Y))
global_raster <- xbt_mercator %>% #rasterize dataset
  mutate(X = st_coordinates(.)[, 1],
         Y = st_coordinates(.)[, 2]) %>% 
  rasterize(r_tmpl, field = "avg_iso_depth", background = NA)

grid <- st_make_grid(xbt_mercator, n = c(10, 10)) %>% 
  st_as_sf() %>% 
  rasterize(r_tmpl)
grid_id <- unique(getValues(grid))

for (g in grid_id) {
  square <- grid
  square[square != g] <- NA
  square <- trim(square)
  lon_range <- as.numeric(extent(square)[1:2])
  lat_range <- as.numeric(extent(square)[3:4])
  d <- xbt_mercator %>% 
    mutate(X = st_coordinates(.)[, 1],
           Y = st_coordinates(.)[, 2]) %>% 
    filter(X > lon_range[1], X < lon_range[2],
           Y > lat_range[1], Y < lat_range[2]) %>% 
    distinct_all() %>% 
    relocate(avg_iso_depth, .after = Y)
  raster_d <- rasterize(as_Spatial(d), square, field = "avg_iso_depth", background = NA)
  d <- d %>% st_drop_geometry()
  m <- gam(avg_iso_depth ~ s(X, Y, k = 500), data = d)
  # scatter plot
  plot(predict(m, d), d$avg_iso_depth,
       xlab = "Interpolated",
       ylab = "Original value",
       frame = FALSE,
       pch = 20,
       col = adjustcolor("grey20", alpha.f = 0.1))
  abline(0, 1, lt = 2)
  fit <- loess(d$avg_iso_depth ~ predict(m, d))
  lines(sort(fit$x), fit$fitted[order(fit$x)], col = "tomato", lw = 2)
  xy <- xyFromCell(raster_d, 1:ncell(raster_d)) %>% 
    as_tibble() %>% 
    transmute(X = x, Y = y)
  preds <- predict(m, xy)
  interpolated <- rasterFromXYZ(cbind(xy, preds), crs = crs(square))
  l <- land %>% 
    projectRaster(square) %>% 
    crop(square)
  interpolated[!is.na(l)] <- NA
  plot(interpolated, col = rev(topo.colors(100)))
}

grid <- global_raster %>% 
  xyFromCell(1:1000) %>%
  as.matrix() %>% 
  list() %>% 
  SpatialMultiPoints() %>% 
  points2grid()


makegrid(global_raster, n = 10)


weights <- d %>% 
  pull(avg_iso_depth) / sum(d$avg_iso_depth, na.rm = TRUE)
m <- gam(avg_iso_depth ~ s(X, Y, k = 1000),
         data = d,
         weights = weights)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 4))
#rasters
plot(r_xbt, col = c("grey80", rev(topo.colors(100))))
plot(preds, col = c("grey80", rev(topo.colors(100))))
plot(d$val, predict(m, d), #scatter plot
     xlab = "Original value",
     ylab = "Interpolated",
     frame = FALSE,
     pch = 20,
     col = adjustcolor("grey20", alpha.f = 0.1))
abline(0, 1, lt = 2)
fit <- loess(predict(m, d) ~ d$val)
lines(sort(fit$x), fit$fitted[order(fit$x)], col = "tomato", lw = 2)# plots



EPSG <- rgdal::make_EPSG()
EPSG <- EPSG[!grepl("longlat", EPSG$prj4), ]
EPSG <- EPSG[grepl("world", EPSG$note, ignore.case = TRUE), ]
xbt_transformed <- xbt #%>% st_transform(crs = EPSG$prj4[1])
xbt_transformed <- xbt_transformed %>% 
  mutate(X = st_coordinates(xbt_transformed)[, 1],
         Y = st_coordinates(xbt_transformed)[, 2]) %>% 
  st_drop_geometry()
# gam
xbt_gam <- gam(avg_iso_depth ~ s(X, Y, k = 200), data = xbt_transformed)
# comparison original vs interpolated values
tibble(Original = xbt_transformed$avg_iso_depth,
       Interpolated = xbt_gam$fitted.values) %>% 
  filter(!is.na(Original)) %>% 
  ggplot() +
  aes(Original, Interpolated) +
  geom_bin2d() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_vline(xintercept = max(foodweb_values), linetype = "dotted") +
  geom_smooth() +
  scale_x_continuous(limits = c(0, 200)) +
  scale_fill_gradient(low = "steelblue", high = "tomato") +
  theme_bw()
interp_xbt <- xyFromCell(r_template, 1:ncell(r_template)) %>% 
  as_tibble() %>% 
  transmute(X = x, Y = y) %>% 
  mutate(Z = predict(xbt_gam, .))
r_interp <- rasterFromXYZ(interp_xbt, crs = crs(xbt))
# plot
plot(mask(r_interp, land, inverse = TRUE))
check_values(xbt$avg_iso_depth, r_interp, p)
foodweb_values <- raster::extract(r_interp, p)
# comparison original vs interpolated values
tibble(Original = getValues(r_xbt),
       Interpolated = getValues(r_interp)) %>% 
  filter(!is.na(Original)) %>% 
  ggplot() +
  aes(Original, Interpolated) +
  geom_bin2d() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_vline(xintercept = max(foodweb_values), linetype = "dotted") +
  geom_smooth() +
  scale_x_continuous(limits = c(0, 200)) +
  scale_fill_gradient(low = "steelblue", high = "tomato") +
  theme_bw()
