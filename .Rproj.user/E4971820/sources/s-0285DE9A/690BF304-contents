library(tidyverse)
library(adespatial)
library(ade4)
library(adegraphics)
library(sf)
library(spdep)
library(vegan)
source("~/Proj/macro-net/env-net-dist/R/utilities.R")
source("~/Proj/R_misc/monkey_theme.R")
source("R/find_mems.R")

# network properties
p <- st_read("data/marine.shp") %>% 
  transmute(foodweb.name = fdwb_nm,
            geographic.location = ggrphc_) %>% 
  left_join(read_csv("~/Proj/macro-net/env-net-dist/data/network_metrics.csv")) %>% 
  filter_all(all_vars(!is.na(.)))
net <- p %>% 
  st_drop_geometry() %>% 
  as_tibble()
# environmental variables
env <- read_csv("data/env.csv")
# Moran's Eigenvector Maps
if ("mems.csv" %in% list.files("data")) {
  mems <- read_csv("data/mems.csv")
} else {
  coords <- st_coordinates(p) %>% 
    as_tibble() %>% 
    transmute(lon = X, lat = Y)
  d <- p %>% 
    st_drop_geometry() %>% 
    bind_cols(coords)
  best_mems <- find_mems(d,
                         weights = "fup", "fdown", "flin",
                         nb = c("pcnm", "gab", "dnear", "del"))
  distance <- st_distance(p)
  # plot SWM in geographic space
  spp <- as_Spatial(p)
  plot(spp, pch = 20, col = adjustcolor("red", alpha.f = 0.5), cex = 1)
  lines(spdep::listw2lines(best_mems$lw.best, spp),
        col = adjustcolor("blue", alpha.f = 0.2),
        cex = 1)
  SWM <- listw2mat(best_mems$lw.best)
  dist <- st_distance(p)
  tibble(Distance = as.vector(dist),
         Weight = as.vector(SWM)) %>% 
    filter(Distance > 0) %>% 
    ggplot() +
    geom_point(aes(Distance, Weight)) +
    scale_x_log10(labels = monkey_scientific) +
    xlab("Distance (km)") +
    theme_bw()
  write_rds(best_mems, "analysis-output/mems.rds")
  # save MEMs of sites to file
  mems <- tibble(foodweb.name = p$foodweb.name) %>% 
    bind_cols(best_mems$best.method$MEM.select)
  write_csv(mems, "data/mems.csv")
}

d <- left_join(net, env) %>% 
  left_join(mems)
d <- d %>% 
  filter_all(all_vars(!is.na(.)))
response <- d[, colnames(net)[3:16]] %>% 
  as_tibble()
predictors <- d[, c(colnames(env)[2:ncol(env)],
                    colnames(mems)[2:ncol(mems)])] %>% 
  as_tibble()
# PCA ---------
pca <- dudi.pca(response,
                scale = TRUE,
                scannf = FALSE,
                nf = 3)
# RDA ----------
rda <- pcaiv(pca,
             predictors,
             scannf = FALSE,
             nf = 2)
# plot rda ----------
bind_rows(
  rda$co %>% 
    tibble() %>% 
    mutate(var = rownames(rda$co),
           type = "food web property"),
  rda$cor[!grepl("MEM", rownames(rda$cor)), ] %>% 
    tibble() %>% 
    mutate(var = rownames(rda$cor)[!grepl("MEM", rownames(rda$cor))],
           type = "environmental predictor") %>% 
    rename_with(~gsub("RS", "Comp", .x), starts_with("RS")),
  rda$cor[grepl("MEM", rownames(rda$cor)), ] %>% 
    tibble() %>% 
    mutate(var = rownames(rda$cor)[grepl("MEM", rownames(rda$cor))],
           type = "MEM") %>% 
    rename_with(~gsub("RS", "Comp", .x), starts_with("RS"))
) %>%
  ggplot() +
  geom_segment(aes(x = 0, xend = Comp1,
                   y = 0, yend = Comp2,
                   col = type),
               arrow = arrow(length = unit(0.03, "npc"))) +
  geom_text(aes(x = Comp1, y = Comp2 * 1.1, label = var)) +
  scale_color_manual(values = c("tomato", "steelblue", "forestgreen")) +
  theme_minimal()
# correlation --------
angles <- bind_rows(
  rda$co %>% 
    tibble() %>% 
    mutate(var = rownames(rda$co),
           type = "response"),
  rda$cor %>% 
    tibble() %>% 
    mutate(var = rownames(rda$cor),
           type = "mem") %>% 
    rename_with(~gsub("RS", "Comp", .x), starts_with("RS"))
) %>% 
  mutate(angle = modify2(Comp1, Comp2, function(x, y) get_angle(x, y))) %>%
  select(var, angle)
correlation <- matrix(NA, nrow = nrow(angles), ncol = nrow(angles))
rownames(correlation) <- angles$var
colnames(correlation) <- angles$var
for (i in seq_len(nrow(correlation))) {
  for (j in seq_len(ncol(correlation))) {
    vals <- angles %>% 
      filter(var %in% c(colnames(correlation)[j], rownames(correlation)[i])) %>% 
      pull(angle)
    if (length(vals) == 1) {
      vals <- rep(vals, 2)
    }
    correlation[i, j] <- cos(diff(vals))
  }
}
# correlatiom plot --------
non_mems <- rownames(correlation)[!grepl("MEM", rownames(correlation))]
corrplot::corrplot.mixed(correlation[non_mems, non_mems],
                         diag = NULL,
                         lower = "number",
                         bg = "grey20",
                         tl.pos = "lt",
                         tl.col = c(rep("steelblue", 14),
                                    rep("tomato", 6),
                                    rep("forestgreen", 14)),
                         number.cex = 0.8)
randtest(rda)
# variation partitioning ----------
vp <- vegan::varpart(pca$tab, predictors[, 1:6], predictors[, 7:18])
plot(vp, bg = c(3, 5), Xnames = c("Environment", "Spatial"))
# scalograms ---------
spp <- right_join(p, d) %>% 
  mutate(X = st_coordinates(.)[, 1],
         Y = st_coordinates(.)[, 2]) %>%
  as_tibble() %>%
  st_as_sf(coords = c("X", "Y"), crs = "EPSG:4326") %>%
  as_Spatial()
best_mems$lw.best
neighs <- dnearneigh(spp, d1 = 0, d2 = 7900) #d2 close enough to best
neighs
listw <- nb2listw(neighs, zero.policy = TRUE, style = "W")
mem_dnear <- mem(listw)
cairo_pdf("figures/scalograms.pdf",
          width = 6,
          height = 4,
          onefile = TRUE)
for (x in names(spp)[3:22]) {
  v <- spp[[x]]
  scalo <- scalogram(v, mem_dnear, nblocks = 116)
  plot(scalo[which(scalo$pvalue <= 0.5)], sub = x)
}
dev.off()
plot(scalogram(spp$E, mem_dnear, nblocks = 20), sub = "S")


mspa_dudi <- mspa(pca, listw, scannf = FALSE, nf = 2)
scatter(mspa_dudi, plot = TRUE)