library(tidyverse)
library(sf)
source("~/Proj/macro-net/env-net-dist/R/utilities.R")

# network properties
p <- st_read("data/marine.shp") %>% 
  transmute(foodweb.name = fdwb_nm,
            geographic.location = ggrphc_) %>% 
  left_join(read_csv("~/Proj/macro-net/env-net-dist/data/network_metrics.csv")) %>% 
  filter_all(all_vars(!is.na(.))) 
net <- p %>% st_drop_geometry()
dist <- st_distance(p)
ans <- tibble()
pb <- progress::progress_bar$new(total = dim(dist)[1])
for (i in seq_len(nrow(dist))) {
  pb$tick()
  for (j in seq_len(ncol(dist))) {
    if (i < j) {
      x <- net[i, ] %>% 
        tibble() %>% 
        select(where(is.numeric))
      y <- net[j, ] %>% 
        tibble() %>% 
        select(where(is.numeric))
      ans <- bind_rows(ans,
                       tibble(distance = dist[i, j],
                              delta = x - y))
    }
  }
}
ans

scaled_delta <- ans$delta %>% 
  as_tibble() %>% 
  mutate_all(scale)

res <- bind_cols(distance = as.numeric(ans$distance), 
                 scaled_delta)

res %>% 
  pivot_longer(cols = 2:ncol(res), names_to = "property", values_to = "value") %>%
  mutate(value = value[,1]) %>%
  ggplot() +
  aes(distance, value) +
  geom_bin2d(aes(fill = log10(..density..))) +
  geom_smooth() +
  scale_fill_viridis_c() +
  scale_x_log10() +
  xlab("Distance (m)") +
  facet_wrap(~property, scales = "free") +
  theme_bw()
