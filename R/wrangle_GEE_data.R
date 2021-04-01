# ####### Summary ########
#
# I combined together all data from GEE and ran a PCA to select variables that
# were not collinear (VIF < 5).
#
# I inspected the selected variables and decided to log10-transform
# avg_chlorophyll-a.
#
# I scaled and centered all variables.
library(tidyverse)
source("~/Proj/macro-net/env-net-dist/R/utilities.R") #for vif()

# combine all GEE files
f <- list.files("data/GEE/GEE-marine", pattern = "csv", full.names = TRUE)
vars <- (str_split(f, "/", simplify = TRUE)[, 4] %>% 
  str_split(., "_", simplify = TRUE))[, 1] %>% 
  unique()
d <- lapply(f, function(x) {
  name <- str_split(x, "/", simplify = TRUE)[4]
  var <- str_split(name, "_", simplify = TRUE)[1]
  stat <- str_split(name, "_", simplify = TRUE)[2]
  year <- str_split(name, "_", simplify = TRUE)[3] %>% 
    gsub("[.]csv", "", .)
  read_csv(x, col_types = cols()) %>% 
    tibble() %>% 
    transmute(foodweb.name = fdwb_nm,
              geographic.location = ggrphc_,
              val = mean,
              stat = stat,
              var = var,
              year = year)
}) %>% bind_rows()
# summarize all years into one value
res <- d %>% 
  group_by(var, stat, foodweb.name, geographic.location) %>% 
  summarize(avg = mean(val, na.rm = TRUE),
            interannual_sd = sd(val, na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(var, foodweb.name, avg, stat, interannual_sd) %>% 
  distinct_all()
# plots to inspect data
res %>% 
  filter(stat == "mean") %>% 
  ggplot() +
  geom_point(aes(avg, interannual_sd),
                  alpha = 0.5) +
  xlab("Average mean") +
  ylab("Annual variation") +
  theme_bw() +
  facet_wrap(~var, scales = "free")
res %>% 
  filter(stat == "sd") %>% 
  ggplot() +
  geom_point(aes(avg, interannual_sd)) +
  xlab("Average seasonality") +
  ylab("Annual variation") +
  theme_bw() +
  facet_wrap(~var, scales = "free")
# missing values
res %>% 
  filter_all(any_vars(is.na(.)))
# PCA to decide which variable retain
res <- res %>% 
  pivot_wider(names_from = c(var, stat), values_from = c(avg, interannual_sd))
pca <- res %>% 
  select(where(is.numeric)) %>% 
  filter_all(all_vars(!is.na(.))) %>% 
  FactoMineR::PCA(scale.unit = TRUE)
# select variables
res %>% 
  select(where(is.numeric)) %>% 
  select(`avg_water-temp-0_mean`,
         `avg_water-temp-0_sd`,
         `interannual_sd_water-temp-0_mean`,
         `avg_chlor-a_mean`,
         `interannual_sd_chlor-a_sd`,
         `interannual_sd_salinity-0_mean`) %>% 
  filter_all(all_vars(!is.na(.))) %>% 
  vif() %>% 
  tibble() %>% 
  mutate(VIF = as.numeric(VIF),
         R2 = round(as.numeric(R2), 2)) %>% 
  arrange(desc(VIF))
# avg_chl log10 transformed
res <- res %>% 
  transmute(foodweb.name,
            avg_water_temp = `avg_water-temp-0_mean`,
            seasonality_water_temp = `avg_water-temp-0_sd`,
            avg_chl = log10(`avg_chlor-a_mean`),
            interann_seasonality_chl = `interannual_sd_chlor-a_sd`,
            interann_avg_salininty = `interannual_sd_salinity-0_mean`,
            interann_avg_temp = `interannual_sd_salinity-0_mean`)
# scale and center
res_scaled <- res %>% 
  mutate_if(is.numeric, function(x) as.vector(scale(x)))
# inspection plot
res_scaled %>% 
  pivot_longer(cols = 2:5, names_to = "var", values_to = "val") %>% 
  ggplot() +
  geom_boxplot(aes(var, val)) +
  theme_bw()
# save to file
write_csv(res, "data/env.csv")
