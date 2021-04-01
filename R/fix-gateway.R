library(tidyverse)
library(magrittr)

gw_marine <- read_csv("~/Proj/gateway/gateway.csv", guess_max = 200000) %>%
  filter(ecosystem.type == "marine")

gw_marine %<>% 
  mutate(
    sampling.start.year = ifelse(grepl("Mendonca", link.citation),
                                 2013,
                                 sampling.start.year),
    sampling.end.year = ifelse(grepl("Mendonca", link.citation),
                               2015,
                               sampling.start.year)
  ) %>%
  # Kefi constists of three source:
  # - one from 1995 to 2005 (Blanhcette et al. 2007)
  # - one from 1997 to 1998 (Boritman et al. 2001)
  # - one from 1998 - 2004 (Wieters 2005)
  # I assign the longest time range here, i.e. 1995-2005.
  mutate(
    sampling.start.year = ifelse(grepl("Kefi", link.citation),
                                 1995,
                                 sampling.start.year),
    sampling.end.year = ifelse(grepl("Kefi", link.citation),
                               2005,
                               sampling.end.year)
  ) %>% 
  # Laffery et al. (2006) is a metaweb
  mutate(
    sampling.start.year = ifelse(grepl("Lafferty", link.citation),
                                 2006,
                                 sampling.start.year),
    sampling.end.year = ifelse(grepl("Lafferty", link.citation),
                               2006,
                               sampling.end.year)
  ) %>% 
  mutate(
    sampling.start.year = ifelse(link.citation == "Jacob et al. (2011)",
                                 2001,
                                 sampling.start.year),
    sampling.end.year = ifelse(link.citation == "Jacob et al. (2011)",
                               2004,
                               sampling.end.year)
  ) %>% 
  # Jacob et al. (2015) is a metaweb
  mutate(
    sampling.start.year = ifelse(link.citation == "Jacob et al. (2015)",
                                 2015,
                                 sampling.start.year),
    sampling.end.year = ifelse(link.citation == "Jacob et al. (2015)",
                               2015,
                               sampling.end.year)
  )


gw_marine %>% 
  group_by(link.citation) %>% 
  summarize(Tmin = min(sampling.start.year),
            Tmax = max(sampling.end.year))

gw_marine %>% 
  filter(link.citation == "Lafferty et al. (2006)") %>% 
  group_by(foodweb.name, study.site) %>% 
  tally()
