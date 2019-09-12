library(sf)
library(raster)
library(fasterize)
library(prioritizr)
library(doParallel)
library(tidyverse)
library(lwgeom)
library(here)

load("tmp.RData")



pu <- data.frame(id = 1:ncell(cost),
                 cost = cost[],
                 status = 0L)

spec <- data.frame(id = 1:nlayers(biod),
                   name = names(biod),
                   stringsAsFactors = FALSE)
biod_df <- as.data.frame(biod)
tmp <- as_tibble(data.frame(pu$id, biod_df))
names(tmp)[1] <- "pu"
names(tmp)[-1] <- spec$id

rij_raw <- tmp %>% gather( "species", "amount", -pu)
rij <- rij_raw %>% filter(!is.na(amount) & amount > 0) 

rij$species <- as.integer(rij$species)

rij <- rij %>% arrange(pu, species)


p1 <- problem(pu, cost_column = "cost", features = spec, rij = rij) %>%
  add_min_set_objective() %>%
  add_relative_targets(0.3) %>%
  # add_binary_decisions() %>%
  add_proportion_decisions() %>%
  add_default_solver(gap = 0.001)

# solve problem
s1 <- solve(p1)