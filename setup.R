library(sf)
library(raster)
library(fasterize)
library(prioritizr)
library(doParallel)
library(tidyverse)

setwd("D:/Work/KBAs/Crit_E/KBA_CriterionE_Brazil")
library(here)

# parallelization
n_cores <- 12
cl <- makeCluster(n_cores)
registerDoParallel(cl)

grd <- st_read(here("grid/Squares_250km_red_Lines.shp"))
proja <- crs(grd)


land <- readRDS(here("grid/gadm36_BRA_0_sf.rds"))
land <- st_transform(land, crs(grd))

base_rast <- raster(ext = extent(grd), res = 10000, crs = crs(grd))

cost <- fasterize(land, base_rast)

setwd("data/")


fls <- list.files(pattern = ".shp$")
ll <- list()
jj <- 1

for(ii in 1:length(fls)){
  tmp_poly <- st_read(fls[ii])
  tmp_rast <-  fasterize(tmp_poly, base_rast)
  
  if(any(!is.na(tmp_rast[]))){
    ll[[jj]] <- tmp_rast
    jj <- jj + 1
  }
  
  rm(tmp_poly, tmp_rast)
}

biod <- stack(ll)

pu <- data.frame(id = 1:ncell(cost),
                 cost = cost[],
                 status = 0L)

spec <- data.frame(id = 1:nlayers(biod),
                   name = names(biod),
                   stringsAsFactors = FALSE)
biod_df <- as.data.frame(biod)
tmp <- as_tibble(data.frame(pu$id, biod_df))
names(tmp)[1] <- "pu"

rij_raw <- tmp %>% gather( "species", "amount", -pu)
rij <- rij_raw %>% filter(!is.na(amount) & amount > 0) 

rij$species <- as.integer(gsub("layer.", "", rij$species))

rij <- rij %>% arrange(pu, species)


p1 <- problem(pu, cost_column = "cost", features = spec, rij = rij) %>%
  add_min_set_objective() %>%
  add_relative_targets(0.3) %>%
  # add_binary_decisions() %>%
  add_proportion_decisions() %>%
  add_default_solver(gap = 0.001)

# solve problem
s1 <- solve(p1)
# plot(s1)

rc1 <- replacement_cost(p1, data.frame(s1$solution_1), threads = n_cores)

rw1 <- rarity_weighted_richness(p1, data.frame(s1$solution_1))

save.image("D:/Work/KBAs/Crit_E/KBA_CriterionE_Canada/tmp.RData")

rw1.r <- rc1.r <- s1.r <- cost
rc1.r[] <- rc1$rc
rw1.r[] <- rw1$rwr
s1.r[] <- s1$solution_1
plot(s1.r)
plot(rc1.r)
plot(rw1.r)

writeRaster(rw1.r, here("/output/rw1.tif"), format="GTiff", overwrite = TRUE)
writeRaster(rc1.r, here("/output/rc1.tif"), format="GTiff", overwrite = TRUE)
writeRaster(s1.r, here("/output/s1.tif"), format="GTiff", overwrite = TRUE)

# set infinite values as 1.09 so we can plot them
rc$rc[rc$rc > 100] <- 1.09

# plot the irreplaceability scores
# planning units that are replaceable are shown in purple, blue, green, and
# yellow, and planning units that are truly irreplaceable are shown in red
spplot(rc1.r, "layer", main = "Irreplaceability", #xlim = c(-0.1, 1.1),
       #ylim = c(-0.1, 1.1), 
       at = c(seq(0, 0.9, 0.1), 1.01, 1.1),
       col.regions = c("#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
                       "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725",
                       "#FF0000"))
spplot(rw1.r, "layer", main = "Irreplaceability", #xlim = c(-0.1, 1.1),
       #ylim = c(-0.1, 1.1), 
       at = c(seq(0, 0.9, 0.1), 1.01, 1.1),
       col.regions = c("#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
                       "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725",
                       "#FF0000"))

spplot(s1.r, "layer", main = "Irreplaceability", #xlim = c(-0.1, 1.1),
       #ylim = c(-0.1, 1.1), 
       at = c(seq(0, 0.9, 0.1), 1.01, 1.1),
       col.regions = c("#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
                       "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725",
                       "#FF0000"))

# clean up
stopCluster(cl)

