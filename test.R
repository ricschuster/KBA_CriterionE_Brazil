library(sf)
library(raster)
library(fasterize)
library(prioritizr)
library(doParallel)
library(tidyverse)
library(lwgeom)
library(here)

load("tmp.RData")

# Code that I ran to produce tmp.RData

#Amph file from here (https://biodiversitymapping.org/wordpress/index.php/home/) to get 10km reference raster
# trm <- raster("Richness_10km_AMPHIBIANS_dec2017_spp_edited_extant_1m_EckertIV_dissolved_Anura_raster.tif")
# base_raster <- raster(trm) 

# #####
# # Mammals
# #####

# source: https://www.iucnredlist.org/resources/spatial-data-download
# setwd("TERRESTRIAL_MAMMALS/")
# 
# mamm <- st_read("TERRESTRIAL_MAMMALS.shp")
# #filter out:
# # presence other than 1 (extant)
# # seasonal 3 (passage) and 4 (seasonal occurrence uncertain)
# # origin 3 (introduced), 4 (vagrant), 5 (origin uncertain)
# mamm_red <- mamm %>% filter(presence == 1, seasonal <= 3, origin <= 2)
# 
# mamm_proj <- mamm_red %>% st_transform(crs = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
# mamm_proj <- mamm_proj %>% st_make_valid()
# mamm_diss <- mamm_proj %>% group_by(binomial) %>% summarise()
# 
# mamm_rast <- raster(mamm_diss, res = 10000)
# 
# for(ii in 1:nrow(mamm_diss)){
#   tmp_rast <- base_raster
#   tmp_rast <- try(fasterize(mamm_diss[ii,], tmp_rast))
#   if(class(tmp_rast) == "RasterLayer" & !all(is.na(tmp_rast[]))){
#     writeRaster(tmp_rast, here("IUCN/Mamm/",paste0(as.character(mamm_diss[ii,]$binomial), ".tif")),
#                 overwrite = TRUE)
#   }
#   rm(tmp_rast)
# }
# source: https://gadm.org/data.html
# land <- readRDS(here("gadm36_BRA_0_sf.rds"))
# land <- st_transform(land, proja)
# 
# base_rast_br <- crop(base_rast, land)
# 
# cost <- fasterize(land, base_rast_br)
# 
# 
# fls <- list.files(path = here("IUCN/Mamm/"), pattern = ".tif$")
# ll <- list()
# jj <- 1
# 
# for(ii in 1:length(fls)){
#   tmp_rast <-  try(raster(fls[ii]) %>% crop(base_rast_br))
#   
#   if(class(tmp_rast) == "RasterLayer" & !all(is.na(values(tmp_rast)))){
#     
#     tmp.df <- data.frame(cost = values(cost),
#                          val = values(tmp_rast))
#     tmp.df.red <- tmp.df %>% filter(!is.na(cost))
#     
#     if(sum(tmp.df.red$val, na.rm = TRUE) >= 10){
#       ll[[jj]] <- tmp_rast
#       jj <- jj + 1
#     }
#     rm(tmp.df, tmp.df.red)
#   }
#   
#   if(!(ii %% 50)){
#     print(ii)
#     flush.console()
#   }
#   rm(tmp_rast)
# }
# 
# biod <- stack(ll)
# 
# save.image("tmp.RData")

####################################################
# Test with rasters
p1 <- problem(cost, biod) %>%
  add_min_set_objective() %>%
  add_relative_targets(0.05) %>%
  add_binary_decisions()# %>%
# add_proportion_decisions() %>%
#add_default_solver(gap = 0.001)

# solve problem
s1 <- solve(p1)
plot(s1)

####################################################
#Test with rij
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

