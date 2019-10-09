library(raster)
library(readxl)
library(dplyr)
library(plyr)
library(tidyverse)
library(sf)
library(raster)
library(RColorBrewer)
library(here)
library(leaflet)
library(png)
library(fields)
library(rgdal)
library(stringr)
library(cluster)
library(smoothr)
library(ggpubr)
library(ggspatial)

load("portfolio_irreplaceability_Targets_Marxan.RData")

ne_land <- read_sf("data/ne-land.gpkg") %>% st_geometry()
ne_country_lines <- read_sf("data/ne-country-lines.gpkg") %>% st_geometry()
ne_state_lines <- read_sf("data/ne-state-lines.gpkg") %>% st_geometry()

stem_crop <- function(x) {
  stopifnot(inherits(x, "Raster"))
  
  # aggregate for faster processing
  x_agg <- raster::aggregate(x, fact = 3)
  
  # extent of non-NA
  x_agg <- stem_to_na(x_agg)
  #x_agg <- raster::trim(x_agg, values = NA)
  #x_ext <- raster::extent(x_agg)
  x_ext <- extent_na(x_agg)
  raster::crop(x, x_ext)
}

extent_na <- function(x) {
  pts <- raster::rasterToPoints(x)
  x_rng <- range(pts[, "x"])
  y_rng <- range(pts[, "y"])
  raster::extent(x_rng[1] - res(x)[1] / 2, x_rng[2] + res(x)[1] / 2,
                 y_rng[1] - res(x)[2] / 2, y_rng[2] + res(x)[2] / 2)
}

stem_to_na <- function(x, value = 0) {
  stopifnot(inherits(x, "Raster"))
  stopifnot(is.numeric(value), length(value) == 1)
  
  if (inherits(x, "RasterLayer")) {
    x[x[] == value] <- NA_real_
  } else {
    for (i in seq.int(raster::nlayers(x))) {
      x[[i]][x[[i]][] == value] <- NA_real_
    }
  }
  return(x)
}

proper <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
}

rrT2 <- stack(tmp_r, rst, rt3, rt4,r.rc.t1, r.rw.t1)
names(rrT2) <- c("Marxan", 
                 "prioritizr", 
                 "pool portfolio", 
                 "shuffle portfolio",
                 "raplecement cost", 
                 "rarity weighted richness")

crs <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
abd_plot3 <- rrT2 %>% 
  stem_to_na() %>% 
  projectRaster(crs = crs, method = 'ngb') %>% 
  #sqrt() %>% 
  stem_crop()

e <- extent(abd_plot3)
#e <- extent(abd_plot)
text_col <- "black"

#weeks <- nms
palette <- c("Greens", "Blues", "YlOrRd", "Reds")
legend_offsets <- c(0.01, 0.06, 0.11, 0.16)
# prepare vector layers
land <- st_transform(ne_land, crs = proj4string(abd_plot3))
country <- st_transform(ne_country_lines, crs = proj4string(abd_plot3))
state <- st_transform(ne_state_lines, crs = proj4string(abd_plot3))
#logo <- readPNG("STEM-logos-lab-ebird-3000.png")

#names(abd_plot3) <- c("abundance", "clusters")
# plot weeks
add_legend4 <- function(title, palette, bump = 0, low_high = FALSE, 
                        text_col = "black") {
  if (low_high) {
    labs <- list(at = c(0, 1), labels = c("low", "high"), line = -1,
                 cex.axis = 1, fg = NA, col.axis = text_col)
    
  } else {
    labs <- list(at = c(0, 1), labels = NA, line = 0, fg = NA)
  }
  fields::image.plot(zlim = c(0, 1), legend.only = TRUE, col = palette(256),
                     legend.width = 1, horizontal = TRUE,
                     smallplot = c(0.70, 0.90, 0.05 + bump, 0.075 + bump),
                     axis.args = labs,
                     legend.args = list(text = proper(title), side = 1,
                                        col = text_col))
}


##Plotting
# here("output", paste0("SNGO_04-25_5cl_50samp_size_20perc", ".png")) %>% 
#   png(width = 3000, height = 3000, res = 300)
# 
# par(mfrow=c(2,2))
rrT2_df <- as.data.frame(rrT2)
rrT2_df$prioritizr <- rrT2_df$prioritizr * 100
rrT2_df$raplecement.cost <- rrT2_df$raplecement.cost * 100
rrT2_df$rarity.weighted.richness <- rrT2_df$rarity.weighted.richness * 100
rr_sel <- sapply(rrT2_df, function(x) sum(x > 90, na.rm = TRUE) )


#A)
show_legend <- c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
legend_head <- c("Selection frequency",
                 "",
                 "Selection frequency",
                 "Selection frequency",
                 "Irreplaceability",
                 "Irreplaceability"
)

for(ii in 1:nlayers(abd_plot3)){
  here::here("figures", paste0(names(abd_plot3[[ii]]), ".png")) %>% 
    png(width = 3000, height = 3000, res = 300)
  
  pal <- brewer.pal(9, palette[3])[2:9]
  pal <- colorRampPalette(pal)
  par(mar = c(0.1, 0.1, 0.1, 0.1), oma = c(0,0,0,0), bg = "white")
  plot(land, col = "grey85", border = NA, xlim = e[1:2], ylim = e[3:4])
  
  plot(abd_plot3[[ii]], add = TRUE, col = pal(256), legend = FALSE, 
       maxpixels = ncell(abd_plot3))
  
  if(show_legend[ii]){
    add_legend4("", pal, 0.10, low_high = TRUE,
                text_col = text_col)
  }
  
  # boundaries
  plot(state, col = "black", lwd = 0.5, lty = 1, add = TRUE)
  plot(country, col = "black", lwd = 1, add = TRUE)
  
  # title
  # plot bounds
  usr <- par("usr")
  xwidth <- usr[2] - usr[1]
  yheight <- usr[4] - usr[3]
  # labels
  
  text(x = usr[1] + 0.70 * xwidth, y = usr[3] + 0.20 * yheight,
       labels = legend_head[ii], pos = 4, font = 1, cex = 1.2, col = text_col)
  
  
  text(x = usr[1] + 0.01 * xwidth, y = usr[3] + 0.97 * yheight,
       labels = names(abd_plot3)[ii], pos = 4, font = 1, cex = 3, col = text_col)
  
  text(x = usr[1] + 0.01 * xwidth, y = usr[3] + 0.92 * yheight,
       labels = paste0(">90%: ", as.numeric(rr_sel[ii])), pos = 4, font = 1, cex = 3, col = text_col)
  
  dev.off()
  
}

