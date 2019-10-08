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


rrT2 <- stack(tmp_r, st, sum(stack(st3)), sum(stack(st4)), rc.t1, rw.t1)
names(rrT2) <- c("Marxan", 
                 "prioritizr", 
                 "pool portfolio", 
                 "shuffle portfolio",
                 "raplecement cost", 
                 "rarity weighted richness")

crs <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
abd_plot3 <- rrT2 %>% 
  stem_to_na() %>% 
  # projectRaster(crs = crs) %>% 
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
  
  
  text(x = usr[1] + 0.05 * xwidth, y = usr[3] + 0.95 * yheight,
       labels = names(abd_plot3)[ii], pos = 4, font = 1, cex = 2, col = text_col)
  
  dev.off()
  
}

#B)
pal <- brewer.pal(9, "YlOrRd")[1:9]
pal <- colorRampPalette(pal)
plot(land, col = "grey85", border = NA, xlim = e[1:2], ylim = e[3:4])
plot(abd_plot3[[3]], add = TRUE, col = pal(256), legend = FALSE, 
     maxpixels = ncell(abd_plot3))

# boundaries
plot(state, col = "black", lwd = 0.5, lty = 1, add = TRUE)
plot(country, col = "black", lwd = 1, add = TRUE)

text(x = usr[1] + 0.6 * xwidth, y = usr[3] + 0.2 * yheight,
     labels = names(abd_plot3)[3], pos = 4, font = 1, cex = 1.2, col = text_col)

#C)
pal <- brewer.pal(5, "Set2")[1:5]
pal <- colorRampPalette(pal)
plot(land, col = "grey85", border = NA, xlim = e[1:2], ylim = e[3:4])
plot(abd_plot3[[4]], add = TRUE, col = pal(256), legend = FALSE, 
     maxpixels = ncell(abd_plot3))

# boundaries
plot(state, col = "black", lwd = 0.5, lty = 1, add = TRUE)
plot(country, col = "black", lwd = 1, add = TRUE)

text(x = usr[1] + 0.6 * xwidth, y = usr[3] + 0.2 * yheight,
     labels = names(abd_plot3)[4], pos = 4, font = 1, cex = 1.2, col = text_col)


#Abundance

pal <- brewer.pal(9, palette[3])[2:9]
pal <- colorRampPalette(pal)

# modified plasma palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)
pal <- colorRampPalette(plasma(10, direction = -1))

#pal <- colorQuantile(pal, values(abd_plot[[weeks[ii]]]), n = 8, #probs = seq(0, 1, length.out = n + 1),
#              na.color = "#808080", alpha = FALSE, reverse = FALSE)
plot(land, col = "grey85", border = NA, xlim = e[1:2], ylim = e[3:4])

plot(abd_plot3[[1]], add = TRUE, col = pal(256), legend = FALSE, 
     maxpixels = ncell(abd_plot3))

# boundaries
plot(state, col = "black", lwd = 0.5, lty = 1, add = TRUE)
plot(country, col = "black", lwd = 1, add = TRUE)

text(x = usr[1] + 0.6 * xwidth, y = usr[3] + 0.2 * yheight,
     labels = names(abd_plot3)[1], pos = 4, font = 1, cex = 1.2, col = text_col)


dev.off()

par(mfrow=c(1,1))






##Variables
here("output", paste0("SNGO_vars", ".png")) %>% 
  png(width = 3000, height = 3000, res = 300)

par(mfrow=c(3,3))

pal <- brewer.pal(9, palette[3])[2:9]
pal <- colorRampPalette(pal)

# modified plasma palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)
pal <- colorRampPalette(plasma(10, direction = -1))

par(mar = c(0.1, 0.1, 0.1, 0.1), oma = c(0,0,0,0), bg = "white")

for(ii in 3:ncol(out.df.df3)){
  
  tmp.r <-mean_stack
  tmp.r[] <- out.df.df3[,ii]
  crs <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
  tmp.r <- tmp.r %>% 
    stem_to_na() %>% 
    projectRaster(crs = crs) %>% 
    #sqrt() %>% 
    stem_crop()
  
  plot(land, col = "grey85", border = NA, xlim = e[1:2], ylim = e[3:4])
  
  # plot(abd_plot3[[1]], add = TRUE, col = pal(256), legend = FALSE, 
  #      maxpixels = ncell(abd_plot3))
  plot(tmp.r, add = TRUE, col = pal(256), legend = FALSE, 
       maxpixels = ncell(tmp.r))
  
  # boundaries
  plot(state, col = "black", lwd = 0.5, lty = 1, add = TRUE)
  plot(country, col = "black", lwd = 1, add = TRUE)
  
  # title
  # plot bounds
  usr <- par("usr")
  xwidth <- usr[2] - usr[1]
  yheight <- usr[4] - usr[3]
  # labels
  
  text(x = usr[1] + 0.8 * xwidth, y = usr[3] + 0.2 * yheight,
       labels = names(out.df.df3)[ii], pos = 4, font = 1, cex = 1.2, col = text_col)
  
}

dev.off()

par(mfrow=c(1,1))




#### FIGURE 1
pll <- c("a)", "b)", "c)", "d)")

f1 <- c("WK_NC_NH", "YR_NC_NH", "WK_NC_HF", "YR_NC_HF")

here("Figures", paste0("Figure1_zoomed", ".png")) %>% 
  png(width = 3000, height = 3000, res = 300)

par(mfrow=c(2,2))
par(mar = c(0.1, 0.1, 0.1, 0.1), oma = c(0,0,0,0), bg = "white")

# plot f1
for (ii in seq_along(f1)) {
  # print map
  plot(land, col = "grey85", border = NA, xlim = e[1:2], ylim = e[3:4])
  
  pal <- brewer.pal(9, palette[3])[2:9]
  pal <- colorRampPalette(pal)

  # boundaries
  plot(state, col = "black", lwd = 0.5, lty = 1, add = TRUE)
  plot(country, col = "black", lwd = 1, add = TRUE)
  
  # title
  # plot bounds
  usr <- par("usr")
  xwidth <- usr[2] - usr[1]
  yheight <- usr[4] - usr[3]
  # labels
  text(x = usr[1] + 0.65 * xwidth, y = usr[3] + 0.9 * yheight,
       labels = title[ii], pos = 4, font = 1, cex = 1.2 * scl, col = text_col)
  
  text(x = usr[1] + 0.01 * xwidth, y = usr[3] + 0.4 * yheight,
       labels = pll[ii], pos = 4, font = 1, cex = 1.5 * scl, col = text_col)
  
  #rasterImage(logo,usr[1] + 0.01 * xwidth, usr[3] + 0.03 * yheight,
  #            usr[1] + 0.38 * xwidth, usr[3] + 0.09 * yheight)
  
}

dev.off()

