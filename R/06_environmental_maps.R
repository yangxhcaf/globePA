#' ---
#' title: "Global Analysis of Protected Areas - Create environmental data maps"
#' author: "RS-eco"
#' ---

rm(list=ls()); invisible(gc())

# Set working directory
workdir <- "C:/Users/admin/Documents/GitHub/globePA"
setwd(workdir)

#Automatically install required packages, which are not yet installed
`%!in%` <- Negate(`%in%`)
if("scico" %!in% installed.packages()[,"Package"]) remotes::install_github("thomasp85/scico")
if("rnaturalearthhires" %!in% installed.packages()[,"Package"]) remotes::install_github("ropensci/rnaturalearthhires/")
packages <- c("sp", "raster", "tidyverse", "patchwork", "RStoolbox", "scico", "rnaturalearthhires")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load packages
l <- sapply(packages, require, character.only = TRUE, quietly=TRUE); rm(packages, l)

# Specify colour scheme
bluered <- rev(scico(255, palette = 'roma'))
bluewhitered <- rev(scico(255, palette = 'romaO'))

# Obtain world map from rnaturalearthhires package
data(countries10)
countries10 <- sf::st_as_sf(countries10)
outline <- sf::st_union(countries10)
outline_ter <- sf::st_crop(outline, extent(-180, 180, -56, 84))

outline_moll <- sf::st_transform(outline, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
outline_ter_moll <- sf::st_transform(outline_ter, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
rm(countries10, outline_ter, outline); invisible(gc())

########################################

# Plot terrestrial maps raw

# Need to read data directly from Worldclim v2 & EarthEnv folder
# as you have specified in the filedir path at the beginning of the script.

# Create maps

bio1 <- raster::raster("extdata/wc2.0_bio_30s_01_wm.nc")
bio1[bio1 == 0] <- NA
bio1 <- bio1 %>% fortify(maxpixels=2000000) %>% tidyr::drop_na(); invisible(gc())
lim_map <- c(min(bio1$layer), max(bio1$layer))
col_val <- scales::rescale(unique(c(seq(min(bio1$layer), 0, length=128), seq(0, max(bio1$layer), length=128))))
p1 <- ggplot() + geom_tile(data=bio1, aes(x=x,y=y,fill=layer)) + 
  scale_fill_gradientn(name="Temp (°C)", colours = bluewhitered, na.value="transparent",
                       values=col_val, limits=lim_map) + 
  geom_sf(data=outline_ter_moll, fill=NA, color = "black", size = 1/.pt) +
  coord_sf() + theme_bw() + ggtitle("(a)") + 
  theme(axis.title=element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        panel.border = element_blank(), plot.title = element_text(vjust = - 7, hjust=0.07))
rm(bio1); invisible(gc())

bio12 <- raster::raster("extdata/wc2.0_bio_30s_12_wm.nc")
bio12[bio12 == 0] <- NA
bio12 <- bio12 %>% fortify(maxpixels=2000000) %>% tidyr::drop_na(); invisible(gc())
#lim_map <- c(0, max(bio12$layer))
#col_val <- scales::rescale(c(0, quantile(bio12$layer, probs=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), max(bio12$layer)))
p2 <- ggplot() + geom_tile(data=bio12, aes(x=x,y=y,fill=layer)) + 
  scale_fill_gradientn(name="Prec (mm)", colours = bluered, na.value="transparent", limits=c(0, NA),
                       values=c(0, 0.025, 0.05, 0.1, 0.15, 0.25, 0.5, 1)) + 
  geom_sf(data=outline_ter_moll, fill=NA, color = "black", size = 1/.pt) +
  coord_sf() + theme_bw() + ggtitle("(b)") + 
  theme(axis.title=element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        panel.border = element_blank(), plot.title = element_text(vjust = - 7, hjust=0.07))
rm(bio12); invisible(gc())

elevation <-  raster::raster("extdata/elevation_1KMmn_GMTEDmn_wm.nc")
elevation[elevation == 0] <- NA
elevation <- elevation %>% fortify(maxpixels=2000000) %>% tidyr::drop_na(); invisible(gc())
p3 <- ggplot() + geom_tile(data=elevation, aes(x=x,y=y,fill=layer)) + 
  scale_fill_gradientn(name="Elevation (m)", colours = bluered, na.value="transparent",
                       values=c(0, 0.1, 0.2, 0.3, 0.6, 0.8, 1)) + 
  geom_sf(data=outline_ter_moll, fill=NA, color = "black", size = 1/.pt) +
  coord_sf() + theme_bw() + ggtitle("(c)") + 
  theme(axis.title=element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        panel.border = element_blank(), plot.title = element_text(vjust = - 7, hjust=0.07))
rm(elevation); invisible(gc())

# Join figures
p <- p1 + p2 + p3 + plot_layout(ncol=1)
ggsave(filename="figures/ter_maps_raw.png", p, width=5.8, height=7, dpi=600)

########################################

# Plot marine maps raw

# Need to read data directly from MARSPEC files!!!

# Create maps

biogeo13 <- raster::raster("extdata/biogeo13_30s_wm.nc") %>% fortify(maxpixels=2000000) %>% tidyr::drop_na()
head(biogeo13); tail(biogeo13)
p1 <- ggplot() + geom_tile(data=biogeo13, aes(x=x,y=y,fill=biogeo13_30s_wm)) + 
  scale_fill_gradientn(name="SST (°C)", colours = bluered, na.value="transparent",
                       breaks=c(0,500,1000,1500,2000,2500), labels=c(0, 5, 10, 15, 20, 25)) + 
  geom_sf(data=outline_moll, fill=NA, color = "black", size = 1/.pt) +
  coord_sf() + theme_bw() + ggtitle("(a)") + 
  theme(axis.title=element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        panel.border = element_blank(), plot.title = element_text(vjust = - 7, hjust=0.07))
rm(biogeo13); invisible(gc())

biogeo08 <- raster::raster("extdata/biogeo08_30s_wm.nc") %>% fortify(maxpixels=2000000) %>% tidyr::drop_na()
p2 <- ggplot() + geom_tile(data=biogeo08, aes(x=x,y=y,fill=biogeo08_30s_wm)) + 
  scale_fill_gradientn(name="SSS (psu)", colours = bluered, na.value="transparent",
                       values=c(0, 0.25, 0.5, 0.65, 0.75, 0.85, 0.9, 0.95, 1),
                       breaks=c(1000, 2000, 3000, 4000), labels=c(10, 20, 30, 40)) + 
  geom_sf(data=outline_moll, fill=NA, color = "black", size = 1/.pt) +
  coord_sf() + theme_bw() + ggtitle("(b)") + 
  theme(axis.title=element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        panel.border = element_blank(), plot.title = element_text(vjust = - 7, hjust=0.07))
rm(biogeo08); invisible(gc)

bathy <- raster::raster("extdata/bathy_30s_wm.nc") %>% fortify(maxpixels=2000000) %>% tidyr::drop_na()
p3 <- ggplot() + geom_tile(data=bathy, aes(x=x,y=y,fill=bathy_30s_wm)) + 
  scale_fill_gradientn(name="Bathy (m)", colours = rev(bluered), 
                       values=c(0, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.9, 0.95, 1),
                       na.value="transparent") + 
  geom_sf(data=outline_moll, fill=NA, color = "black", size = 1/.pt) +
  coord_sf() + theme_bw() + ggtitle("(c)") + 
  theme(axis.title=element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        panel.border = element_blank(), plot.title = element_text(vjust = - 7, hjust=0.07))
rm(bathy); invisible(gc())

# Join figures
p <- p1 + p2 + p3 + plot_layout(ncol=1)
ggsave(filename="figures/mar_maps_raw.png", p, width=5.8, height=8, dpi=600)
