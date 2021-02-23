#' ---
#' title: "Global Analysis of Protected Areas - Create figures for method scheme"
#' author: "RS-eco"
#' ---

# Set working directory
workdir <- "C:/Users/admin/Documents/GitHub/globePA"
setwd(workdir)

# Set file directory
filedir <- "F:/Data/"

#Automatically install required packages, which are not yet installed
packages <- c("sf", "raster", "tidyverse", "fasterize", "RStoolbox", "scico")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load packages
l <- sapply(packages, require, character.only = TRUE, quietly=TRUE); rm(packages, l)

# Specify colour scheme
bluered <- rev(scico(255, palette = 'roma'))

########################################

#download.file("https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_Jan2020_DEU-shapefile.zip",
#              destfile=paste0(filedir, "WDPA/WDPA_Jan2020_DEU-shapefile.zip"))
pa_poly <- sf::st_read(dsn=paste0(filedir, "/WDPA/WDPA_Jan2020_DEU-shapefile-polygons.shp"))
pa_bg <- pa_poly[grep(pa_poly$NAME, pattern="Berchtesgaden"),]

plot(st_geometry(pa_bg[4,]))
pa_sub <- pa_bg[4,]

pa_sub_moll <- sf::st_transform(pa_sub, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
plot(st_geometry(pa_sub_moll))

r <- raster::raster(nrow=14969, ncol=44012, xmn=-18045016, xmx=18044824, ymn=-6492085, ymx=8776295,
                    crs=" +proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
r_sub <- raster::crop(r, pa_sub_moll)
r_sub[] <- 1

r_sub_high <- raster::disaggregate(r_sub, fact=10)
r_pa_sub_high <- fasterize::fasterize(pa_sub_moll, r_sub_high)
r_pa_sub <- raster::aggregate(r_pa_sub_high, 10, fun=sum)
plot(r_pa_sub)

df_pa_sub <- as.data.frame(raster::rasterToPoints(r_pa_sub)) 
df_pa_all <- tidyr::expand(df_pa_sub, x, y) %>% left_join(df_pa_sub) %>% replace_na(list(layer = 0))

ggplot() + 
  geom_sf(data=pa_sub_moll, fill=NA, colour="black") + 
  geom_tile(data=df_pa_all, aes(x, y), fill="transparent", colour="black") + 
  coord_sf(ndiscr=F) + theme_minimal() + 
  theme(axis.title=element_blank(), axis.text = element_blank())
ggsave(filename="figures/pa_example_shape.png", dpi=1000, height=3.5, width=4)

########################################

ggplot() + 
  geom_tile(data=df_pa_all, aes(x, y, fill = layer)) + 
  geom_sf(data = pa_sub_moll, fill=NA, colour="black") + 
  geom_tile(data=df_pa_all, aes(x, y), fill="transparent", colour="black") + 
  scale_x_continuous() + scale_y_continuous() + 
  scale_fill_gradientn(name="Protection (%)", colors=bluered) + 
  theme_minimal() + coord_sf(ndiscr=F) + 
  theme(axis.title=element_blank(), axis.text = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
ggsave(filename="figures/pa_example_coverage.png", dpi=1000, height=3, width=5.2)

########################################

# Plot terrestrial maps raw of example PA

# Need to read data directly from Worldclim v2 & EarthEnv files!!!

# Create maps
bio1 <- raster::raster("extdata/wc2.0_bio_30s_01_wm.nc")
bio1_sub <- raster::crop(bio1, pa_sub_moll, snap="out")

# Plot map
ggplot() + ggR(bio1_sub, geom_raster=T, ggLayer=T) + 
  scale_fill_gradientn(name="Temp (°C)", colours = bluered, na.value="transparent",
                       values=c(0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) + 
  geom_sf(data = pa_sub_moll, fill=NA, colour="black") + 
  #geom_tile(data=df_pa_all, aes(x, y), fill="transparent", colour="black") + 
  theme_minimal() + coord_sf(ndiscr=F) + 
  theme(axis.title=element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
ggsave(paste0("figures/pa_example_bio1.png"), height=3, width=5, dpi=1000)

########################################

# Plot terrestrial maps cut for example PA

# Create maps
bio1 <- raster::raster("extdata/wc2.0_bio_30s_01_perc.nc")
bio1_sub <- raster::crop(bio1, pa_sub_moll, snap="out")

m_ee <- readRDS("data/summary_wc_perc_optim.rds")

# Round first and last value to include all ranges!
m_ee[1,2] <- floor(m_ee[1,2])
m_ee[102,2] <- floor(m_ee[102,2])
m_ee[203,2] <- floor(m_ee[203,2])
m <- m_ee %>% group_split(var)

bio1_df <- as.data.frame(m[[1]]) %>% 
  rowwise() %>% mutate(val=paste0("(", round(x, 2), ",", round(y, 2), "]")) %>%
  select(-c(x,y,var)) %>% as.data.frame()

# Subsitute values
bio1_sub <- raster::subs(bio1_sub, bio1_df, by="z")

# Plot map
# Specify colour scheme
bluered2 <- rev(scico(17, palette = 'roma'))

ggplot() + ggR(bio1_sub, geom_raster=T, ggLayer=T) + 
  geom_sf(data = pa_sub_moll, fill=NA, colour="black") + 
  geom_tile(data=df_pa_all, aes(x, y), fill="transparent", colour="black") + 
  scale_fill_manual(name="Temp (°C)", values=bluered2, na.value="transparent", na.translate=F) + 
  scale_x_continuous() + scale_y_continuous() + 
  theme_minimal() + coord_sf(ndiscr=F) + 
  theme(axis.title=element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), legend.key.height = unit(0.4, "cm"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=10))
ggsave(paste0("figures/pa_example_bio1_perc.png"), height=3, width=5, dpi=1000)

########################################

# Create table

perc_prot <- raster::crop(raster::raster("extdata/pa_cov_total.nc"), pa_sub_moll)
perc_prot <- raster::resample(perc_prot, bio1_sub)
bio1_dat <- stack(bio1_sub, perc_prot)

df_bio1_sub <- as.data.frame(raster::rasterToPoints(bio1_dat))
df_bio1_sub <- df_bio1_sub %>% left_join(bio1_df, by=c("val"="z"))
df_bio1_sub[99:110,]

########################################

# Plot of temperature

# Read and prepare data
ter_dat <- readRDS("data/summary_ind_ter_perc.rds")
head(ter_dat)
colnames(ter_dat) <- c("path", "var", "I-II", "III-IV", "Not-designated", "Total", "V-VI", "n")
ter_dat$sum <- rowSums(ter_dat[,c("I-II", "III-IV", "V-VI", "Not-designated")], na.rm=T)

#####

# How can sum of individual protection be smaller than total?
ter_dat[round(ter_dat$sum, 1) < round(ter_dat$Total,1),]

#####

ter_dat$`I-II` <- ifelse(ter_dat$sum > ter_dat$Total, 
                         ifelse(ter_dat$`I-II` > ter_dat$Total, ter_dat$Total, ter_dat$`I-II`), 
                         ter_dat$`I-II`)
ter_dat$`III-IV` <- ifelse(ter_dat$sum > ter_dat$Total, 
                           ifelse(ter_dat[,c("I-II")] == ter_dat$Total, 0,
                                  ifelse(rowSums(ter_dat[,c("I-II", "III-IV")], na.rm=T) >= ter_dat$Total,
                                         ter_dat$Total-ter_dat$`I-II`,
                                         ter_dat$`III-IV`)), 
                           ter_dat$`III-IV`)
ter_dat$`V-VI` <- ifelse(ter_dat$sum > ter_dat$Total, 
                         ifelse(rowSums(ter_dat[,c("I-II", "III-IV")], na.rm=T) == ter_dat$Total, 0, 
                                ifelse(rowSums(ter_dat[,c("I-II", "III-IV", "V-VI")], na.rm=T) >= ter_dat$Total,
                                       ter_dat$Total-rowSums(ter_dat[,c("I-II", "III-IV")], na.rm=T),
                                       ter_dat$`V-VI`)), 
                         ter_dat$`V-VI`)
ter_dat$`Not-designated` <- ifelse(ter_dat$sum > ter_dat$Total, 
                                   ifelse(rowSums(ter_dat[,c("I-II", "III-IV", "V-VI")], na.rm=T) == ter_dat$Total, 0,
                                          ifelse(rowSums(ter_dat[,c("I-II", "III-IV", "V-VI", "Not-designated")], na.rm=T) >= ter_dat$Total,
                                                 ter_dat$Total-rowSums(ter_dat[,c("I-II", "III-IV", "V-VI")], na.rm=T),
                                                 ter_dat$`Not-designated`)), 
                                   ter_dat$`Not-designated`)

# Summary
ter_dat <- ter_dat %>% dplyr::select(-c(Total, sum)) %>% 
  tidyr::gather(iucn_cat, perc, -c(path, var, n)) %>% 
  mutate(iucn_cat = factor(iucn_cat, levels=c("Not-designated", "V-VI", "III-IV", "I-II"),
                           labels=c("Non-designated", "V-VI", "III-IV", "I-II"))) %>% drop_na()

m_ee <- readRDS("data/summary_wc_perc_optim.rds")

# Round first and last value to include all ranges!
m_ee[1,2] <- floor(m_ee[1,2])
m_ee[102,2] <- floor(m_ee[102,2])
m_ee[203,2] <- floor(m_ee[203,2])
m <- m_ee %>% group_split(var)

ter_dat %>% filter(path=="bio01_perc") %>% 
  left_join(as.data.frame(m[[1]]), by=c("var"="z")) %>% tidyr::drop_na() %>% 
  ggplot(aes(x = var, y=perc)) + 
  geom_bar(stat="identity", position="stack", width=1) +
  labs(x="Annual mean temperature (°C)", y="% area protected") +
  scale_x_continuous(breaks=c(1,25,50,75,98), labels=round(rowMeans(m[[1]][c(1,25,50,75,98), 2:3]),1), 
                     limits=c(0,101), expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), 
                     limits=c(0,45)) +
  theme_bw() + theme(axis.title=element_text(size=16),
                     axis.text = element_text(size=12))
ggsave(paste0("figures/pa_example_bio1_bar.png"), height=4, width=5, dpi=1000)

########################################

# Plot terrestrial protection map

# Create maps
bio1 <- raster::raster("extdata/wc2.0_bio_30s_01_perc.nc")
bio1 <- raster::crop(bio1, pa_sub_moll, snap="out")

ggplot() + ggR(bio1, geom_raster=T, ggLayer=T) + 
  scale_fill_gradientn(colours = bluered, name = "Protection (%)",
                       na.value="transparent") + 
  geom_sf(data = pa_sub_moll, fill=NA, colour="black") + 
  geom_tile(data=df_pa_all, aes(x, y), fill="transparent", colour="black") + 
  scale_x_continuous() + scale_y_continuous() + 
  theme_minimal() + coord_sf(ndiscr=F) + 
  theme(axis.title=element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
ggsave("figures/pa_example_bio1_prot_map.png", height=3, width=4.5, dpi=600)

########################################