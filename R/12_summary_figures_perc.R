#' ---
#' title: "Global Analysis of Protected Areas - Create environmental summary figures"
#' author: "RS-eco"
#' ---

rm(list=ls()); gc()

#Automatically install required packages, which are not yet installed
packages <- c("tidyverse", "patchwork", "ggpubr", "ggpmisc")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) instTotal.packages(new.packages); rm(new.packages)

# Load packages
l <- sapply(packages, require, character.only = TRUE, quietly=TRUE); rm(packages, l)

# Set working directory
workdir <- "C:/Users/admin/Documents/GitHub/globePA/"
setwd(workdir)

########################################

# Plot of temperature, precipitation/salinity and elevation/depth

# Read and prepare data
ter_dat <- readRDS("data/summary_ind_ter_perc.rds")
head(ter_dat)
colnames(ter_dat) <- c("path", "var", "I-II", "III-IV",  "V-VI", "Not-designated", "Total", "n")
ter_dat$sum <- rowSums(ter_dat[,c("I-II", "III-IV", "V-VI", "Not-designated")], na.rm=T)

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

ter_dat <- ter_dat %>% dplyr::select(-c(Total, sum)) %>% 
  tidyr::gather(iucn_cat, perc, -c(path, var, n)) %>% 
  mutate(iucn_cat = factor(iucn_cat, levels=c("Not-designated", "V-VI", "III-IV", "I-II"),
                           labels=c("Non-designated", "V-VI", "III-IV", "I-II"))) %>% drop_na()

# Check number of cells!
unique(ter_dat$n)

ter_dat %>% filter(path=="bio12_perc", iucn_cat == "I-II") %>% 
  ungroup() %>% dplyr::select(perc) %>% summary()

# Summary
ter_dat %>% group_by(path, var) %>% summarise(sum=sum(perc)) %>% 
  summarise(max(sum))

# Number of cells
ter_dat %>% group_by(path) %>% summarise(total_cells=sum(n)/4)

# Plot frequeny plots of envdata
ter_dat %>% ggplot(aes(x = var, y=perc, fill=iucn_cat)) + 
  geom_area(stat="identity", position="stack") + facet_wrap(.~ path, scales="free") + 
  labs(x="", y="% of area protected ") + theme_bw() + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=expansion(mult=c(0,.1)))

ter_dat %>% group_by(path)

####################

#' ### Hyper-geometric distribution

# Total area & Total area protected
(tot_sum <- ter_dat %>% group_by(path, var) %>% summarise(sum=sum(n)/n(), prot_cells=sum((perc*n))) %>%
   ungroup() %>% group_by(path) %>% summarise(sum=sum(sum), prot_cells=sum(prot_cells)) %>%
   mutate(prop_prot=prot_cells/sum))

# Global area climate
(clim_area <- ter_dat %>% group_by(path,var) %>% summarise(area_clim=sum(n)/n()) %>% 
    left_join(tot_sum))

# Expected proportion climate  
(exp_val <- clim_area %>% mutate(prop_clim = area_clim/sum*100) %>% 
    mutate(exp = (prop_clim*prop_prot),
           exp_aichi = 15*prop_clim,
           var_exp = (prop_clim*prop_prot*(1-prop_clim)*(1-prop_prot)/(sum-1))))

# Need to multiply by 100 to get perc value rather than proportion

# Need to multiply by bin size to get values up to 100 %
exp_val %>% group_by(path) %>% summarise(sum(prop_clim))
exp_val %>% group_by(path) %>% summarise(sum(exp))
exp_val %>% group_by(path) %>% summarise(sum(exp_aichi))

####################

m_ee <- readRDS("data/summary_wc_perc_optim.rds")

# Round first and last value to include all ranges!
m_ee[1,2] <- floor(m_ee[1,2])
m_ee[102,2] <- floor(m_ee[102,2])
m_ee[203,2] <- floor(m_ee[203,2])
m <- m_ee %>% group_split(var)

bio01_dat <- ter_dat %>% filter(path=="bio01_perc") %>% left_join(exp_val) %>% 
  left_join(as.data.frame(m[[1]]), by=c("var"="z"))

###
# NOTE: If some summaries still have NAs after left_join, check if re-classification is correct
###

# goodness of fit test
test <- bio01_dat %>% group_by(path, var) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  summarise(chisq = round(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$statistic, 2),
            df = chisq.test(x=.$perc,p=.$exp, rescale.p=T)$parameter,
            pvalue = signif(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$p.value, 2))

p1 <- bio01_dat %>% ggplot() + 
  geom_bar(aes(x = var, y=perc, fill=iucn_cat), width=1, stat="identity", position="stack") +
  geom_line(aes(x=var, y=exp), colour="black") + geom_text_npc(npcx=0.04,npcy=0.95, label="(a)") + 
  geom_text_npc(npcx=0.5, npcy=0.95, size=2, label=paste0("X-squared = ", test$chisq,", \n df = ", test$df, ", \n p-value = ", test$pvalue)) + 
  labs(x="Annual mean temp. (°C)", y="% of area protected") + # Annual mean temperature is too long
  scale_x_continuous(breaks=c(1,25,50,75,98), labels=round(rowMeans(m[[1]][c(1,25,50,75,98), 2:3]),1), 
                     limits=c(0,101), expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), limits=c(0,35)) + # , breaks=c(0, 1, 5, 15, 30)) +
  scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"), 
                    guide = guide_legend(reverse = TRUE)) + 
  theme_bw() + theme(legend.position="none", panel.grid.minor = element_blank()) #+ coord_trans(y="sqrt")
p1

bio01_dat %>% group_by(x,y) %>% summarise(perc=sum(perc), exp=mean(exp)) %>% 
  filter(perc < exp) %>% view()

bio12_dat <- ter_dat %>% filter(path=="bio12_perc") %>% left_join(exp_val) %>% 
  left_join(as.data.frame(m[[4]]), by=c("var"="z"))

# goodness of fit test
test <- bio12_dat %>% group_by(path, var) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  summarise(chisq = round(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$statistic, 2),
            df = chisq.test(x=.$perc,p=.$exp, rescale.p=T)$parameter,
            pvalue = signif(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$p.value, 2))

p2 <- bio12_dat %>% ggplot() + 
  geom_bar(aes(x = var, y=perc, fill=iucn_cat), width=1, stat="identity", position="stack") +
  geom_line(aes(x=var, y=exp), colour="black") + 
  geom_text_npc(npcx=0.5, npcy=0.95, size=2, label=paste0("X-squared = ", test$chisq,", \n df = ", test$df, ", \n p-value = ", test$pvalue)) + 
  geom_text_npc(npcx=0.04,npcy=0.95, label="(b)") + labs(x="Annual precipitation (mm)") +
  scale_x_continuous(breaks=c(1,25,50,75,98), labels=round(rowMeans(m[[4]][c(1,25,50,75,98),2:3]), 0), 
                     limits=c(0,101), expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), limits=c(0,35)) + #, breaks=c(0, 1, 5, 15, 30)) +
  scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"), 
                    guide = guide_legend(reverse = TRUE)) + 
  theme_bw() + theme(legend.position="none", axis.title.y=element_blank(), 
                     panel.grid.minor = element_blank()) #+ coord_trans(y="sqrt")
p2

bio12_dat %>% group_by(x,y) %>% summarise(perc=sum(perc), exp=mean(exp)) %>% 
  filter(perc < exp) %>% view()

m <- readRDS("data/summary_earthenv_perc_optim.rds")

alt_dat <- ter_dat %>% filter(path=="elevation_perc") %>% left_join(exp_val) %>% 
  left_join(as.data.frame(m), by=c("var"="z"))

# goodness of fit test
test <- alt_dat %>% group_by(path, var) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  summarise(chisq = round(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$statistic, 2),
            df = chisq.test(x=.$perc,p=.$exp, rescale.p=T)$parameter,
            pvalue = signif(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$p.value, 2))

p3 <- alt_dat %>% ggplot() + 
  geom_bar(aes(x = var, y=perc, fill=iucn_cat), width=1, stat="identity", position="stack") +
  geom_line(aes(x=var, y=exp), colour="black") + geom_text_npc(npcx=0.04,npcy=0.95, label="(c)") + 
  geom_text_npc(npcx=0.5, npcy=0.95, size=2, label=paste0("X-squared = ", test$chisq,", \n df = ", test$df, ", \n p-value = ", test$pvalue)) + 
  labs(x="Elevation (m)", y="", fill="IUCN") +
  scale_x_continuous(breaks=c(1,25,50,75,98), labels=round(rowMeans(m[c(1,25,50,75,98),2:3]), 0), 
                     limits=c(0,101), expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), limits=c(0,35)) + #, breaks=c(0, 1, 5, 15, 30)) +
  scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"), 
                    guide = guide_legend(reverse = TRUE)) + 
  theme_bw() + theme(legend.position="bottom", axis.title.y=element_blank(), 
                     panel.grid.minor = element_blank()) #+ coord_trans(y="sqrt")
p3

alt_dat %>% group_by(x,y) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  filter(perc < exp) %>% view()

#########################

# Plot of environmental variables for marine areas

# Read and prepare data

marspec_dat <- readRDS("data/summary_ind_mar_perc.rds")
head(marspec_dat)
unique(marspec_dat$path)
marspec_dat %>% filter(path == "biogeo17_perc") %>% group_by(var) %>% group_keys() %>% unlist()

colnames(marspec_dat) <- c("path", "var", "I-II", "III-IV", "V-VI", "Not-designated", "Total", "n")
marspec_dat$sum <- rowSums(marspec_dat[,c("I-II", "III-IV", "V-VI", "Not-designated")], na.rm=T)

marspec_dat$`I-II` <- ifelse(marspec_dat$sum > marspec_dat$Total, 
                             ifelse(marspec_dat$`I-II` > marspec_dat$Total, marspec_dat$Total, marspec_dat$`I-II`), 
                             marspec_dat$`I-II`)
marspec_dat$`III-IV` <- ifelse(marspec_dat$sum > marspec_dat$Total, 
                               ifelse(marspec_dat[,c("I-II")] == marspec_dat$Total, 0,
                                      ifelse(rowSums(marspec_dat[,c("I-II", "III-IV")], na.rm=T) >= marspec_dat$Total,
                                             marspec_dat$Total-marspec_dat$`I-II`,
                                             marspec_dat$`III-IV`)), 
                               marspec_dat$`III-IV`)
marspec_dat$`V-VI` <- ifelse(marspec_dat$sum > marspec_dat$Total, 
                             ifelse(rowSums(marspec_dat[,c("I-II", "III-IV")], na.rm=T) == marspec_dat$Total, 0, 
                                    ifelse(rowSums(marspec_dat[,c("I-II", "III-IV", "V-VI")], na.rm=T) >= marspec_dat$Total,
                                           marspec_dat$Total-rowSums(marspec_dat[,c("I-II", "III-IV")], na.rm=T),
                                           marspec_dat$`V-VI`)), 
                             marspec_dat$`V-VI`)
marspec_dat$`Not-designated` <- ifelse(marspec_dat$sum > marspec_dat$Total, 
                                       ifelse(rowSums(marspec_dat[,c("I-II", "III-IV", "V-VI")], na.rm=T) == marspec_dat$Total, 0,
                                              ifelse(rowSums(marspec_dat[,c("I-II", "III-IV", "V-VI", "Not-designated")], na.rm=T) >= marspec_dat$Total,
                                                     marspec_dat$Total-rowSums(marspec_dat[,c("I-II", "III-IV", "V-VI")], na.rm=T),
                                                     marspec_dat$`Not-designated`)), 
                                       marspec_dat$`Not-designated`)

###

# Why are some areas 0, but have a value for protection???
# only the case for biogeo08

# => Re-check this is still the case!

###

marspec_dat <- marspec_dat %>% dplyr::select(-c(Total, sum)) %>%
  tidyr::gather(iucn_cat, perc, -c(path, var, n)) %>%
  mutate(iucn_cat = factor(iucn_cat, levels=c("Not-designated", "V-VI", "III-IV", "I-II"),
                           labels=c("Non-designated", "V-VI", "III-IV", "I-II"))) %>% drop_na()

# Summary should be maximum of 100!!!
marspec_dat %>% group_by(path, var) %>% summarise(sum=sum(perc)) %>% 
  summarise(max(sum))

# Number of cells
# Need to divide area by 4, as we have 4 categories!!!
marspec_dat %>% group_by(path) %>% summarise(total_cells=sum(n)/4)

# Plot frequeny plots of envdata
marspec_dat %>% 
  ggplot(aes(x = var, y=perc, fill=iucn_cat)) + 
  geom_area(stat="identity", position="stack") + facet_wrap(.~ path, scales="free") + 
  labs(x="", y="% of area protected") + theme_bw() + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=expansion(mult=c(0,.1)))

####################

#' ### Hyper-geometric distribution

# Total area & Total area protected
(tot_sum <- marspec_dat %>% group_by(path, var) %>% summarise(sum=sum(n)/n(), prot_cells=sum((perc*n))) %>%
   ungroup() %>% group_by(path) %>% summarise(sum=sum(sum), prot_cells=sum(prot_cells)) %>%
   mutate(prop_prot=prot_cells/sum))

# Global area climate
(clim_area <- marspec_dat %>% group_by(path,var) %>% summarise(area_clim=sum(n)/n()) %>% 
    left_join(tot_sum))

# Expected proportion climate  
(exp_val <- clim_area %>% mutate(prop_clim = area_clim/sum*100) %>% 
    mutate(exp = (prop_clim*prop_prot),
           exp_aichi = 15*prop_clim,
           var_exp = (prop_clim*prop_prot*(1-prop_clim)*(1-prop_prot)/(sum-1))))

# Need to multiply by 100 to get perc value rather than proportion

# Need to multiply by bin size to get values up to 100 %
exp_val %>% group_by(path) %>% summarise(sum(prop_clim))
exp_val %>% group_by(path) %>% summarise(sum(exp))
exp_val %>% group_by(path) %>% summarise(sum(exp_aichi))

####

m_ee <- readRDS("data/summary_marspec_perc_optim.rds")

# Round first and last value to include all ranges!
m_ee[1,2] <- m_ee[1,2]-1
m_ee[185,2] <- m_ee[185,2]-1
m_ee[335,2] <- m_ee[335,2]-1
m_ee[431,2] <- floor(m_ee[431,2])
m_ee[528,2] <- floor(m_ee[528,2])
m <- m_ee %>% group_split(var)

sapply(m, nrow)

####################

biogeo13_dat <- marspec_dat %>% filter(path=="biogeo13_perc") %>% left_join(exp_val) %>% 
  left_join(as.data.frame(m[[5]]), by=c("var"="z"))

# goodness of fit test
test <- biogeo13_dat %>% group_by(path, var) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  summarise(chisq = round(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$statistic, 2),
            df = chisq.test(x=.$perc,p=.$exp, rescale.p=T)$parameter,
            pvalue = signif(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$p.value, 2))

p4 <- biogeo13_dat %>% ggplot() + geom_bar(aes(x = var, y=perc, fill=iucn_cat), width=1, stat="identity", position="stack") +
  geom_line(aes(x=var, y=exp), colour="black") + geom_text_npc(npcx=0.04,npcy=0.95, label="(d)") + 
  geom_text_npc(npcx=0.5, npcy=0.95, size=2, label=paste0("X-squared = ", test$chisq,", \n df = ", test$df, ", \n p-value = ", test$pvalue)) + 
  labs(x="Mean annual SST (°C)", y="% of area protected") + 
  scale_x_continuous(breaks=c(1,24,48,72,94), labels=round(rowMeans(m[[5]][c(1,24,48,72,94),2:3])/100, 1), 
                     limits=c(0,97), expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), limits=c(0,25)) + #, breaks=c(0, 1, 5, 15, 30)) +
  scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"), 
                    guide = guide_legend(reverse = TRUE)) + 
  theme_bw() + theme(legend.position="none", panel.grid.minor = element_blank()) #+ coord_trans(y="sqrt")
p4

biogeo13_dat %>% group_by(x,y) %>% summarise(perc=sum(perc), exp=mean(exp))

biogeo08_dat <- marspec_dat %>% filter(path=="biogeo08_perc") %>% left_join(exp_val) %>% 
  left_join(as.data.frame(m[[2]]), by=c("var"="z"))

# goodness of fit test
test <- biogeo08_dat %>% group_by(path, var) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  summarise(chisq = round(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$statistic, 2),
            df = chisq.test(x=.$perc,p=.$exp, rescale.p=T)$parameter,
            pvalue = signif(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$p.value, 2))

p5 <- biogeo08_dat %>% ggplot() + geom_bar(aes(x = var, y=perc, fill=iucn_cat), width=1, stat="identity", position="stack") +
  geom_text_npc(npcx=0.5, npcy=0.95, size=2, label=paste0("X-squared = ", test$chisq,", \n df = ", test$df, ", \n p-value = ", test$pvalue)) + 
  geom_line(aes(x=var, y=exp), colour="black") + geom_text_npc(npcx=0.04,npcy=0.95, label="(e)") + 
  labs(x="Mean annual SSS (psu)") + 
  scale_x_continuous(breaks=c(1,21,42,63,82), labels=round(rowMeans(m[[2]][c(1,21,42,63,82),2:3])/100, 1), 
                     limits=c(0,85), expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), limits=c(0,25)) + #, breaks=c(0, 1, 5, 15, 30)) +
  scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"), 
                    guide = guide_legend(reverse = TRUE)) + 
  theme_bw() + theme(legend.position="none", axis.title.y=element_blank(), 
                     panel.grid.minor = element_blank()) #+ coord_trans(y="sqrt")

p5
biogeo08_dat %>% group_by(x,y) %>% summarise(perc=sum(perc), exp=mean(exp)) %>% 
  filter(perc < exp) %>% view()

bathy_dat <- marspec_dat %>% filter(path=="bathy_perc") %>% left_join(exp_val) %>% 
  left_join(as.data.frame(m[[1]]), by=c("var"="z"))

# goodness of fit test
test <- bathy_dat %>% group_by(path, var) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  summarise(chisq = round(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$statistic, 2),
            df = chisq.test(x=.$perc,p=.$exp, rescale.p=T)$parameter,
            pvalue = signif(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$p.value, 2))

p6 <- bathy_dat %>% 
  ggplot() + geom_bar(aes(x = var, y=perc, fill=iucn_cat), width=1, stat="identity", position="stack") +
  geom_line(aes(x=var, y=exp), colour="black") + geom_text_npc(npcx=0.04,npcy=0.95, label="(f)") + 
  geom_text_npc(npcx=0.5, npcy=0.95, size=2, label=paste0("X-squared = ", test$chisq,", \n df = ", test$df, ", \n p-value = ", test$pvalue)) + 
  scale_x_continuous(breaks=c(1,25,50,75,98), labels=round(rowMeans(m[[1]][c(1,25,50,75,98),2:3]), 0), 
                     limits=c(0,101), expand=c(0,0)) + labs(x="Bathymetry (m)", y="") + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), limits=c(0,25)) + #, breaks=c(0, 1, 5, 15, 30)) +
  scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"), 
                    guide = guide_legend(reverse = TRUE)) + 
  theme_bw() + theme(legend.position="none", axis.title.y=element_blank(),
                     panel.grid.minor = element_blank()) #+ coord_trans(y="sqrt")

p6
bathy_dat %>% group_by(x,y) %>% summarise(perc=sum(perc), exp=mean(exp)) %>% 
  filter(perc < exp) %>% view()

p7 <- ggpubr::as_ggplot(ggpubr::get_legend(p3))

p <- p1 + p2 + {p3 + theme(legend.position="none")} + p4 + p5 + p6 + plot_spacer() + p7 + plot_spacer() + 
  plot_layout(ncol=3, heights=c(4,4,1))
ggsave("figures/Figure2_perc.png", p, dpi=1000, width=8, height=5)
ggsave("figures/Figure2.pdf", p, dpi=1000, width=8, height=5)
