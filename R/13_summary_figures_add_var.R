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
colnames(ter_dat) <- c("path", "var", "I-II", "III-IV", "V-VI", "Not-designated", "Total", "n")

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

###
# NOTE: If some summaries still have NAs after left_join, check if re-classification is correct
###

####################

ter_dat %>% group_by(path, var) %>% group_keys

#' ### Hyper-geometric distribution

# Total area & Total area protected
(tot_sum <- ter_dat %>% group_by(path, var) %>% summarise(sum=sum(n)/4, prot_cells=sum((perc*n))) %>%
    ungroup() %>% group_by(path) %>% summarise(sum=sum(sum), prot_cells=sum(prot_cells)) %>%
    mutate(prop_prot=prot_cells/sum))

# Global area climate
(clim_area <- ter_dat %>% group_by(path,var) %>% summarise(area_clim=mean(n)) %>% 
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
sapply(m, nrow)

bio04_dat <- ter_dat %>% full_join(exp_val) %>% filter(path=="bio04_perc") %>% 
  full_join(as.data.frame(m[[2]]), by=c("var"="z"))

# goodness of fit test
test <- bio04_dat %>% group_by(path, var) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  summarise(chisq = round(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$statistic, 2),
            df = chisq.test(x=.$perc,p=.$exp, rescale.p=T)$parameter,
            pvalue = signif(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$p.value, 2))

p1 <- bio04_dat %>% ggplot(aes(x = var, y=perc, fill=iucn_cat)) + 
  geom_bar(stat="identity", position="stack", width=1) +
  geom_line(aes(x=var, y=exp)) + 
  geom_text_npc(npcx=0.04,npcy=0.95, label="(a)") + 
  geom_text_npc(npcx=0.5, npcy=0.95, size=2, label=paste0("X-squared = ", test$chisq,", \n df = ", 
                                                          test$df, ", \n p-value = ", test$pvalue)) + 
  labs(x="Temperature seasonality", y="% of area protected") +
  scale_x_continuous(breaks=c(1,25,50,75,98), labels=round(rowMeans(m[[2]][c(1,25,50,75,98),2:3]), 0), 
                     limits=c(0,101), expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), limits=c(0,35)) + 
  scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"), 
                    guide = guide_legend(reverse = TRUE)) + 
  theme_bw() + theme(legend.position="none", panel.grid.minor = element_blank())

bio07_dat <- ter_dat %>% full_join(exp_val) %>% filter(path=="bio07_perc") %>% 
  full_join(as.data.frame(m[[3]]), by=c("var"="z"))

# goodness of fit test
test <- bio07_dat %>% group_by(path, var) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  summarise(chisq = round(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$statistic, 2),
            df = chisq.test(x=.$perc,p=.$exp, rescale.p=T)$parameter,
            pvalue = signif(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$p.value, 2))

p2 <- bio07_dat %>% 
  ggplot(aes(x = var, y=perc, fill=iucn_cat)) + 
  geom_bar(stat="identity", position="stack", width=1) +
  geom_line(aes(x=var, y=exp)) + 
  geom_text_npc(npcx=0.04,npcy=0.95, label="(b)") + 
  geom_text_npc(npcx=0.5, npcy=0.95, size=2, label=paste0("X-squared = ", test$chisq,", \n df = ", 
                                                          test$df, ", \n p-value = ", test$pvalue)) + 
  labs(x="Temperature annual range", y="% of area protected") +
  scale_x_continuous(breaks=c(1,25,50,75,98), labels=round(rowMeans(m[[3]][c(1,25,50,75,98),2:3]), 1), 
                     limits=c(0,101), expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), limits=c(0,35)) + 
  scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"), 
                    guide = guide_legend(reverse = TRUE)) + 
  theme_bw() + theme(legend.position="none", axis.title.y=element_blank(), 
                     panel.grid.minor = element_blank())

bio15_dat <- ter_dat %>% full_join(exp_val) %>% filter(path=="bio15_perc") %>% 
  full_join(as.data.frame(m[[5]]), by=c("var"="z"))

# goodness of fit test
test <- bio15_dat %>% group_by(path, var) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  summarise(chisq = round(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$statistic, 2),
            df = chisq.test(x=.$perc,p=.$exp, rescale.p=T)$parameter,
            pvalue = signif(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$p.value, 2))

p3 <- bio15_dat %>% 
  mutate(var2=rowMeans(cbind(x,y))) %>%
  ggplot(aes(x = var, y=perc, fill=iucn_cat)) + 
  geom_bar(stat="identity", position="stack", width=1) +
  geom_line(aes(x=var, y=exp)) + 
  geom_text_npc(npcx=0.04,npcy=0.95, label="(c)") + 
  geom_text_npc(npcx=0.5, npcy=0.95, size=2, label=paste0("X-squared = ", test$chisq,", \n df = ", 
                                                          test$df, ", \n p-value = ", test$pvalue)) + 
  scale_x_continuous(breaks=c(1,25,50,75,98), labels=round(rowMeans(m[[5]][c(1,25,50,75,98),2:3]), 0), 
                     limits=c(0,101), expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), limits=c(0,35)) +
  scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"), 
                    guide = guide_legend(reverse = TRUE)) + 
  labs(x="Precipitation seasonality", y="", fill="IUCN") + 
  theme_bw() + theme(legend.position = "bottom", panel.grid.minor = element_blank())

leg <- ggpubr::as_ggplot(ggpubr::get_legend(p3))

p <- p1 + p2 + {p3 + theme(legend.position="none")} + plot_spacer() + leg + plot_spacer() + 
  plot_layout(ncol=3, heights=c(4,1))
ggsave("figures/ter_sum_add_var.png", p, dpi=1000, width=8, height=3)

#########################

# Plot of environmental variables for marine areas

# Read and prepare data

marspec_dat <- readRDS("data/summary_ind_mar_perc.rds")
head(marspec_dat)
colnames(marspec_dat) <- c("path", "var", "I-II", "III-IV",  "V-VI", "Not-designated", "Total","n")
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
(tot_sum <- marspec_dat %>% group_by(path, var) %>% summarise(sum=sum(n)/4, prot_cells=sum((perc*n))) %>%
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
seq(1, 97, length.out=5)

####################
biogeo16_dat <- marspec_dat %>% full_join(exp_val) %>% filter(path=="biogeo16_perc") %>% 
  left_join(as.data.frame(m[[6]]), by=c("var"="z"))

# goodness of fit test
test <- biogeo16_dat %>% group_by(path, var) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  summarise(chisq = round(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$statistic, 2),
            df = chisq.test(x=.$perc,p=.$exp, rescale.p=T)$parameter,
            pvalue = signif(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$p.value, 2))

p1 <- biogeo16_dat %>% ggplot() + geom_bar(aes(x = var, y=perc, fill=iucn_cat), width=1, stat="identity", position="stack") +
  geom_line(aes(x=var, y=exp), colour="black") + geom_text_npc(npcx=0.04,npcy=0.95, label="(a)") + 
  geom_text_npc(npcx=0.5, npcy=0.95, size=2, label=paste0("X-squared = ", test$chisq,", \n df = ", test$df, ", \n p-value = ", test$pvalue)) + 
  labs(x="Annual range in SST", y="% of area protected") + 
  scale_x_continuous(breaks=c(1,25,49,73,95), labels=round(rowMeans(m[[6]][c(1,25,49,73,95),2:3]), 1), 
                     limits=c(0,98), expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), limits=c(0,25)) + #, breaks=c(0, 1, 5, 15, 30)) + 
  scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"), 
                    guide = guide_legend(reverse = TRUE)) + 
  theme_bw() + theme(legend.position="none", panel.grid.minor = element_blank()) #+ coord_trans(y="sqrt")

biogeo17_dat <- marspec_dat %>% full_join(exp_val) %>% filter(path=="biogeo17_perc") %>% 
  left_join(as.data.frame(m[[7]]), by=c("var"="z"))

# goodness of fit test
test <- biogeo17_dat %>% group_by(path, var) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  summarise(chisq = round(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$statistic, 2),
            df = chisq.test(x=.$perc,p=.$exp, rescale.p=T)$parameter,
            pvalue = signif(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$p.value, 2))

p2 <- biogeo17_dat %>% 
  ggplot() + geom_bar(aes(x = var, y=perc, fill=iucn_cat), width=1, stat="identity", position="stack") +
  geom_line(aes(x=var, y=exp), colour="black") + geom_text_npc(npcx=0.04,npcy=0.95, label="(b)") + 
  geom_text_npc(npcx=0.5, npcy=0.95, size=2, label=paste0("X-squared = ", test$chisq,", \n df = ", 
                                                          test$df, ", \n p-value = ", test$pvalue)) + 
  labs(x="Annual variance in SST", y="% of area protected") + 
  scale_x_continuous(breaks=c(1,25,49,74,96), labels=round(rowMeans(m[[7]][c(1,25,49,74,96),2:3]), 0), 
                     limits=c(0,99), expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), limits=c(0,25)) +
  scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"), 
                    guide = guide_legend(reverse = TRUE)) + 
  theme_bw() + theme(legend.position="none", axis.title.y=element_blank(), 
                     panel.grid.minor = element_blank())

biogeo11_dat <- marspec_dat %>% full_join(exp_val) %>% filter(path=="biogeo11_perc")  %>% 
  left_join(as.data.frame(m[[3]]), by=c("var"="z"))

# goodness of fit test
test <- biogeo11_dat %>% group_by(path, var) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  summarise(chisq = round(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$statistic, 2),
            df = chisq.test(x=.$perc,p=.$exp, rescale.p=T)$parameter,
            pvalue = signif(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$p.value, 2))

p3 <- biogeo11_dat %>% 
  ggplot() + geom_bar(aes(x = var, y=perc, fill=iucn_cat), width=1, stat="identity", position="stack") +
  geom_line(aes(x=var, y=exp), colour="black") + geom_text_npc(npcx=0.04,npcy=0.95, label="(c)") + 
  geom_text_npc(npcx=0.5, npcy=0.95, size=2, label=paste0("X-squared = ", test$chisq,", \n df = ", 
                                                          test$df, ", \n p-value = ", test$pvalue)) + 
  labs(x="Annual range in SSS", y="% of area protected", fill="IUCN") + 
  scale_x_continuous(breaks=c(1,18,35,52,67), labels=round(rowMeans(m[[3]][c(1,18,35,52,67),2:3]), 1), 
                     limits=c(0,70), expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), limits=c(0,25)) + 
  scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"), 
                    guide = guide_legend(reverse = TRUE)) + theme_bw() + 
  theme(legend.position="bottom", panel.grid.minor = element_blank())

biogeo12_dat <- marspec_dat %>% full_join(exp_val) %>% filter(path=="biogeo12_perc") %>% 
  full_join(as.data.frame(m[[4]]), by=c("var"="z"))

# goodness of fit test
test <- biogeo12_dat %>% group_by(path, var) %>% summarise(perc=sum(perc), exp=mean(exp)) %>%
  summarise(chisq = round(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$statistic, 2),
          df = chisq.test(x=.$perc,p=.$exp, rescale.p=T)$parameter,
          pvalue = signif(chisq.test(x=.$perc,p=.$exp, rescale.p=T)$p.value, 2))

p4 <- biogeo12_dat %>% ggplot() + geom_bar(aes(x = var, y=perc, fill=iucn_cat), width=1, stat="identity", position="stack") +
  geom_line(aes(x=var, y=exp), colour="black") + geom_text_npc(npcx=0.04,npcy=0.95, label="(d)") + 
  geom_text_npc(npcx=0.5, npcy=0.95, size=2, label=paste0("X-squared = ", test$chisq,", \n df = ", 
                                                         test$df, ", \n p-value = ", test$pvalue)) + 
 labs(x="Annual variance in SSS", y="% of area protected") + 
 scale_x_continuous(breaks=c(1,21,41,61,81), labels=round(rowMeans(m[[4]][c(1,21,41,61,81),2:3]), 0), 
                   limits=c(0,82), expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), limits=c(0,25)) + 
  scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"), 
                    guide = guide_legend(reverse = TRUE)) + theme_bw() + 
 theme(legend.position="none", axis.title.y=element_blank(), 
       panel.grid.minor = element_blank())

p5 <- ggpubr::as_ggplot(ggpubr::get_legend(p3))

p <- (p1 | p2 ) / ({p3 + theme(legend.position="none")} | p4) / (plot_spacer() + p5 + plot_spacer()) + plot_layout(heights=c(4,4,1))
ggsave("figures/mar_sum_add_var.png", p, dpi=1000, width=8, height=6)
