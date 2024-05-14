#=========================##
## PERMANOVAS SIMPER NMDS ##
#=========================#

library(dplyr)
library(tidyr)
library(data.table)
library(tidyverse)
library(codyn)
library(vegan)
library(lme4)


#Other packages to load
library(emmeans)
library(lmerTest)
library(MuMIn)
library(AICcmodavg)

#DOWNLOADING PAIRWISE ADONIS
#install.packages("devtools")
library("devtools")

#remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)


#Use spcomp_long.final dataframe

spcomp_long.final<-read.csv("data/data_filtered_names.csv")
spcomp_long.final$Trt<- factor(spcomp_long.final$Trt, levels=c('con', 'chr', 'int'))
spcomp_long.final$Species <- NULL

#Make the following items into a factor:
spcomp_long.final<-spcomp_long.final %>% 
  dplyr::mutate_at(c('Site','Year','Trt','Plot','Block','Spcode'),as.factor)

#all years one model; not separated by period 
spcomp_long.final.drt<- spcomp_long.final %>%
  dplyr::filter(Year %in% c(2013,2014,2015,2016,2017)) 

spcomp_long.final.rec<- spcomp_long.final %>%
  dplyr::filter(Year %in% c(2017,2018,2019,2020,2021)) 


#separate by site and period for matrix 

spcomp_long.drt_CHY <- spcomp_long.final.drt %>% dplyr::filter(Site=="CHY")
spcomp_long.drt_KNZ <- spcomp_long.final.drt %>% dplyr::filter(Site=="KNZ")
spcomp_long.drt_HYS <- spcomp_long.final.drt %>% dplyr::filter(Site=="HYS")
spcomp_long.drt_SGS <- spcomp_long.final.drt %>% dplyr::filter(Site=="SGS")

spcomp_long.rec_CHY <- spcomp_long.final.rec %>% dplyr::filter(Site=="CHY")
spcomp_long.rec_KNZ <- spcomp_long.final.rec %>% dplyr::filter(Site=="KNZ")
spcomp_long.rec_HYS <- spcomp_long.final.rec %>% dplyr::filter(Site=="HYS")
spcomp_long.rec_SGS <- spcomp_long.final.rec %>% dplyr::filter(Site=="SGS")

#########################
#### VEGAN ANALYSES #####
#########################

CHY_wide.drt<-spread(spcomp_long.drt_CHY,Spcode,avg.cover)
KNZ_wide.drt<-spread(spcomp_long.drt_KNZ,Spcode,avg.cover)
HYS_wide.drt<-spread(spcomp_long.drt_HYS,Spcode,avg.cover)
SGS_wide.drt<-spread(spcomp_long.drt_SGS,Spcode,avg.cover)

CHY_wide.rec<-spread(spcomp_long.rec_CHY,Spcode,avg.cover)
KNZ_wide.rec<-spread(spcomp_long.rec_KNZ,Spcode,avg.cover)
HYS_wide.rec<-spread(spcomp_long.rec_HYS,Spcode,avg.cover)
SGS_wide.rec<-spread(spcomp_long.rec_SGS,Spcode,avg.cover)

#Change NAs to zero
CHY_wide.drt[is.na(CHY_wide.drt)] <- 0
KNZ_wide.drt[is.na(KNZ_wide.drt)] <- 0
HYS_wide.drt[is.na(HYS_wide.drt)] <- 0
SGS_wide.drt[is.na(SGS_wide.drt)] <- 0

CHY_wide.rec[is.na(CHY_wide.rec)] <- 0
KNZ_wide.rec[is.na(KNZ_wide.rec)] <- 0
HYS_wide.rec[is.na(HYS_wide.rec)] <- 0
SGS_wide.rec[is.na(SGS_wide.rec)] <- 0


#community matrix 

CHY.comm.drt<-CHY_wide.drt[,6:54] 
KNZ.comm.drt<-KNZ_wide.drt[,6:50] 
HYS.comm.drt<-HYS_wide.drt[,6:76]
SGS.comm.drt<-SGS_wide.drt[,6:44]

CHY.comm.rec<-CHY_wide.rec[,6:72] 
KNZ.comm.rec<-KNZ_wide.rec[,6:63] 
HYS.comm.rec<-HYS_wide.rec[,6:82]
SGS.comm.rec<-SGS_wide.rec[,6:41]

#PERMANOVA
adonis2(CHY.comm.drt~CHY_wide.drt$Year*CHY_wide.drt$Trt, method="bray", permutations=1000, Blocks(CHY_wide.drt$Block))
adonis2(KNZ.comm.drt~KNZ_wide.drt$Year*KNZ_wide.drt$Trt, method="bray", permutations=1000, Blocks(KNZ_wide.drt$Block))
adonis2(HYS.comm.drt~HYS_wide.drt$Year*HYS_wide.drt$Trt, method="bray", permutations=1000, Blocks(HYS_wide.drt$Block))
adonis2(SGS.comm.drt~SGS_wide.drt$Year*SGS_wide.drt$Trt, method="bray", permutations=1000, Blocks(SGS_wide.drt$Block))

adonis2(CHY.comm.rec~CHY_wide.rec$Year*CHY_wide.rec$Trt, method="bray", permutations=1000, Blocks(CHY_wide.rec$Block))
adonis2(KNZ.comm.rec~KNZ_wide.rec$Year*KNZ_wide.rec$Trt, method="bray", permutations=1000, Blocks(KNZ_wide.rec$Block))
adonis2(HYS.comm.rec~HYS_wide.rec$Year*HYS_wide.rec$Trt, method="bray", permutations=1000, Blocks(HYS_wide.rec$Block))
adonis2(SGS.comm.rec~SGS_wide.rec$Year*SGS_wide.rec$Trt, method="bray", permutations=1000, Blocks(SGS_wide.rec$Block))


factor.drt.CHY<-CHY_wide.drt[,1:5]
factor.drt.KNZ<-KNZ_wide.drt[,1:5]
factor.drt.HYS<-HYS_wide.drt[,1:5]
factor.drt.SGS<-SGS_wide.drt[,1:5]

factor.rec.CHY<-CHY_wide.rec[,1:5]
factor.rec.KNZ<-KNZ_wide.rec[,1:5]
factor.rec.HYS<-HYS_wide.rec[,1:5]
factor.rec.SGS<-SGS_wide.rec[,1:5]

#Pairwise comparisons 

pairwise.adonis2(CHY.comm.drt~Trt, data=factor.drt.CHY, sim.method="bray",p.adjust.m = "bonferroni")
pairwise.adonis2(CHY.comm.rec~Trt, data=factor.rec.CHY, sim.method="bray",p.adjust.m = "bonferroni")

pairwise.adonis2(KNZ.comm.drt~Trt, data=factor.drt.KNZ, sim.method="bray",p.adjust.m = "bonferroni")
pairwise.adonis2(KNZ.comm.rec~Trt, data=factor.rec.KNZ, sim.method="bray",p.adjust.m = "bonferroni")

pairwise.adonis2(HYS.comm.drt~Trt, data=factor.drt.HYS, sim.method="bray",p.adjust.m = "bonferroni")
pairwise.adonis2(HYS.comm.rec~Trt, data=factor.rec.HYS, sim.method="bray",p.adjust.m = "bonferroni")

pairwise.adonis2(SGS.comm.drt~Trt, data=factor.drt.SGS, sim.method="bray",p.adjust.m = "bonferroni")
pairwise.adonis2(SGS.comm.rec~Trt, data=factor.rec.SGS, sim.method="bray",p.adjust.m = "bonferroni")

#SGS with interation 
perm <- how(nperm = 999)
setBlocks(perm) <- with(factor.drt.SGS, factor.drt.SGS$plot)#this is a workaround to essentially include plot as a random effect in the permanovas
factor.drt.SGS <- unite(factor.drt.SGS, comb, c("Year", "Trt"), remove = FALSE)
mod.SGS <- pairwise.adonis2(SGS.comm.drt~comb, data=factor.SGS, sim.method="bray",p.adjust.m = "bonferroni", permutations = perm)

#Drought 
#2013
mod.SGS$`2013_con_vs_2013_chr`[1,5]
mod.SGS$`2013_con_vs_2013_int`[1,5]
mod.SGS$`2013_chr_vs_2013_int`[1,5]
#2014
mod.SGS$`2014_con_vs_2014_chr`[1,5]
mod.SGS$`2014_con_vs_2014_int`[1,5]
mod.SGS$`2014_chr_vs_2014_int`[1,5]
#2015
mod.SGS$`2015_con_vs_2015_chr`[1,5]
mod.SGS$`2015_con_vs_2015_int`[1,5]
mod.SGS$`2015_chr_vs_2015_int`[1,5]
#2016
mod.SGS$`2016_con_vs_2016_chr`[1,5]
mod.SGS$`2016_con_vs_2016_int`[1,5]
mod.SGS$`2016_chr_vs_2016_int`[1,5]
#2017
mod.SGS$`2017_con_vs_2017_chr`[1,5]
mod.SGS$`2017_con_vs_2017_int`[1,5]
mod.SGS$`2017_chr_vs_2017_int`[1,5]

#Recovery 
perm <- how(nperm = 999)
setBlocks(perm) <- with(factor.rec.SGS, factor.rec.SGS$plot)#this is a workaround to essentially include plot as a random effect in the permanovas
factor.rec.SGS <- unite(factor.rec.SGS, comb, c("Year", "Trt"), remove = FALSE)
mod.SGS.rec <- pairwise.adonis2(SGS.comm.rec~comb, data=factor.rec.SGS, sim.method="bray",p.adjust.m = "bonferroni", permutations = perm)

#2018
mod.SGS.rec$`2018_con_vs_2018_chr`[1,5]
mod.SGS.rec$`2018_con_vs_2018_int`[1,5]
mod.SGS.rec$`2018_chr_vs_2018_int`[1,5]
#2019
mod.SGS$`2019_con_vs_2019_chr`[1,5]
mod.SGS$`2019_con_vs_2019_int`[1,5]
mod.SGS$`2019_chr_vs_2019_int`[1,5]
#2020
mod.SGS$`2020_con_vs_2020_chr`[1,5]
mod.SGS$`2020_con_vs_2020_int`[1,5]
mod.SGS$`2020_chr_vs_2020_int`[1,5]
#2021
mod.SGS$`2021_con_vs_2021_chr`[1,5]
mod.SGS$`2021_con_vs_2021_int`[1,5]
mod.SGS$`2021_chr_vs_2021_int`[1,5]

##############
####SIMPER####
##############

#all years model
#SGS
#sig trt yrs: 2015-2021 
sgswide_sigyrs <- SGS_wide %>%
  dplyr::filter(Year %in% c(2015,2016,2017,2018,2019,2020,2021)) 

SGS.comm_sigyrs<-sgswide_sigyrs[,6:49]

(sim<-with(sgswide_sigyrs,simper(SGS.comm_sigyrs,Trt))) 
summary(sim,ordered=TRUE,digits=3)

#CHY 
#sig trt yrs: 2014-2021

chywide_sigyrs <- CHY_wide %>%
  dplyr::filter(Year %in% c(2015,2016,2017,2018,2019,2020,2021)) 

CHY.comm_sigyrs<-chywide_sigyrs[,6:76]

(sim<-with(chywide_sigyrs,simper(CHY.comm_sigyrs,Trt))) 
summary(sim,ordered=TRUE,digits=3)

#HYS
#sig trt yrs: 2015-2019
hyswide_sigyrs <- HYS_wide %>%
  dplyr::filter(Year %in% c(2015,2016,2017,2018,2019)) 

HYS.comm_sigyrs<-hyswide_sigyrs[,6:91]

(sim<-with(hyswide_sigyrs,simper(HYS.comm_sigyrs,Trt))) 
summary(sim,ordered=TRUE,digits=3)

#KNZ
#sig trt years: 2016-2018, 2020-2021
knzwide_sigyrs <- KNZ_wide %>%
  dplyr::filter(Year %in% c(2016,2017,2018,2020,2021)) 

KNZ.comm_sigyrs<-knzwide_sigyrs[,6:65]
(sim<-with(knzwide_sigyrs,simper(KNZ.comm_sigyrs,Trt))) 
summary(sim,ordered=TRUE,digits=3)

### drought model SIMPER ###
#all years model
#SGS
#sig trt yrs: 2015-2021 
sgswide_sigyrs <- SGS_wide %>%
  dplyr::filter(Year %in% c(2015,2016,2017,2018,2019,2020,2021)) 

SGS.comm_sigyrs<-sgswide_sigyrs[,6:49]

(sim<-with(sgswide_sigyrs,simper(SGS.comm_sigyrs,Trt))) 
summary(sim,ordered=TRUE,digits=3)

#CHY 
#sig trt yrs: 2014-2021

chywide_sigyrs <- CHY_wide %>%
  dplyr::filter(Year %in% c(2015,2016,2017,2018,2019,2020,2021)) 

CHY.comm_sigyrs<-chywide_sigyrs[,6:76]

(sim<-with(chywide_sigyrs,simper(CHY.comm_sigyrs,Trt))) 
summary(sim,ordered=TRUE,digits=3)

#HYS
#sig trt yrs: 2015-2019
hyswide_sigyrs <- HYS_wide %>%
  dplyr::filter(Year %in% c(2015,2016,2017,2018,2019)) 

HYS.comm_sigyrs<-hyswide_sigyrs[,6:91]

(sim<-with(hyswide_sigyrs,simper(HYS.comm_sigyrs,Trt))) 
summary(sim,ordered=TRUE,digits=3)

#KNZ
#sig trt years: 2016-2018, 2020-2021
knzwide_sigyrs <- KNZ_wide %>%
  dplyr::filter(Year %in% c(2016,2017,2018,2020,2021)) 

KNZ.comm_sigyrs<-knzwide_sigyrs[,6:65]
(sim<-with(knzwide_sigyrs,simper(KNZ.comm_sigyrs,Trt))) 
summary(sim,ordered=TRUE,digits=3)

########### NMDS PANEL FIGS ##############

#NMDS 
CHY.mds<-metaMDS(CHY.comm, distance = "bray", k = 3, maxit = 999, trymax = 500)
KNZ.mds<-metaMDS(KNZ.comm, distance = "bray", k = 3, maxit = 999, trymax = 500)
HYS.mds<-metaMDS(HYS.comm, distance = "bray", k = 3, maxit = 999, trymax = 500)
SGS.mds<-metaMDS(SGS.comm, distance = "bray", k = 3, maxit = 999, trymax = 500)

plot(CHY.mds$points)
plot(KNZ.mds$points)
plot(HYS.mds$points)
plot(SGS.mds$points)


CHY.mds_xy<- data.frame(CHY.mds$points)
CHY.mds_xy$Site <- CHY_wide$Site
CHY.mds_xy$Year <- CHY_wide$Year
CHY.mds_xy$Plot <- CHY_wide$Plot
CHY.mds_xy$Trt <- CHY_wide$Trt
#CHY.mds_xy$Block <- CHYcomm_wide$Block

KNZ.mds_xy<- data.frame(KNZ.mds$points)
KNZ.mds_xy$Site <- KNZ_wide$Site
KNZ.mds_xy$Year <- KNZ_wide$Year
KNZ.mds_xy$Plot <- KNZ_wide$Plot
KNZ.mds_xy$Trt <- KNZ_wide$Trt
#KNZ.mds_xy$Block <- KNZcomm_wide$Block

HYS.mds_xy<- data.frame(HYS.mds$points)
HYS.mds_xy$Site <- HYS_wide$Site
HYS.mds_xy$Year <- HYS_wide$Year
HYS.mds_xy$Plot <- HYS_wide$Plot
HYS.mds_xy$Trt <- HYS_wide$Trt
#HYS.mds_xy$Block <- HYScomm_wide$Block

SGS.mds_xy<- data.frame(SGS.mds$points)
SGS.mds_xy$Site <- SGS_wide$Site
SGS.mds_xy$Year <- SGS_wide$Year
SGS.mds_xy$Plot <- SGS_wide$Plot
SGS.mds_xy$Trt <- SGS_wide$Trt
#SGS.mds_xy$Block <- SGScomm_wide$Block

allNMDS<- rbind(SGS.mds_xy, CHY.mds_xy,HYS.mds_xy,KNZ.mds_xy)
allNMDS$Trt <- factor(allNMDS$Trt, levels=c('con', 'chr', 'int'))
allNMDS$Site <- factor(allNMDS$Site, levels=c('SGS', 'CHY', 'HYS','KNZ'))

#all sites panel 

NMDS_all<-ggplot(subset(allNMDS,Year %in% c("2013","2017","2021")), aes(x=MDS1, y=MDS2, color=Trt,shape=Trt)) + 
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Chronic", "Control","Intense"))+
  scale_shape_manual(name="Treatment",values=c(16,18,15),labels = c("Chronic", "Control","Intense"))+
  facet_grid(Site~Year, scales = 'free_y')+
  geom_point(size=3)+
  stat_ellipse(show.legend = FALSE)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), legend.position = "top",legend.title = element_blank(),
                   strip.text = element_text(size = 14))
NMDS_all
ggsave(filename = "NMDS_all.jpeg", plot = NMDS_all, bg = "transparent", width =  9, height = 6, units = "in", dpi = 600)

#sgs

NMDS_SGS<-ggplot(subset(SGS.mds_xy,Year %in% c("2013","2017","2021")), aes(x=MDS1, y=MDS3, color=Trt,shape=Trt)) + 
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Chronic", "Control","Intense"))+
  scale_shape_manual(name="Treatment",values=c(16,18,15),labels = c("Chronic", "Control","Intense"))+
  ggtitle("a) SGS")+
  facet_wrap(~Year, ncol = 3, scales = 'free_y')+
  geom_point(size=3)+
  stat_ellipse(show.legend = FALSE)+
  #scale_x_continuous(breaks = c(-1.25,0,1.25),limits = c(-1.30,1.25))+
  scale_y_continuous(breaks = c(-1.0,0,1.0),expand = c(0, 0), limits = c(-1.0,1.25))+
  theme_classic()+
  theme(legend.position = "none",legend.title = element_blank(),panel.grid.major = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        plot.title = element_text(size=20),
        axis.text.y=element_text(size=11, color="black"),
        axis.text.x=element_text(size=11, color="black"), 
        strip.background = element_rect(fill = "lightgrey"),strip.text = element_text(size = 16),
        panel.background = element_blank())
NMDS_SGS

#CHY
NMDS_CHY<-ggplot(subset(CHY.mds_xy,Year %in% c("2013","2017","2021")), aes(x=MDS1, y=MDS3, color=Trt,shape=Trt)) + 
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Chronic", "Control","Intense"))+
  scale_shape_manual(name="Treatment",values=c(16,18,15),labels = c("Chronic", "Control","Intense"))+
  ggtitle("b) HPG")+
  facet_wrap(~Year, ncol=3,scales = 'free_y')+
  geom_point(size=3)+
  stat_ellipse(show.legend = FALSE)+
  #scale_x_continuous(breaks = c(-1.25,0,1.25),limits = c(-1.25,1.25))+
  scale_y_continuous(breaks = c(-1.0,0,1.0),expand = c(0, 0), limits = c(-1.25,1.30))+  
  theme_classic()+
  theme(legend.position = "none",legend.title = element_blank(),panel.grid.major = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        plot.title = element_text(size=20),
        axis.text.y=element_text(size=11, color="black"),
        axis.text.x=element_text(size=11, color="black"),strip.background = element_blank(),strip.text = element_blank(),
        panel.background = element_blank())
NMDS_CHY

#HYS
NMDS_HYS<-ggplot(subset(HYS.mds_xy,Year %in% c("2013","2017","2021")), aes(x=MDS1, y=MDS3, color=Trt,shape=Trt)) + 
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Chronic", "Control","Intense"))+
  scale_shape_manual(name="Treatment",values=c(16,18,15),labels = c("Chronic", "Control","Intense"))+
  ggtitle("c) HYS")+
  facet_wrap(~Year, ncol=3, scales = 'free_y')+
  geom_point(size=3)+
  stat_ellipse(show.legend = FALSE)+
  #scale_x_continuous(breaks = c(-1.0,0,1.0),limits = c(-1.0,1.0))+
  scale_y_continuous(breaks = c(-1.0,0,1.0),expand = c(0, 0), limits = c(-1.25,1.30))+
  theme_classic()+
  theme(legend.position = "none",legend.title = element_blank(),panel.grid.major = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        plot.title = element_text(size=20),
        axis.text.y=element_text(size=11, color="black"),
        axis.text.x=element_text(size=11, color="black"),strip.background = element_blank(),strip.text = element_blank(),
        panel.background = element_blank())
NMDS_HYS

#KNZ
NMDS_KNZ<-ggplot(subset(KNZ.mds_xy,Year %in% c("2013","2017","2021")), aes(x=MDS1, y=MDS3, color=Trt,shape=Trt)) + 
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Chronic", "Control","Intense"))+
  scale_shape_manual(name="Treatment",values=c(16,18,15),labels = c("Chronic", "Control","Intense"))+
  ggtitle("d) KNZ")+
  facet_wrap(~Year, ncol=3, scales = 'free_y')+
  geom_point(size=3)+
  stat_ellipse(show.legend = FALSE)+
  #scale_x_continuous(breaks = c(-1.0,0,1.0),limits = c(-1.0,1.0))+
  scale_y_continuous(breaks = c(-1.0,0,1.0),expand = c(0, 0), limits = c(-1.25,1.30))+
  theme_classic()+
  theme(legend.position = "none",legend.title = element_blank(),panel.grid.major = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        plot.title = element_text(size=20),
        axis.text.y=element_text(size=11, color="black"),
        axis.text.x=element_text(size=11, color="black"),strip.background = element_blank(),strip.text = element_blank(),
        panel.background = element_blank())
NMDS_KNZ

NMDS_panel_2<-ggarrange(NMDS_SGS,NMDS_CHY,NMDS_HYS,NMDS_KNZ, ncol=1, nrow=4, heights = c(1.11,1,1,1))
NMDS_panel_2

ggsave(filename = "NMDS_panel_2.pdf", plot = NMDS_panel_2, bg = "transparent", width =  9, height = 6, units = "in", dpi = 600)

