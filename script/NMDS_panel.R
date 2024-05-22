
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





spcomp_long.final<-read.csv("data/data_filtered_names.csv")
spcomp_long.final$Trt <- factor(spcomp_long.final$Trt, levels=c('con', 'chr', 'int'))
spcomp_long.final$Species <- NULL


#Make the following items into a factor:
spcomp_long.final<-spcomp_long.final %>% 
  dplyr::mutate_at(c('Site','Year','Trt','Plot','Block'),as.factor)

#all years one model; not separated by period 
spcomp_long.final<- spcomp_long.final %>%
  dplyr::filter(Year %in% c(2013,2014,2015,2016,2017,2018,2019,2020,2021)) 


#separate by site for matrix 

spcomp_long.final_CHY <- spcomp_long.final %>% dplyr::filter(Site=="CHY")
spcomp_long.final_KNZ <- spcomp_long.final %>% dplyr::filter(Site=="KNZ")
spcomp_long.final_HYS <- spcomp_long.final %>% dplyr::filter(Site=="HYS")
spcomp_long.final_SGS <- spcomp_long.final %>% dplyr::filter(Site=="SGS")

#########################
#### VEGAN ANALYSES #####
#########################

CHY_wide<-spread(spcomp_long.final_CHY,Spcode,avg.cover)
KNZ_wide<-spread(spcomp_long.final_KNZ,Spcode,avg.cover)
HYS_wide<-spread(spcomp_long.final_HYS,Spcode,avg.cover)
SGS_wide<-spread(spcomp_long.final_SGS,Spcode,avg.cover)


#Change NAs to zero
CHY_wide[is.na(CHY_wide)] <- 0
KNZ_wide[is.na(KNZ_wide)] <- 0
HYS_wide[is.na(HYS_wide)] <- 0
SGS_wide[is.na(SGS_wide)] <- 0

str(CHY_wide)

#community matrix 

CHY.comm<-CHY_wide[,6:76] #Columns with species code
KNZ.comm<-KNZ_wide[,6:64] #Columns with species code
HYS.comm<-HYS_wide[,6:90]
SGS.comm<-SGS_wide[,6:49]


#PERMANOVA
adonis2(CHY.comm~CHY_wide$Year*CHY_wide$Trt, method="bray", permutations=1000, Blocks(CHY_wide$Block))
adonis2(KNZ.comm~KNZ_wide$Year*KNZ_wide$Trt, method="bray", permutations=1000, Blocks(KNZ_wide$Block))
adonis2(HYS.comm~HYS_wide$Year*HYS_wide$Trt, method="bray", permutations=1000, Blocks(HYS_wide$Block))
adonis2(SGS.comm~SGS_wide$Year*SGS_wide$Trt, method="bray", permutations=1000, Blocks(SGS_wide$Block))

factor.CHY<-CHY_wide[,1:5]
factor.KNZ<-KNZ_wide[,1:5]
factor.HYS<-HYS_wide[,1:5]
factor.SGS<-SGS_wide[,1:5]

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
                   panel.grid.minor = element_blank(), legend.position = "top",
                   strip.text = element_text(size = 14))
NMDS_all
ggsave(filename = "NMDS_all.jpeg", plot = NMDS_all, bg = "transparent", width =  9, height = 6, units = "in", dpi = 600)
