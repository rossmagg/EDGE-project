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

#Pairwise comparisons for each site, no interation 

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
mod.SGS <- pairwise.adonis2(SGS.comm.drt~comb, data=factor.drt.SGS, sim.method="bray",p.adjust.m = "bonferroni", permutations = perm)

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

# SIMPER BY PERIOD 
#SGS
(sim<-with(SGS_wide.drt,simper(SGS.comm.drt,Trt))) 
summary(sim,ordered=TRUE,digits=2)

(sim<-with(SGS_wide.rec,simper(SGS.comm.rec,Trt))) 
summary(sim,ordered=TRUE,digits=2)

#CHY 

(sim<-with(CHY_wide.drt,simper(CHY.comm.drt,Trt))) 
summary(sim,ordered=TRUE,digits=2)

(sim<-with(CHY_wide.rec,simper(CHY.comm.rec,Trt))) 
summary(sim,ordered=TRUE,digits=2)


#HYS
(sim<-with(HYS_wide.drt,simper(HYS.comm.drt,Trt))) 
summary(sim,ordered=TRUE,digits=2)

(sim<-with(HYS_wide.rec,simper(HYS.comm.rec,Trt))) 
summary(sim,ordered=TRUE,digits=2)


#KNZ
(sim<-with(KNZ_wide.drt,simper(KNZ.comm.drt,Trt))) 
summary(sim,ordered=TRUE,digits=2)

(sim<-with(KNZ_wide.rec,simper(KNZ.comm.rec,Trt))) 
summary(sim,ordered=TRUE,digits=2)


