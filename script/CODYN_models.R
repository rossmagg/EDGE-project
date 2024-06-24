##codyn mixed models##

library(ggplot2)
library(tidyverse)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)
library(Matrix)
library(codyn)
#RAC change

spcomp_long.final<-read.csv("data/spdata.plot_names.csv")

#Make the following items into a factor:
spcomp_long.final<-spcomp_long.final %>% 
  dplyr::mutate_at(c('Site','Spcode','Year','Plot','Block','Trt'),as.factor)


spcomp_long.final<- spcomp_long.final %>%
  dplyr::filter(Year %in% c(2013,2014,2015,2016,2017,2018,2019,2020,2021)) 

sites<-unique(spcomp_long.final$Site)
racchange<-data.frame() #Make dataframe called racchange

plotinfo<-spcomp_long.final %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Site,Plot,Trt,Block) %>% 
  unique ()


#need to put in console i=1 to test if it works
for(i in 1:length(sites)){
  sub<-spcomp_long.final %>% 
    filter(Site==sites[i])
  out<-RAC_change(sub, time.var = "Year", species.var = "Spcode", abundance.var = "avg.cover", replicate.var = "Plot", reference.time = 2013)
  out$Site<-sites[i]  
  racchange<-rbind(racchange,out)
}

racchangetrt<- merge(racchange, plotinfo, fix.by=c("Site","Plot"))
racchangetrt$Year.Year2<-str_c(racchangetrt$Year,'.',racchangetrt$Year2)

racchangetrt.drought <- racchangetrt %>% dplyr::filter(Year.Year2 %in% c("2013.2014","2013.2015","2013.2016","2013.2017"))

racchangtrt.recovery <- racchangetrt %>% dplyr::filter(Year.Year2 %in% c("2013.2018","2013.2019","2013.2020","2013.2021"))


### recovery compared to 2017 
spcomp_long.final17<- spcomp_long.final %>%
  dplyr::filter(Year %in% c(2017,2018,2019,2020,2021)) 

sites<-unique(spcomp_long.final$Site)
racchange17<-data.frame() #Make dataframe called racchange

plotinfo<-spcomp_long.final17 %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Site,Plot,Trt,Block) %>% 
  unique ()


#need to put in console i=1 to test if it works
for(i in 1:length(sites)){
  sub<-spcomp_long.final17 %>% 
    filter(Site==sites[i])
  out<-RAC_change(sub, time.var = "Year", species.var = "Spcode", 
                  abundance.var = "avg.cover", replicate.var = "Plot", 
                  reference.time = 2017)
  out$Site<-sites[i]  
  racchange17<-rbind(racchange17,out)
}

racchangetrt17<- merge(racchange17, plotinfo, fix.by=c("Site","Plot"))

#Create a new column showing the year comparisons
racchangetrt17$Year.Year2<-str_c(racchangetrt17$Year,'.',racchangetrt17$Year2)

###### RAC change MODELS #####
#### SGS ####
#rich full
RiC.SGS<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="SGS"))
plot(RiC.SGS)
anova(RiC.SGS, ddf="Kenward-Roger")
emmeans(RiC.SGS, pairwise ~Trt|Year.Year2)

#rich drought
Dr.RiC.SGS<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="SGS"))
plot(Dr.RiC.SGS)
anova(Dr.RiC.SGS, ddf="Kenward-Roger")
emmeans(Dr.RiC.SGS, pairwise ~Trt|Year.Year2)

#rich recovery-2013
Re.RiC.SGS<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="SGS"))
plot(Re.RiC.SGS)
anova(Re.RiC.SGS, ddf="Kenward-Roger")

#rich recovery-2017
Re17.RiC.SGS<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="SGS"))
plot(Re17.RiC.SGS)
anova(Re17.RiC.SGS, ddf="Kenward-Roger")


#even full 
EvC.SGS<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="SGS"))
plot(EvC.SGS)
anova(EvC.SGS, ddf="Kenward-Roger")
emmeans(EvC.SGS, pairwise ~Trt|Year.Year2)

#even drought
Dr.EvC.SGS<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="SGS"))
plot(Dr.EvC.SGS)
anova(Dr.EvC.SGS, ddf="Kenward-Roger")
emmeans(Dr.EvC.SGS, pairwise ~Trt|Year.Year2)

#even recovery-2013
Re.EvC.SGS<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="SGS"))
plot(Re.EvC.SGS)
anova(Re.EvC.SGS, ddf="Kenward-Roger")
emmeans(Re.EvC.SGS, pairwise ~Trt|Year.Year2)

#even recovery-2017
Re17.EvC.SGS<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="SGS"))
plot(Re17.EvC.SGS)
anova(Re17.EvC.SGS, ddf="Kenward-Roger")
emmeans(Re17.EvC.SGS, pairwise ~Trt|Year.Year2)

#rank change full
RaC.SGS<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="SGS"))
plot(RaC.SGS)
anova(RaC.SGS, ddf="Kenward-Roger")
emmeans(RaC.SGS, pairwise ~Trt|Year.Year2)

#rank drought
Dr.RaC.SGS<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="SGS"))
plot(Dr.RaC.SGS)
anova(Dr.RaC.SGS, ddf="Kenward-Roger")

#rank  recovery-2013
Re.RaC.SGS<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="SGS"))
plot(Re.RaC.SGS)
anova(Re.RaC.SGS, ddf="Kenward-Roger")
emmeans(Re.RaC.SGS, pairwise ~Trt|Year.Year2)

#rank recovery-2017
Re17.RaC.SGS<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="SGS"))
plot(Re17.RaC.SGS)
anova(Re17.RaC.SGS, ddf="Kenward-Roger")
emmeans(Re17.RaC.SGS, pairwise ~Trt|Year.Year2)

#GAINS full
GaC.SGS<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="SGS"))
plot(GaC.SGS)
anova(GaC.SGS, ddf="Kenward-Roger")

#gains drought
Dr.GaC.SGS<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="SGS"))
plot(Dr.GaC.SGS)
anova(Dr.GaC.SGS, ddf="Kenward-Roger")

#gains recovery-2013
Re.GaC.SGS<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="SGS"))
plot(Re.GaC.SGS)
anova(Re.GaC.SGS, ddf="Kenward-Roger")

#gains recovery-2017
Re17.GaC.SGS<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="SGS"))
plot(Re17.GaC.SGS)
anova(Re17.GaC.SGS, ddf="Kenward-Roger")

#LOSSES
LoC.SGS<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="SGS"))
plot(LoC.SGS)
anova(LoC.SGS, ddf="Kenward-Roger")
emmeans(LoC.SGS, pairwise ~Trt|Year.Year2)

#losses drought
Dr.LoC.SGS<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="SGS"))
plot(Dr.LoC.SGS)
anova(Dr.LoC.SGS, ddf="Kenward-Roger")
emmeans(Dr.LoC.SGS, pairwise ~Trt|Year.Year2)

#losses recovery-2013
Re.LoC.SGS<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="SGS"))
plot(Re.LoC.SGS)
anova(Re.LoC.SGS, ddf="Kenward-Roger")
emmeans(Re.LoC.SGS, pairwise ~Trt|Year.Year2)

#losses recovery-2017
Re17.LoC.SGS<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="SGS"))
plot(Re17.LoC.SGS)
anova(Re17.LoC.SGS, ddf="Kenward-Roger")

#### CHY ####
#rich full
RiC.CHY<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="CHY"))
plot(RiC.CHY)
anova(RiC.CHY, ddf="Kenward-Roger")

#rich drought
Dr.RiC.CHY<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="CHY"))
plot(Dr.RiC.CHY)
anova(Dr.RiC.CHY, ddf="Kenward-Roger")

#rich recovery-2013
Re.RiC.CHY<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="CHY"))
plot(Re.RiC.CHY)
anova(Re.RiC.CHY, ddf="Kenward-Roger")

#rich recovery-2017
Re17.RiC.CHY<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="CHY"))
plot(Re17.RiC.CHY)
anova(Re17.RiC.CHY, ddf="Kenward-Roger")

#even full 
EvC.CHY<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="CHY"))
plot(EvC.CHY)
anova(EvC.CHY, ddf="Kenward-Roger")
emmeans(EvC.CHY, pairwise ~Trt|Year.Year2)

#even drought
Dr.EvC.CHY<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="CHY"))
plot(Dr.EvC.CHY)
anova(Dr.EvC.CHY, ddf="Kenward-Roger")

#even recovery-2013
Re.EvC.CHY<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="CHY"))
plot(Re.EvC.CHY)
anova(Re.EvC.CHY, ddf="Kenward-Roger")
emmeans(Re.EvC.CHY, pairwise ~Trt|Year.Year2)

#even recovery-2017
Re17.EvC.CHY<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="CHY"))
plot(Re17.EvC.CHY)
anova(Re17.EvC.CHY, ddf="Kenward-Roger")
emmeans(Re17.EvC.CHY, pairwise ~Trt|Year.Year2)

#rank change full
RaC.CHY<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="CHY"))
plot(RaC.CHY)
anova(RaC.CHY, ddf="Kenward-Roger")
emmeans(RaC.CHY, pairwise ~Trt|Year.Year2)

#rank drought
Dr.RaC.CHY<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="CHY"))
plot(Dr.RaC.CHY)
anova(Dr.RaC.CHY, ddf="Kenward-Roger")

#rank  recovery-2013
Re.RaC.CHY<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="CHY"))
plot(Re.RaC.CHY)
anova(Re.RaC.CHY, ddf="Kenward-Roger")
emmeans(Re.RaC.CHY, pairwise ~Trt|Year.Year2)

#rank recovery-2017
Re17.RaC.CHY<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="CHY"))
plot(Re17.RaC.CHY)
anova(Re17.RaC.CHY, ddf="Kenward-Roger")
emmeans(Re17.RaC.CHY, pairwise ~Trt|Year.Year2)

#GAINS full
GaC.CHY<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="CHY"))
plot(GaC.CHY)
anova(GaC.CHY, ddf="Kenward-Roger")

#gains drought
Dr.GaC.CHY<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="CHY"))
plot(Dr.GaC.CHY)
anova(Dr.GaC.CHY, ddf="Kenward-Roger")

#gains recovery-2013
Re.GaC.CHY<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="CHY"))
plot(Re.GaC.CHY)
anova(Re.GaC.CHY, ddf="Kenward-Roger")

#gains recovery-2017
Re17.GaC.CHY<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="CHY"))
plot(Re17.GaC.CHY)
anova(Re17.GaC.CHY, ddf="Kenward-Roger")

#LOSSES
LoC.CHY<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="CHY"))
plot(LoC.CHY)
anova(LoC.CHY, ddf="Kenward-Roger")

#losses drought
Dr.LoC.CHY<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="CHY"))
plot(Dr.LoC.CHY)
anova(Dr.LoC.CHY, ddf="Kenward-Roger")

#losses recovery-2013
Re.LoC.CHY<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="CHY"))
plot(Re.LoC.CHY)
anova(Re.LoC.CHY, ddf="Kenward-Roger")

#losses recovery-2017
Re17.LoC.CHY<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="CHY"))
plot(Re17.LoC.CHY)
anova(Re17.LoC.CHY, ddf="Kenward-Roger")

#### HYS ####
#rich full
RiC.HYS<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="HYS"))
plot(RiC.HYS)
anova(RiC.HYS, ddf="Kenward-Roger")
emmeans(RiC.HYS, pairwise ~Trt|Year.Year2)

#rich drought
Dr.RiC.HYS<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="HYS"))
plot(Dr.RiC.HYS)
anova(Dr.RiC.HYS, ddf="Kenward-Roger")
emmeans(Dr.RiC.HYS, pairwise ~Trt|Year.Year2)

#rich recovery-2013
Re.RiC.HYS<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="HYS"))
plot(Re.RiC.HYS)
anova(Re.RiC.HYS, ddf="Kenward-Roger")
emmeans(Re.RiC.HYS, pairwise ~Trt|Year.Year2)

#rich recovery-2017
Re17.RiC.HYS<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="HYS"))
plot(Re17.RiC.HYS)
anova(Re17.RiC.HYS, ddf="Kenward-Roger")
emmeans(Re17.RiC.HYS, pairwise ~Trt|Year.Year2)

#even full 
EvC.HYS<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="HYS"))
plot(EvC.HYS)
anova(EvC.HYS, ddf="Kenward-Roger")
emmeans(EvC.HYS, pairwise ~Trt|Year.Year2)

#even drought
Dr.EvC.HYS<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="HYS"))
plot(Dr.EvC.HYS)
anova(Dr.EvC.HYS, ddf="Kenward-Roger")
emmeans(Dr.EvC.HYS, pairwise ~Trt|Year.Year2)

#even recovery-2013
Re.EvC.HYS<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="HYS"))
plot(Re.EvC.HYS)
anova(Re.EvC.HYS, ddf="Kenward-Roger")
emmeans(Re.EvC.HYS, pairwise ~Trt|Year.Year2)

#even recovery-2017
Re17.EvC.HYS<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="HYS"))
plot(Re17.EvC.HYS)
anova(Re17.EvC.HYS, ddf="Kenward-Roger")


#rank change full
RaC.HYS<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="HYS"))
plot(RaC.HYS)
anova(RaC.HYS, ddf="Kenward-Roger")
emmeans(RaC.HYS, pairwise ~Trt|Year.Year2)

#rank drought
Dr.RaC.HYS<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="HYS"))
plot(Dr.RaC.HYS)
anova(Dr.RaC.HYS, ddf="Kenward-Roger")
emmeans(Dr.RaC.HYS, pairwise ~Trt|Year.Year2)

#rank  recovery-2013
Re.RaC.HYS<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="HYS"))
plot(Re.RaC.HYS)
anova(Re.RaC.HYS, ddf="Kenward-Roger")
emmeans(Re.RaC.HYS, pairwise ~Trt|Year.Year2)

#rank recovery-2017
Re17.RaC.HYS<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="HYS"))
plot(Re17.RaC.HYS)
anova(Re17.RaC.HYS, ddf="Kenward-Roger")
emmeans(Re17.RaC.HYS, pairwise ~Trt|Year.Year2)

#GAINS full
GaC.HYS<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="HYS"))
plot(GaC.HYS)
anova(GaC.HYS, ddf="Kenward-Roger")

#gains drought
Dr.GaC.HYS<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="HYS"))
plot(Dr.GaC.HYS)
anova(Dr.GaC.HYS, ddf="Kenward-Roger")

#gains recovery-2013
Re.GaC.HYS<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="HYS"))
plot(Re.GaC.HYS)
anova(Re.GaC.HYS, ddf="Kenward-Roger")

#gains recovery-2017
Re17.GaC.HYS<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="HYS"))
plot(Re17.GaC.HYS)
anova(Re17.GaC.HYS, ddf="Kenward-Roger")

#LOSSES
LoC.HYS<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="HYS"))
plot(LoC.HYS)
anova(LoC.HYS, ddf="Kenward-Roger")


#losses drought
Dr.LoC.HYS<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="HYS"))
plot(Dr.LoC.HYS)
anova(Dr.LoC.HYS, ddf="Kenward-Roger")
emmeans(Dr.LoC.HYS, pairwise ~Trt|Year.Year2)

#losses recovery-2013
Re.LoC.HYS<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="HYS"))
plot(Re.LoC.HYS)
anova(Re.LoC.HYS, ddf="Kenward-Roger")
emmeans(Re.LoC.HYS, pairwise ~Trt|Year.Year2)

#losses recovery-2017
Re17.LoC.HYS<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="HYS"))
plot(Re17.LoC.HYS)
anova(Re17.LoC.HYS, ddf="Kenward-Roger")
emmeans(Re17.LoC.HYS, pairwise ~Trt|Year.Year2)

#### KNZ ####
#rich full
RiC.KNZ<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="KNZ"))
plot(RiC.KNZ)
anova(RiC.KNZ, ddf="Kenward-Roger")

#rich drought
Dr.RiC.KNZ<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="KNZ"))
plot(Dr.RiC.KNZ)
anova(Dr.RiC.KNZ, ddf="Kenward-Roger")

#rich recovery-2013
Re.RiC.KNZ<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="KNZ"))
plot(Re.RiC.KNZ)
anova(Re.RiC.KNZ, ddf="Kenward-Roger")

#rich recovery-2017
Re17.RiC.KNZ<- lmer(richness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="KNZ"))
plot(Re17.RiC.KNZ)
anova(Re17.RiC.KNZ, ddf="Kenward-Roger")

#even full 
EvC.KNZ<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="KNZ"))
plot(EvC.KNZ)
anova(EvC.KNZ, ddf="Kenward-Roger")

#even drought
Dr.EvC.KNZ<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="KNZ"))
plot(Dr.EvC.KNZ)
anova(Dr.EvC.KNZ, ddf="Kenward-Roger")

#even recovery-2013
Re.EvC.KNZ<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="KNZ"))
plot(Re.EvC.KNZ)
anova(Re.EvC.KNZ, ddf="Kenward-Roger")

#even recovery-2017
Re17.EvC.KNZ<- lmer(evenness_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="KNZ"))
plot(Re17.EvC.KNZ)
anova(Re17.EvC.KNZ, ddf="Kenward-Roger")

#rank change full
RaC.KNZ<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="KNZ"))
plot(RaC.KNZ)
anova(RaC.KNZ, ddf="Kenward-Roger")

#rank drought
Dr.RaC.KNZ<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="KNZ"))
plot(Dr.RaC.KNZ)
anova(Dr.RaC.KNZ, ddf="Kenward-Roger")

#rank  recovery-2013
Re.RaC.KNZ<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="KNZ"))
plot(Re.RaC.KNZ)
anova(Re.RaC.KNZ, ddf="Kenward-Roger")

#rank recovery-2017
Re17.RaC.KNZ<- lmer(rank_change~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="KNZ"))
plot(Re17.RaC.KNZ)
anova(Re17.RaC.KNZ, ddf="Kenward-Roger")
emmeans(Re17.RaC.KNZ, pairwise ~Trt|Year.Year2)

#GAINS full
GaC.KNZ<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="KNZ"))
plot(GaC.KNZ)
anova(GaC.KNZ, ddf="Kenward-Roger")
emmeans(GaC.KNZ, pairwise ~Trt|Year.Year2)

#gains drought
Dr.GaC.KNZ<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="KNZ"))
plot(Dr.GaC.KNZ)
anova(Dr.GaC.KNZ, ddf="Kenward-Roger")

#gains recovery-2013
Re.GaC.KNZ<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="KNZ"))
plot(Re.GaC.KNZ)
anova(Re.GaC.KNZ, ddf="Kenward-Roger")

#gains recovery-2017
Re17.GaC.KNZ<- lmer(gains~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="KNZ"))
plot(Re17.GaC.KNZ)
anova(Re17.GaC.KNZ, ddf="Kenward-Roger")
emmeans(Re17.GaC.KNZ, pairwise ~Trt|Year.Year2)

#LOSSES
LoC.KNZ<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt, Site=="KNZ"))
plot(LoC.KNZ)
anova(LoC.KNZ, ddf="Kenward-Roger")

#losses drought
Dr.LoC.KNZ<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt.drought, Site=="KNZ"))
plot(Dr.LoC.KNZ)
anova(Dr.LoC.KNZ, ddf="Kenward-Roger")

#losses recovery-2013
Re.LoC.KNZ<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangtrt.recovery, Site=="KNZ"))
plot(Re.LoC.KNZ)
anova(Re.LoC.KNZ, ddf="Kenward-Roger")

#losses recovery-2017
Re17.LoC.KNZ<- lmer(losses~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(racchangetrt17, Site=="KNZ"))
plot(Re17.LoC.KNZ)
anova(Re17.LoC.KNZ, ddf="Kenward-Roger")

########### Abundance change models ##############

comp_zeros <- read.csv("data/spcomp_zeros_trt.csv")
comp_zeros$Species = NULL
colnames(comp_zeros)[colnames(comp_zeros) == "max.cover"] <- "cover" 
comp_zeros<- comp_zeros %>%
  dplyr::filter(Year %in% c(2013,2014,2015,2016,2017,2018,2019,2020,2021)) 

comp.wide <- comp_zeros %>%
  pivot_wider(names_from = Year, values_from = cover)
comp.wide[is.na(comp.wide)] <- 0

comp.zeros_final <- comp.wide %>%
  pivot_longer(cols = "2013":"2020",names_to = "Year", values_to = "cover")

comp.aveplot <- comp.zeros_final %>%
  dplyr::group_by(Site,Year,Plot,Spcode) %>% 
  dplyr::summarise(avg.cover = mean(cover))%>%
  dplyr::as_tibble()

comp.wide_plot <- comp.aveplot %>%
  pivot_wider(names_from = Plot, values_from = avg.cover)
comp.wide_plot[is.na(comp.wide_plot)] <- 0
comp.aveplot_zeros<- comp.wide_plot %>%
  pivot_longer(cols = "1":"30",names_to = "Plot", values_to = "avg.cover")

trt <- read.csv("data/trt_info_2.csv")
trt<-trt %>% 
  dplyr::mutate_at(c('Site','Plot','Block','Trt'),as.factor)
comp.aveplot_trt<- merge(trt,comp.aveplot_zeros, fix.by=c("Site", "Plot"))

comp.species <- comp.aveplot_trt %>% dplyr::filter(Spcode %in% c("45","46","47","48","49","55","130","137","197","317","95","44", "9","7","89","221","235","234", "204","50","53","194","51","14","88","96","388","94","56","57","325"))
comp.species$Spcode <- factor(comp.species$Spcode, levels = c("45","46","47","48","49","55","130","137","197","317","95","44","9","7","89","221","235","234", "204","50","53","194","51","14","88","96","388","94","56","57","325"), 
                              labels = c("BODA","BOGR","BOHI","BRJA","BRTE","CAEL","HECO","KOPY","PASM","VUOC","ELEL","BOCU","ANTE","ANGE","DIOL","SCSC","SPAS","SONU","POPR","CABR","CAME","PAVI","CAHE","ARPU","DIAC","ERSP","PACA","ELSP","CAFI","CAGR","CABL"))

comp.path <- comp.species %>%
  mutate(Path = case_when(
    Spcode == "BODA" ~ "C4",
    Spcode == "BOGR" ~ "C4",
    Spcode == "VUOC" ~ "C3",
    Spcode == "HECO" ~ "C3",
    Spcode == "BRTE" ~ "C3",
    Spcode == "CAEL" ~ "C3", 
    Spcode == "KOPY" ~ "C3",
    Spcode == "PASM" ~ "C3",
    Spcode == "BOCU" ~ "C4",
    Spcode == "ANTE" ~ "C4",
    Spcode == "ANGE" ~ "C4",
    Spcode == "SPAS" ~ "C4", 
    Spcode == "SONU" ~ "C4",
    Spcode == "SCSC" ~ "C4",
    Spcode == "DIOL" ~ "C3",
    Spcode == "POPR" ~ "C3",
    Spcode == "PAVI" ~ "C4",
    Spcode == "CAME" ~ "C3",
    Spcode == "CABR" ~ "C3",
    Spcode == "CAHE" ~ "C3",
    Spcode == "ELEL" ~ "C3",
    Spcode == "BOHI" ~ "C4",
    Spcode == "BRJA" ~ "C3",
    Spcode == "ARPU" ~ "C4", 
    Spcode == "DIAC" ~ "C3",
    Spcode == "ERSP" ~ "C4",
    Spcode == "PACA" ~ "C4",
    Spcode == "ELSP" ~ "C3",
    Spcode == "CAFI" ~ "C3",
    Spcode == "CAGR" ~ "C3",
    Spcode == "CABL" ~ "C3"))

comppath.plotsum <- comp.path %>%
  dplyr::group_by(Site,Year,Block,Trt,Plot,Path) %>% 
  dplyr::summarise(pathsum = sum(avg.cover))

comppath.plotsum$Trt <- factor(comppath.plotsum$Trt, levels=c('con', 'chr', 'int'))
comppath.plotsum$Site <- factor(comppath.plotsum$Site, levels=c('SGS', 'CHY', 'HYS','KNZ'))

sites<-unique(comppath.plotsum$Site)
abunchng<-data.frame() #Make dataframe called abunchng

#Making a dataset with treatment info
plotinfo<-comppath.plotsum %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Site,Plot, Trt,Block) %>% 
  unique()

#need to put in console i=1 to test if it works
for(i in 1:length(sites)){
  sub<-comppath.plotsum %>% 
    filter(Site==sites[i])
  out<-abundance_change(
    df=sub,
    time.var="Year",
    species.var="Path",
    abundance.var="pathsum",
    replicate.var = "Plot",
    reference.time = "2013"
  )
  out$Site<-sites[i]  
  abunchng<-rbind(abunchng,out)
}

#Merge dataframe with trt info
abunchngtrt<-abunchng %>% 
  left_join(plotinfo)

#Create a new column showing the year comparisons
abunchngtrt$Year.Year2<-str_c(abunchngtrt$Year,'-',abunchngtrt$Year2)

#change to wide 
abunchngtrt_wide <- abunchngtrt %>%
  pivot_wider(names_from = Path, values_from = change)

abunchngtrt_wide.drought <- abunchngtrt_wide %>% dplyr::filter(Year.Year2 %in% c("2013-2014","2013-2015","2013-2016","2013-2017"))

abunchngtrt_wide.recovery <- abunchngtrt_wide %>% dplyr::filter(Year.Year2 %in% c("2013-2018","2013-2019","2013-2020","2013-2021"))

#### recovery change compared to 2017 ####
comppath.plotsum17 <- comppath.plotsum %>% dplyr::filter(Year %in% c("2017","2018","2019","2020","2021"))

abunchng17<-data.frame()
#need to put in console i=1 to test if it works
for(i in 1:length(sites)){
  sub<-comppath.plotsum17 %>% 
    filter(Site==sites[i])
  out<-abundance_change(
    df=sub,
    time.var="Year",
    species.var="Path",
    abundance.var="pathsum",
    replicate.var = "Plot",
    reference.time = "2017"
  )
  out$Site<-sites[i]  
  abunchng17<-rbind(abunchng17,out)
}

#Merge dataframe with trt info
abunchngtrt17<-abunchng17 %>% 
  left_join(plotinfo)

#Create a new column showing the year comparisons
abunchngtrt17$Year.Year2<-str_c(abunchngtrt17$Year,'-',abunchngtrt17$Year2)

#change to wide 
abunchngtrt_wide17 <- abunchngtrt17 %>%
  pivot_wider(names_from = Path, values_from = change)

##C3 MODELS
#SGS
#full
SGS.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide, Site=="SGS"))
plot(SGS.C3)
anova(SGS.C3, ddf="Kenward-Roger")
emmeans(SGS.C3, pairwise ~Trt|Year.Year2)

#drought
Dr.SGS.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.drought, Site=="SGS"))
plot(Dr.SGS.C3)
anova(Dr.SGS.C3, ddf="Kenward-Roger")

#recov-2013
Re.SGS.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.recovery, Site=="SGS"))
plot(Re.SGS.C3)
anova(Re.SGS.C3, ddf="Kenward-Roger")
emmeans(Re.SGS.C3, pairwise ~Trt|Year.Year2)

#recov-2017

Re17.SGS.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide17, Site=="SGS"))
plot(Re17.SGS.C3)
anova(Re17.SGS.C3, ddf="Kenward-Roger")
emmeans(Re17.SGS.C3, pairwise ~Trt|Year.Year2)

#CHY
#full
CHY.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide, Site=="CHY"))
plot(CHY.C3)
anova(CHY.C3, ddf="Kenward-Roger")
emmeans(CHY.C3, pairwise ~Trt|Year.Year2)

#drought
Dr.CHY.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.drought, Site=="CHY"))
plot(Dr.CHY.C3)
anova(Dr.CHY.C3, ddf="Kenward-Roger")
emmeans(Dr.CHY.C3, pairwise ~Trt|Year.Year2)

#recov-2013
Re.CHY.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.recovery, Site=="CHY"))
plot(Re.CHY.C3)
anova(Re.CHY.C3, ddf="Kenward-Roger")

#recov-2017

Re17.CHY.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide17, Site=="CHY"))
plot(Re17.CHY.C3)
anova(Re17.CHY.C3, ddf="Kenward-Roger")

#HYS
#full
HYS.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide, Site=="HYS"))
plot(HYS.C3)
anova(HYS.C3, ddf="Kenward-Roger")
emmeans(HYS.C3, pairwise ~Trt|Year.Year2)

#drought
Dr.HYS.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.drought, Site=="HYS"))
plot(Dr.HYS.C3)
anova(Dr.HYS.C3, ddf="Kenward-Roger")
emmeans(Dr.HYS.C3, pairwise ~Trt|Year.Year2)

#recov-2013
Re.HYS.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.recovery, Site=="HYS"))
plot(Re.HYS.C3)
anova(Re.HYS.C3, ddf="Kenward-Roger")

#recov-2017

Re17.HYS.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide17, Site=="HYS"))
plot(Re17.HYS.C3)
anova(Re17.HYS.C3, ddf="Kenward-Roger")
emmeans(Re17.HYS.C3, pairwise ~Trt|Year.Year2)

#KNZ
#full
KNZ.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide, Site=="KNZ"))
plot(KNZ.C3)
anova(KNZ.C3, ddf="Kenward-Roger")
emmeans(KNZ.C3, pairwise ~Trt|Year.Year2)

#drought
Dr.KNZ.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.drought, Site=="KNZ"))
plot(Dr.KNZ.C3)
anova(Dr.KNZ.C3, ddf="Kenward-Roger")

#recov-2013
Re.KNZ.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.recovery, Site=="KNZ"))
plot(Re.KNZ.C3)
anova(Re.KNZ.C3, ddf="Kenward-Roger")
emmeans(Re.KNZ.C3, pairwise ~Trt|Year.Year2)

#recov-2017

Re17.KNZ.C3<- lmer(C3~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide17, Site=="KNZ"))
plot(Re17.KNZ.C3)
anova(Re17.KNZ.C3, ddf="Kenward-Roger")
emmeans(Re17.KNZ.C3, pairwise ~Trt|Year.Year2)

##C4 MODELS
#SGS
#SGS
#full
SGS.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide, Site=="SGS"))
plot(SGS.C4)
anova(SGS.C4, ddf="Kenward-Roger")
emmeans(SGS.C4, pairwise ~Trt|Year.Year2)

#drought
Dr.SGS.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.drought, Site=="SGS"))
plot(Dr.SGS.C4)
anova(Dr.SGS.C4, ddf="Kenward-Roger")
emmeans(Dr.SGS.C4, pairwise ~Trt|Year.Year2)

#recov-2013
Re.SGS.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.recovery, Site=="SGS"))
plot(Re.SGS.C4)
anova(Re.SGS.C4, ddf="Kenward-Roger")
emmeans(Re.SGS.C4, pairwise ~Trt|Year.Year2)

#recov-2017

Re17.SGS.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide17, Site=="SGS"))
plot(Re17.SGS.C4)
anova(Re17.SGS.C4, ddf="Kenward-Roger")
emmeans(Re17.SGS.C4, pairwise ~Trt|Year.Year2)

#CHY
#full
CHY.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide, Site=="CHY"))
plot(CHY.C4)
anova(CHY.C4, ddf="Kenward-Roger")
emmeans(CHY.C4, pairwise ~Trt|Year.Year2)

#drought
Dr.CHY.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.drought, Site=="CHY"))
plot(Dr.CHY.C4)
anova(Dr.CHY.C4, ddf="Kenward-Roger")
emmeans(Dr.CHY.C4, pairwise ~Trt|Year.Year2)

#recov-2013
Re.CHY.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.recovery, Site=="CHY"))
plot(Re.CHY.C4)
anova(Re.CHY.C4, ddf="Kenward-Roger")
emmeans(Re.CHY.C4, pairwise ~Trt|Year.Year2)

#recov-2017

Re17.CHY.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide17, Site=="CHY"))
plot(Re17.CHY.C4)
anova(Re17.CHY.C4, ddf="Kenward-Roger")
emmeans(Re17.CHY.C4, pairwise ~Trt|Year.Year2)

#HYS
#full
HYS.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide, Site=="HYS"))
plot(HYS.C4)
anova(HYS.C4, ddf="Kenward-Roger")
emmeans(HYS.C4, pairwise ~Trt|Year.Year2)

#drought
Dr.HYS.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.drought, Site=="HYS"))
plot(Dr.HYS.C4)
anova(Dr.HYS.C4, ddf="Kenward-Roger")
emmeans(Dr.HYS.C4, pairwise ~Trt|Year.Year2)

#recov-2013
Re.HYS.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.recovery, Site=="HYS"))
plot(Re.HYS.C4)
anova(Re.HYS.C4, ddf="Kenward-Roger")

#recov-2017

Re17.HYS.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide17, Site=="HYS"))
plot(Re17.HYS.C4)
anova(Re17.HYS.C4, ddf="Kenward-Roger")
emmeans(Re17.HYS.C4, pairwise ~Trt|Year.Year2)

#KNZ
#full
KNZ.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide, Site=="KNZ"))
plot(KNZ.C4)
anova(KNZ.C4, ddf="Kenward-Roger")
emmeans(KNZ.C4, pairwise ~Trt|Year.Year2)

#drought
Dr.KNZ.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.drought, Site=="KNZ"))
plot(Dr.KNZ.C4)
anova(Dr.KNZ.C4, ddf="Kenward-Roger")
emmeans(Dr.KNZ.C4, pairwise ~Trt|Year.Year2)

#recov-2013
Re.KNZ.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide.recovery, Site=="KNZ"))
plot(Re.KNZ.C4)
anova(Re.KNZ.C4, ddf="Kenward-Roger")
emmeans(Re.KNZ.C4, pairwise ~Trt|Year.Year2)

#recov-2017

Re17.KNZ.C4<- lmer(C4~Trt*Year.Year2+(1|Plot)+(1|Block), data = subset(abunchngtrt_wide17, Site=="KNZ"))
plot(Re17.KNZ.C4)
anova(Re17.KNZ.C4, ddf="Kenward-Roger")
emmeans(Re17.KNZ.C4, pairwise ~Trt|Year.Year2)
