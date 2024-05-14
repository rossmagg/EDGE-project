#################################################
## RAC and Abundance change averages by period ##
#################################################

library(dplyr)
library(tidyr)
library(data.table)
library(tidyverse)
library(plyr)
library(codyn)
library(vegan)
library(lme4)

#Other packages to load
library(emmeans)
library(lmerTest)
library(MuMIn)
library(AICcmodavg)

spcomp_long.final<-read.csv("data/spdata.plot_names.csv")

#Make the following items into a factor:
spcomp_long.final<-spcomp_long.final %>% 
  dplyr::mutate_at(c('Site','Spcode','Year','Plot','Block','Trt'),as.factor)

spcomp_long.final$Site <- factor(spcomp_long.final$Site, levels=c('SGS', 'CHY', 'HYS','KNZ'))
spcomp_long.final$Trt <- factor(spcomp_long.final$Trt, levels=c('con', 'chr', 'int'))


spcomp_long.final.drt<- spcomp_long.final %>%
  dplyr::filter(Year %in% c(2013,2014,2015,2016,2017)) 

spcomp_long.final.rec<- spcomp_long.final %>%
  dplyr::filter(Year %in% c(2017,2018,2019,2020,2021)) 

### RAC CHANGE during drought; compared to 2013 ###

sites<-unique(spcomp_long.final.drt$Site)
racchange<-data.frame() #Make dataframe called racchange

plotinfo<-spcomp_long.final.drt %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Site,Plot,Trt,Block) %>% 
  unique ()


#need to put in console i=1 to test if it works
for(i in 1:length(sites)){
  sub<-spcomp_long.final.drt %>% 
    filter(Site==sites[i])
  out<-RAC_change(sub, time.var = "Year", species.var = "Spcode", 
                  abundance.var = "avg.cover", replicate.var = "Plot", reference.time = 2013)
  out$Site<-sites[i]  
  racchange<-rbind(racchange,out)
}

racchangetrt<- merge(racchange, plotinfo, fix.by=c("Site","Plot"))

#Create a new column showing the year comparisons
racchangetrt$Year.Year2<-str_c(racchangetrt$Year,'.',racchangetrt$Year2)

racchangetrt.long<-racchangetrt %>% gather(Change.type,Value, richness_change:losses)

### Abundace change drought ###

#first calculate average cover and then convert to average cover of c3 and c4 graminoids 
comp_zeros <- read.csv("data/spcomp_zeros_trt.csv")
comp_zeros<-comp_zeros%>% 
  dplyr::mutate_at(c('Site','Spcode','Year','Plot','Block','Trt'),as.factor)
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

comppath.plotsum.drt<- comppath.plotsum %>%
  dplyr::filter(Year %in% c(2013,2014,2015,2016,2017)) 

comppath.plotsum.rec<- comppath.plotsum %>%
  dplyr::filter(Year %in% c(2017,2018,2019,2020,2021)) 

## calculate abundance change of total graminoid C3 and C4 for drought period 

sites<-unique(comppath.plotsum.drt$Site)
abunchng<-data.frame() #Make dataframe called abunchng

#Making a dataset with treatment info
plotinfo<-comppath.plotsum.drt %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Site,Plot, Trt) %>% 
  unique()

#need to put in console i=1 to test if it works
for(i in 1:length(sites)){
  sub<-comppath.plotsum.drt %>% 
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

#rename columns so they match with racchange

names(abunchngtrt)[names(abunchngtrt) == "Path"] <- "Change.type"
names(abunchngtrt)[names(abunchngtrt) == "change"] <- "Value"
racchangetrt.long$Block <- NULL

str(racchangetrt.long)

#merge dfs 
RAC_abund <- rbind(racchangetrt.long, abunchngtrt)
RAC_abund$Change.type<-as.factor(RAC_abund$Change.type)
str(RAC_abund)

#calculate aves 

RAC_abund_ave<-RAC_abund%>%
  dplyr::group_by(Site,Trt,Change.type)%>%
  dplyr::summarize(mean_change = mean(Value,na.rm=T), n = n(),sd = sd(Value,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_change - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_change + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

RAC_abund_ave$Trt <- factor(RAC_abund_ave$Trt, levels=c('con', 'chr', 'int'))
RAC_abund_ave$Change.type <- factor(RAC_abund_ave$Change.type, levels=c('richness_change','evenness_change', 'rank_change',"gains","losses","C3","C4"), 
                                    labels = c("Richness Chg.","Evenness Chg.","Reordering","Sp. Gains","Sp. Losses","C3 Abund. Chg.","C4 Abund. Chg."))

#figure
RAC_abund_drt_fig<-ggplot(RAC_abund_ave, aes(x=Site, y=mean_change,color=Trt))+
  ylab("Average change during drought")+
  xlab("Site")+
  facet_wrap(~Change.type, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3.5,position=position_dodge(.60))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.60), width=0, size=1.25, show.legend = TRUE)+
  scale_x_discrete(labels=c('SGS', 'HPG', 'HYS','KNZ'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                        axis.text.y = element_text(size=12, color = "black"),
                        strip.text = element_text(size=14),
                        axis.title.x = element_text(size=14),
                        axis.title.y = element_text(size=14),
                        legend.position = "top",legend.text = element_text(size=11))
RAC_abund_drt_fig
ggsave(filename = "RAC_abund_drt_fig.pdf", plot = RAC_abund_drt_fig, bg = "transparent", width =  6, height = 9, units = "in", dpi = 600)

########### Change during recovery; reference year 2017 ##############

sites<-unique(spcomp_long.final.rec$Site)
racchange.rec<-data.frame() #Make dataframe called racchange

plotinfo<-spcomp_long.final.rec %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Site,Plot,Trt,Block) %>% 
  unique ()


#need to put in console i=1 to test if it works
for(i in 1:length(sites)){
  sub<-spcomp_long.final.rec %>% 
    filter(Site==sites[i])
  out<-RAC_change(sub, time.var = "Year", species.var = "Spcode", 
                  abundance.var = "avg.cover", replicate.var = "Plot", 
                  reference.time = 2017)
  out$Site<-sites[i]  
  racchange.rec<-rbind(racchange.rec,out)
}

racchangetrt.rec<- merge(racchange, plotinfo, fix.by=c("Site","Plot"))

#Create a new column showing the year comparisons
racchangetrt.rec$Year.Year2<-str_c(racchangetrt.rec$Year,'.',racchangetrt.rec$Year2)

racchangetrt.long.rec<-racchangetrt %>% gather(Change.type,Value, richness_change:losses)

## abundance change 
sites<-unique(comppath.plotsum.rec$Site)
abunchng.rec<-data.frame() #Make dataframe called abunchng

#Making a dataset with treatment info
plotinfo<-comppath.plotsum.rec %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Site,Plot, Trt) %>% 
  unique()

#need to put in console i=1 to test if it works
for(i in 1:length(sites)){
  sub<-comppath.plotsum.rec %>% 
    filter(Site==sites[i])
  out.rec<-abundance_change(
    df=sub,
    time.var="Year",
    species.var="Path",
    abundance.var="pathsum",
    replicate.var = "Plot",
    reference.time = "2017"
  )
  out.rec$Site<-sites[i]  
  abunchng.rec<-rbind(abunchng.rec,out.rec)
}

#Merge dataframe with trt info
abunchngtrt.rec<-abunchng.rec %>% 
  left_join(plotinfo)

#Create a new column showing the year comparisons
abunchngtrt.rec$Year.Year2<-str_c(abunchngtrt.rec$Year,'-',abunchngtrt.rec$Year2)

#rename columns so they match with racchange

names(abunchngtrt.rec)[names(abunchngtrt.rec) == "Path"] <- "Change.type"
names(abunchngtrt.rec)[names(abunchngtrt.rec) == "change"] <- "Value"
racchangetrt.long.rec$Block <- NULL

str(racchangetrt.long)

#merge dfs 
RAC_abund.rec <- rbind(racchangetrt.long.rec, abunchngtrt.rec)
RAC_abund.rec$Change.type<-as.factor(RAC_abund.rec$Change.type)
str(RAC_abund.rec)

#calculate aves 

RAC_abund_ave.rec<-RAC_abund.rec%>%
  dplyr::group_by(Site,Trt,Change.type)%>%
  dplyr::summarize(mean_change = mean(Value,na.rm=T), n = n(),sd = sd(Value,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_change - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_change + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

RAC_abund_ave.rec$Trt <- factor(RAC_abund_ave.rec$Trt, levels=c('con', 'chr', 'int'))
RAC_abund_ave.rec$Change.type <- factor(RAC_abund_ave.rec$Change.type, levels=c('richness_change','evenness_change', 'rank_change',"gains","losses","C3","C4"), 
                                    labels = c("Richness Chg.","Evenness Chg.","Reordering","Sp. Gains","Sp. Losses","C3 Abund. Chg.","C4 Abund. Chg."))

#figure
RAC_abund_drt_fig.rec<-ggplot(RAC_abund_ave.rec, aes(x=Site, y=mean_change,color=Trt))+
  ylab("Average change during recovery")+
  xlab("Site")+
  facet_wrap(~Change.type, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3.5,position=position_dodge(.60))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.60), width=0, size=1.25, show.legend = TRUE)+
  scale_x_discrete(labels=c('SGS', 'HPG', 'HYS','KNZ'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                        axis.text.y = element_text(size=12, color = "black"),
                        strip.text = element_text(size=14),
                        axis.title.x = element_text(size=14),
                        axis.title.y = element_text(size=14),
                        legend.position = "top",legend.text = element_text(size=11))
RAC_abund_drt_fig.rec
ggsave(filename = "RAC_abund_drt_fig.rec.pdf", plot = RAC_abund_drt_fig.rec, bg = "transparent", width =  6, height = 9, units = "in", dpi = 600)



