## CODYN ANALYSES ##

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

#DOWNLOADING PAIRWISE ADONIS
#install.packages("devtools")
library(devtools)

#remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

spcomp_long.final<-read.csv("data/spdata.plot_names.csv")

#Make the following items into a factor:
spcomp_long.final<-spcomp_long.final %>% 
  dplyr::mutate_at(c('Site','Spcode','Year','Plot','Block','Trt'),as.factor)

spcomp_long.final$Site <- factor(spcomp_long.final$Site, levels=c('SGS', 'CHY', 'HYS','KNZ'))
spcomp_long.final$Trt <- factor(spcomp_long.final$Trt, levels=c('con', 'chr', 'int'))


spcomp_long.final<- spcomp_long.final %>%
  dplyr::filter(Year %in% c(2013,2014,2015,2016,2017,2018,2019,2020,2021)) 

#RAC CHANGE 
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
  out<-RAC_change(sub, time.var = "Year", species.var = "Spcode", 
                  abundance.var = "avg.cover", replicate.var = "Plot", reference.time = 2013)
  out$Site<-sites[i]  
  racchange<-rbind(racchange,out)
}

racchangetrt<- merge(racchange, plotinfo, fix.by=c("Site","Plot"))

#Switch to long format to calculate mean and CI
racchangetrt.long<-racchangetrt %>% gather(Change.type,Value, richness_change:losses)


#Create a new column showing the year comparisons
racchangetrt.long$Year.Year2<-str_c(racchangetrt$Year,'.',racchangetrt$Year2)

#Take average and CI
racave<-racchangetrt.long%>%
  dplyr::group_by(Site, Year.Year2,Trt,Change.type)%>%
  dplyr::summarize(mean_change = mean(Value,na.rm=T), n = n(),sd = sd(Value,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_change - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_change + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()
racave$Trt <- factor(racave$Trt, levels=c('con', 'chr', 'int'))


#racave$Site<-as.character(racave$Site)

drought_racave <- racave %>% dplyr::filter(Year.Year2 %in% c("2013.2014","2013.2015","2013.2016","2013.2017"))

recovery_racave <- racave %>% dplyr::filter(Year.Year2 %in% c("2013.2018","2013.2019","2013.2020","2013.2021"))


##### wrap by site, sep by metric #######

#richness change drought 
Rich_RC_drought<-ggplot(subset(drought_racave,Change.type=="richness_change"),aes(x=Year.Year2 , y=mean_change,color=Trt))+
  ylab("Change from pre-treatment")+
  xlab("Years of Treatment")+
  ggtitle("richness change")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
Rich_RC_drought

ggsave(filename = "Rich_RC_drought.jpeg", plot = Rich_RC_drought, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)

#evenness change drought 
even_RC_drought<-ggplot(subset(drought_racave,Change.type=="evenness_change"),aes(x=Year.Year2 , y=mean_change,color=Trt))+
  ylab("Change from pre-treatment")+
  xlab("Years of Treatment")+
  ggtitle("evenness change")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
even_RC_drought

ggsave(filename = "even_RC_drought.jpeg", plot = even_RC_drought, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)

#rank change drought 
rank_RC_drought<-ggplot(subset(drought_racave,Change.type=="rank_change"),aes(x=Year.Year2 , y=mean_change,color=Trt))+
  ylab("Change from pre-treatment")+
  xlab("Years of Treatment")+
  ggtitle("rank change")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
rank_RC_drought

ggsave(filename = "rank_RC_drought.jpeg", plot = rank_RC_drought, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)


#gains change drought 
gains_RC_drought<-ggplot(subset(drought_racave,Change.type=="gains"),aes(x=Year.Year2 , y=mean_change,color=Trt))+
  ylab("Change from pre-treatment")+
  xlab("Years of Treatment")+
  ggtitle("gains")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
gains_RC_drought

ggsave(filename = "gains_RC_drought.jpeg", plot = gains_RC_drought, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)

#losses change drought 
losses_RC_drought<-ggplot(subset(drought_racave,Change.type=="losses"),aes(x=Year.Year2 , y=mean_change,color=Trt))+
  ylab("Change from pre-treatment")+
  xlab("Years of Treatment")+
  ggtitle("losses")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
losses_RC_drought

ggsave(filename = "losses_RC_drought.jpeg", plot = losses_RC_drought, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)

## figure with average change ##

racchangetrt.long_drought <- racchangetrt.long %>% dplyr::filter(Year.Year2 %in% c("2013.2014","2013.2015","2013.2016","2013.2017"))

racave_ave<-racchangetrt.long_drought%>%
  dplyr::group_by(Site,Trt,Change.type)%>%
  dplyr::summarize(mean_change = mean(Value,na.rm=T), n = n(),sd = sd(Value,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_change - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_change + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()
racave_ave$Trt <- factor(racave_ave$Trt, levels=c('con', 'chr', 'int'))


Drought_RACave_fig<-ggplot(racave_ave, aes(x=Site, y=mean_change,color=Trt))+
  ylab("Mean change")+
  xlab("Site")+
  facet_wrap(~Change.type,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.45))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.45), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('SGS', 'HPG', 'HYS','KNZ'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
Drought_RACave_fig
ggsave(filename = "Drought_RACave_fig.jpeg", plot = Drought_RACave_fig, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)


########### RAC change recovery compared to 2017 ##############
sites<-unique(spcomp_long.final$Site)
racchange17<-data.frame() #Make dataframe called racchange

plotinfo<-spcomp_long.final %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Site,Plot,Trt,Block) %>% 
  unique ()


#need to put in console i=1 to test if it works
for(i in 1:length(sites)){
  sub<-spcomp_long.final %>% 
    filter(Site==sites[i])
  out<-RAC_change(sub, time.var = "Year", species.var = "Spcode", 
                  abundance.var = "avg.cover", replicate.var = "Plot", 
                  reference.time = 2017)
  out$Site<-sites[i]  
  racchange17<-rbind(racchange17,out)
}

racchangetrt17<- merge(racchange17, plotinfo, fix.by=c("Site","Plot"))

#Switch to long format to calculate mean and CI
racchangetrt17.long<-racchangetrt17 %>% gather(Change.type,Value, richness_change:losses)


#Create a new column showing the year comparisons
racchangetrt17.long$Year.Year2<-str_c(racchangetrt17$Year,'.',racchangetrt17$Year2)

#Take average and CI
racave17<-racchangetrt17.long%>%
  dplyr::group_by(Site, Year.Year2,Trt,Change.type)%>%
  dplyr::summarize(mean_change = mean(Value,na.rm=T), n = n(),sd = sd(Value,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_change - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_change + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()
racave17$Trt <- factor(racave$Trt, levels=c('con', 'chr', 'int'))

recovery_racave17 <- racave17 %>% dplyr::filter(Year.Year2 %in% c("2017.2018","2017.2019","2017.2020","2017.2021"))

#richness change recovery 2017
Rich_RC17_recovery<-ggplot(subset(recovery_racave17,Change.type=="richness_change"),aes(x=Year.Year2 , y=mean_change,color=Trt))+
  ylab("Change from end of treatment")+
  xlab("Years of Recovery")+
  ggtitle("richness change")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
Rich_RC17_recovery

ggsave(filename = "Rich_RC17_recovery.jpeg", plot = Rich_RC17_recovery, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)

#evenness change recovery 
even_RC17_recovery<-ggplot(subset(recovery_racave17,Change.type=="evenness_change"),aes(x=Year.Year2 , y=mean_change,color=Trt))+
  ylab("Change from end of treatment")+
  xlab("Years of Recovery")+
  ggtitle("evenness change")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
even_RC17_recovery

ggsave(filename = "even_RC17_recovery.jpeg", plot = even_RC17_recovery, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)

#rank change recovery 
rank_RC17_recovery<-ggplot(subset(recovery_racave17,Change.type=="rank_change"),aes(x=Year.Year2 , y=mean_change,color=Trt))+
  ylab("Change from end of treatment")+
  xlab("Years of Recovery")+
  ggtitle("rank change")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
rank_RC17_recovery

ggsave(filename = "rank_RC17_recovery.jpeg", plot = rank_RC17_recovery, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)

#gains change recovery 
gains_RC17_recovery<-ggplot(subset(recovery_racave17,Change.type=="gains"),aes(x=Year.Year2 , y=mean_change,color=Trt))+
  ylab("Change from end of treatment")+
  xlab("Years of Recovery")+
  ggtitle("gains")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
gains_RC17_recovery

ggsave(filename = "gains_RC17_recovery.jpeg", plot = gains_RC17_recovery, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)

#losses change recovery 
losses_RC17_recovery<-ggplot(subset(recovery_racave17,Change.type=="losses"),aes(x=Year.Year2 , y=mean_change,color=Trt))+
  ylab("Change from end of treatment")+
  xlab("Years of Recovery")+
  ggtitle("losses")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
losses_RC17_recovery

ggsave(filename = "losses_RC17_recovery.jpeg", plot = losses_RC17_recovery, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)


## figure with average change ##

racchangetrt.long_recovery <- racchangetrt17.long %>% dplyr::filter(Year.Year2 %in% c("2017.2018","2017.2019","2017.2020","2017.2021"))

racave_ave_recov<-racchangetrt.long_recovery%>%
  dplyr::group_by(Site,Trt,Change.type)%>%
  dplyr::summarize(mean_change = mean(Value,na.rm=T), n = n(),sd = sd(Value,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_change - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_change + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

racave_ave_recov$Trt <- factor(racave_ave$Trt, levels=c('con', 'chr', 'int'))


#losses change recovery 
Recovery_RACave_fig<-ggplot(racave_ave_recov, aes(x=Site, y=mean_change,color=Trt))+
  ylab("Mean change")+
  xlab("Site")+
  facet_wrap(~Change.type,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('SGS', 'HPG', 'HYS','KNZ'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
Recovery_RACave_fig

ggsave(filename = "Recovery_RACave_fig.jpeg", plot = Drought_RACave_fig, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)


## figure with average change ##

racchangetrt.long_drought <- racchangetrt.long %>% dplyr::filter(Year.Year2 %in% c("2013.2014","2013.2015","2013.2016","2013.2017"))

racave_ave<-racchangetrt.long_drought%>%
  dplyr::group_by(Site,Trt,Change.type)%>%
  dplyr::summarize(mean_change = mean(Value,na.rm=T), n = n(),sd = sd(Value,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_change - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_change + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()
racave_ave$Trt <- factor(racave_ave$Trt, levels=c('con', 'chr', 'int'))


#losses change drought 
Drought_RACave_fig<-ggplot(racave_ave, aes(x=Site, y=mean_change,color=Trt))+
  ylab("Mean change")+
  xlab("Site")+
  facet_wrap(~Change.type,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('SGS', 'HPG', 'HYS','KNZ'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
Drought_RACave_fig
ggsave(filename = "Drought_RACave_fig.jpeg", plot = Drought_RACave_fig, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)




######################################################################

#using diff dataset for abundance. Just looking at C3/C4 
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

#ABUNDANCE CHANGE====
#Calculates the abundance change for species in a replicate between two time points. Changes are
#on abundance values provided, if relative data is used, then changes in relative abundance will be
#calculated.

#"For loop" to calculate abundance change
#Extract unique values for Site (Bison, Cattle)

sites<-unique(comppath.plotsum$Site)
abunchng<-data.frame() #Make dataframe called abunchng

#Making a dataset with treatment info
plotinfo<-comppath.plotsum %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Site,Plot, Trt) %>% 
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

#Average and CI by time point comparison
#can remove trt so that n=>1 for all spp, but then it's not super useful
abunchng.ave<-abunchngtrt%>%
  dplyr::group_by(Site, Year.Year2,Trt,Path)%>%
  dplyr::summarize(mean_chng = mean(change,na.rm=T), n = n(),sd = sd(change,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_chng - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_chng + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

drought_ACave <- abunchng.ave %>% dplyr::filter(Year.Year2 %in% c("2013-2014","2013-2015","2013-2016","2013-2017"))


### wrap by Site, sep by path 
#C3
C3_AC_Drought<-ggplot(subset(drought_ACave,Path=="C3"),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pre-treatment")+
  xlab("Years of Treatment")+
  ggtitle("C3")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
C3_AC_Drought

ggsave(filename = "C3_AC_Drought.jpeg", plot = C3_AC_Drought, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)


#C4
C4_AC_Drought<-ggplot(subset(drought_ACave,Path=="C4"),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pre-treatment")+
  xlab("Years of Treatment")+
  ggtitle("C4")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
C4_AC_Drought

ggsave(filename = "C4_AC_Drought.jpeg", plot = C4_AC_Drought, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)


## C3 and C4 aves not by year 

abundchngtrt_drt <- abunchngtrt %>% dplyr::filter(Year.Year2 %in% c("2013-2014","2013-2015","2013-2016","2013-2017"))

abundchg_main_ave<-abundchngtrt_drt%>%
  dplyr::group_by(Site,Trt,Path)%>%
  dplyr::summarize(mean_change = mean(change,na.rm=T), n = n(),sd = sd(change,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_change - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_change + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

drt_abundave_fig<-ggplot(abundchg_main_ave, aes(x=Site, y=mean_change,color=Trt))+
  ylab("Mean change")+
  xlab("Site")+
  facet_wrap(~Path,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('SGS', 'HPG', 'HYS','KNZ'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
drt_abundave_fig


#### ABUND and RAC change drought datasets combined

#change datasets to wide format 



#### recovery change compared to 2017 ####

abunchng17<-data.frame()
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
abunchngtrt17 <- abunchngtrt17 %>% dplyr::filter(Year.Year2 %in% c("2017-2018","2017-2019","2017-2020","2017-2021"))

#Average and CI by time point comparison
#can remove trt so that n=>1 for all spp, but then it's not super useful
abunchng.ave17<-abunchngtrt17%>%
  dplyr::group_by(Site, Year.Year2,Trt,Path)%>%
  dplyr::summarize(mean_chng = mean(change,na.rm=T), n = n(),sd = sd(change,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_chng - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_chng + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

C3_AC_recov17<-ggplot(subset(recovery_ACave17,Path=="C3"),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from end of treatment")+
  xlab("Years of Recovery")+
  ggtitle("C3")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
C3_AC_recov17

ggsave(filename = "C3_AC_recov17.jpeg", plot = C3_AC_recov17, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)


C4_AC_recov17<-ggplot(subset(recovery_ACave17,Path=="C4"),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from end of treatment")+
  xlab("Years of Recovery")+
  ggtitle("C4")+
  facet_wrap(~Site,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
C4_AC_recov17

ggsave(filename = "C4_AC_recov17.jpeg", plot = C4_AC_recov17, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)

## C3 and C4 aves not by year 

abundchgave17_ave<-abunchngtrt17%>%
  dplyr::group_by(Site,Trt,Path)%>%
  dplyr::summarize(mean_change = mean(change,na.rm=T), n = n(),sd = sd(change,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_change - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_change + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()


#abund change recov compared to 2017
recov_abundave_fig<-ggplot(abundchgave17_ave, aes(x=Site, y=mean_change,color=Trt))+
  ylab("Mean change")+
  xlab("Site")+
  facet_wrap(~Path,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('SGS', 'HPG', 'HYS','KNZ'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
recov_abundave_fig
ggsave(filename = "recov_abundave_fig.jpeg", plot = recov_abundave_fig, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)














