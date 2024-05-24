
### RAC change by period ###

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

RAC_drt_ave<-racchangetrt.long%>%
  dplyr::group_by(Site,Trt,Change.type)%>%
  dplyr::summarize(mean_change = mean(Value,na.rm=T), n = n(),sd = sd(Value,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_change - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_change + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

RAC_drt_ave$Trt <- factor(RAC_drt_ave$Trt, levels=c('con', 'chr', 'int'))
RAC_drt_ave$Change.type <- factor(RAC_drt_ave$Change.type, levels=c('richness_change','evenness_change', 'rank_change',"gains","losses"), 
                                    labels = c("Richness","Evenness","Reordering","Sp. Gains","Sp. Losses"))

#figure
RAC_drt_fig<-ggplot(RAC_drt_ave, aes(x=Site, y=mean_change,color=Trt))+
  ylab("Average change during drought")+
  xlab("Site")+
  geom_hline(yintercept = 0, color="grey")+
  facet_wrap(~Change.type, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.40))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.40), width=0, linewidth=1.25, show.legend = TRUE)+
  scale_x_discrete(labels=c('SGS', 'HPG', 'HYS','KNZ'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                        axis.text.y = element_text(size=12, color = "black"),
                        strip.text = element_text(size=12),
                        axis.title.x = element_text(size=14),
                        axis.title.y = element_text(size=14),
                        legend.position = "top",
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())
RAC_drt_fig
ggsave(filename = "RAC_drt_fig.pdf", plot = RAC_drt_fig, bg = "transparent", width =  6, height = 9, units = "in", dpi = 600)

### RAC change during recovery, compared to 2017 ###

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
  out17<-RAC_change(sub, time.var = "Year", species.var = "Spcode", 
                  abundance.var = "avg.cover", replicate.var = "Plot", 
                  reference.time = 2017)
  out17$Site<-sites[i]  
  racchange.rec<-rbind(racchange.rec,out17)
}

racchangetrt.rec<- merge(racchange.rec, plotinfo, fix.by=c("Site","Plot"))

#Create a new column showing the year comparisons
racchangetrt.rec$Year.Year2<-str_c(racchangetrt.rec$Year,'.',racchangetrt.rec$Year2)

racchangetrt.long.rec<-racchangetrt.rec %>% gather(Change.type,Value, richness_change:losses)

RAC_rec_ave<-racchangetrt.long.rec%>%
  dplyr::group_by(Site,Trt,Change.type)%>%
  dplyr::summarize(mean_change = mean(Value,na.rm=T), n = n(),sd = sd(Value,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_change - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_change + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

RAC_rec_ave$Trt <- factor(RAC_rec_ave$Trt, levels=c('con', 'chr', 'int'))
RAC_rec_ave$Change.type <- factor(RAC_rec_ave$Change.type, levels=c('richness_change','evenness_change', 'rank_change',"gains","losses"), 
                                  labels = c("Richness","Evenness","Reordering","Sp. Gains","Sp. Losses"))

RAC_rec_fig<-ggplot(RAC_rec_ave, aes(x=Site, y=mean_change,color=Trt))+
  ylab("Average change during recovery")+
  xlab("Site")+
  geom_hline(yintercept = 0, color="grey")+
  facet_wrap(~Change.type, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.40))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.40), width=0, linewidth=1.25, show.legend = TRUE)+
  scale_x_discrete(labels=c('SGS', 'HPG', 'HYS','KNZ'))+
  theme_classic()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                        axis.text.y = element_text(size=12, color = "black"),
                        strip.text = element_text(size=12),
                        axis.title.x = element_text(size=14),
                        axis.title.y = element_text(size=14),
                        legend.position = "top",
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())
RAC_rec_fig
ggsave(filename = "RAC_rec_fig.pdf", plot = RAC_rec_fig, bg = "transparent", width =  6, height = 9, units = "in", dpi = 600)

