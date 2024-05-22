### Species abundance changes by period ###
#sp selected based on SIMPER analysis 

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
comp.aveplot_trt$Trt <- factor(comp.aveplot_trt$Trt, levels=c('con', 'chr', 'int'))
comp.aveplot_trt$Site <- factor(comp.aveplot_trt$Site, levels=c('SGS', 'CHY', 'HYS','KNZ'))


#drought and recov df
comp.aveplot.drt<- comp.aveplot_trt %>%
  dplyr::filter(Year %in% c(2013,2014,2015,2016,2017)) 

comp.aveplot.rec<- comp.aveplot_trt %>%
  dplyr::filter(Year %in% c(2017,2018,2019,2020,2021)) 

####### species abudances by site ########
#### abund change by spp, not path ##### 

sites<-unique(comp.aveplot.drt$Site)
abunchng_allspp.drt<-data.frame() #Make dataframe called abunchng

#Making a dataset with treatment info
plotinfo<-comp.aveplot.drt %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Site,Plot, Trt) %>% 
  unique()

#need to put in console i=1 to test if it works
for(i in 1:length(sites)){
  sub<-comp.aveplot.drt%>% 
    filter(Site==sites[i])
  out<-abundance_change(
    df=sub,
    time.var="Year",
    species.var="Spcode",
    abundance.var="avg.cover",
    replicate.var = "Plot",
    reference.time = "2013"
  )
  out$Site<-sites[i]  
  abunchng_allspp.drt<-rbind(abunchng_allspp.drt,out)
}

#Merge dataframe with trt info
abunchng_allspp.drt_trt<-abunchng_allspp.drt %>% 
  left_join(plotinfo)

#Create a new column showing the year comparisons
abunchng_allspp.drt_trt$Year.Year2<-str_c(abunchng_allspp.drt_trt$Year,'-',abunchng_allspp.drt_trt$Year2)

#Average and CI by time point comparison
#can remove trt so that n=>1 for all spp, but then it's not super useful
abunchng_allspp.drt.ave<-abunchng_allspp.drt_trt%>%
  dplyr::group_by(Site, Year.Year2,Trt,Spcode)%>%
  dplyr::summarize(mean_chng = mean(change,na.rm=T), n = n(),sd = sd(change,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_chng - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_chng + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()


#drought fig 

#sep by site and spp of interest based on SIMPER analysis and rename so unique for each site 
drt.SGS <- abunchng_allspp.drt.ave %>% dplyr::filter(Site=='SGS', Spcode %in% c("46","317","55"))
drt.SGS$Spcode <- factor(drt.SGS$Spcode, levels = c("46","317","55"),
                              labels = c("Bouteloua gracilis","Vulpia octoflora","Carex eleocharis"))

drt.CHY <- abunchng_allspp.drt.ave %>% dplyr::filter(Site=='CHY',Spcode %in% c("46","317","197"))
drt.CHY$Spcode <- factor(drt.CHY$Spcode, levels = c("46","317","197"),
                         labels = c("Bouteloua gracilis","Vulpia octoflora","Pascopyrum smithii"))


drt.HYS <- abunchng_allspp.drt.ave %>% dplyr::filter(Site=='HYS',Spcode %in% c("44","48","197","235"))
drt.HYS$Spcode <- factor(drt.HYS$Spcode, levels = c("44","48","197","235"),
                         labels = c("Bouteloua curtipendula","Bromus japonicus",
                                    "Pascopyrum smithii","Sporobolus asper"))


drt.KNZ <- abunchng_allspp.drt.ave %>% dplyr::filter(Site=='KNZ', Spcode %in% c("7", "221","234"))
drt.KNZ$Spcode <- factor(drt.KNZ$Spcode, levels = c("7","221","234"),
                        labels = c("Andropogon gerardii","Schizachryium scoparium",
                                   "Sorghastrum nutans"))

#remerge dfs
#drt.allsites <- rbind(drt.SGS,drt.CHY,drt.HYS,drt.KNZ)

SGS_drt.fig<-ggplot(drt.SGS, aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pretreatment")+
  xlab("Years of Treatment")+
  facet_wrap(~Spcode, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.40))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.40), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+ 
  theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color = "black"),
        strip.text = element_text(size=12, face="italic"),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "none",legend.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
SGS_drt.fig

CHY_drt.fig<-ggplot(drt.CHY, aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pretreatment")+
  xlab("Years of Treatment")+
  facet_wrap(~Spcode, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.40))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.40), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+ 
  theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color = "black"),
        strip.text = element_text(size=12, face="italic"),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "none",legend.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
CHY_drt.fig

HYS_drt.fig<-ggplot(drt.HYS, aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pretreatment")+
  xlab("Years of Treatment")+
  facet_wrap(~Spcode, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.40))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.40), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+ 
  theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color = "black"),
        strip.text = element_text(size=12, face="italic"),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "none",legend.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
HYS_drt.fig

KNZ_drt.fig<-ggplot(drt.KNZ, aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pretreatment")+
  xlab("Years of Treatment")+
  facet_wrap(~Spcode, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.40))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.40), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+ 
  theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color = "black"),
        strip.text = element_text(size=12, face="italic"),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "none",legend.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
KNZ_drt.fig

all.drt<-ggarrange(SGS_drt.fig,CHY_drt.fig,HYS_drt.fig,KNZ_drt.fig, ncol=2,labels = c('a)', 'b)','c)','d)'))
all.drt

ggsave(filename = "all.drt.pdf", plot = all.drt, bg = "transparent", width =  11, height = 10, units = "in", dpi = 600)


### recovery ###

sites<-unique(comp.aveplot.rec$Site)
abunchng_allspp.rec<-data.frame() #Make dataframe called abunchng

#Making a dataset with treatment info
plotinfo<-comp.aveplot.rec %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Site,Plot, Trt) %>% 
  unique()

#need to put in console i=1 to test if it works
for(i in 1:length(sites)){
  sub<-comp.aveplot.rec%>% 
    filter(Site==sites[i])
  out<-abundance_change(
    df=sub,
    time.var="Year",
    species.var="Spcode",
    abundance.var="avg.cover",
    replicate.var = "Plot",
    reference.time = "2017"
  )
  out$Site<-sites[i]  
  abunchng_allspp.rec<-rbind(abunchng_allspp.rec,out)
}

#Merge dataframe with trt info
abunchng_allspp.rec_trt<-abunchng_allspp.rec %>% 
  left_join(plotinfo)

#Create a new column showing the year comparisons
abunchng_allspp.rec_trt$Year.Year2<-str_c(abunchng_allspp.rec_trt$Year,'-',abunchng_allspp.rec_trt$Year2)

#Average and CI by time point comparison
#can remove trt so that n=>1 for all spp, but then it's not super useful
abunchng_allspp.rec.ave<-abunchng_allspp.rec_trt%>%
  dplyr::group_by(Site, Year.Year2,Trt,Spcode)%>%
  dplyr::summarize(mean_chng = mean(change,na.rm=T), n = n(),sd = sd(change,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_chng - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_chng + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()


rec.SGS <- abunchng_allspp.rec.ave %>% dplyr::filter(Site=='SGS', Spcode %in% c("49","317","95"))
rec.SGS$Spcode <- factor(rec.SGS$Spcode, levels = c("49","317","95"),
                         labels = c("BRJA","Vulpia octoflora", "ELEL"))

rec.CHY <- abunchng_allspp.rec.ave %>% dplyr::filter(Site=='CHY',Spcode %in% c("137","130","55","197"))
rec.CHY$Spcode <- factor(rec.CHY$Spcode, levels = c("137","130","55","197"),
                         labels = c("KOMA","HECO","CAEL","PASM"))


rec.HYS <- abunchng_allspp.rec.ave %>% dplyr::filter(Site=='HYS',Spcode %in% c("235","221","5","44"))
rec.HYS$Spcode <- factor(rec.HYS$Spcode, levels = c("235","221","5","44"),
                         labels = c("SPAS","SCSC",
                                    "AMPS","BOCU"))


rec.KNZ <- abunchng_allspp.rec.ave %>% dplyr::filter(Site=='KNZ', Spcode %in% c("7", "221","234","5"))
rec.KNZ$Spcode <- factor(rec.KNZ$Spcode, levels = c("7","221","234","5"),
                         labels = c("Andropogon gerardii","Schizachryium scoparium",
                                    "Sporobolus asper","AMPS"))

#remerge dfs
#rec.allsites <- rbind(rec.SGS,rec.CHY,rec.HYS,rec.KNZ)

SGS_rec.fig<-ggplot(rec.SGS, aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pretreatment")+
  xlab("Years of Treatment")+
  facet_wrap(~Spcode, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.40))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.40), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+ 
  theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color = "black"),
        strip.text = element_text(size=12, face="italic"),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "none",legend.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
SGS_rec.fig

CHY_rec.fig<-ggplot(rec.CHY, aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pretreatment")+
  xlab("Years of Treatment")+
  facet_wrap(~Spcode, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.40))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.40), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+ 
  theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color = "black"),
        strip.text = element_text(size=12, face="italic"),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "none",legend.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
CHY_rec.fig

HYS_rec.fig<-ggplot(rec.HYS, aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pretreatment")+
  xlab("Years of Treatment")+
  facet_wrap(~Spcode, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.40))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.40), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+ 
  theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color = "black"),
        strip.text = element_text(size=12, face="italic"),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "none",legend.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
HYS_rec.fig

KNZ_rec.fig<-ggplot(rec.KNZ, aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pretreatment")+
  xlab("Years of Treatment")+
  facet_wrap(~Spcode, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.40))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.40), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+ 
  theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color = "black"),
        strip.text = element_text(size=12, face="italic"),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "none",legend.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
KNZ_rec.fig

all.rec<-ggarrange(SGS_rec.fig,CHY_rec.fig,HYS_rec.fig,KNZ_rec.fig, ncol=2,labels = c('a)', 'b)','c)','d)'))
all.rec

ggsave(filename = "all.rec.pdf", plot = all.rec, bg = "transparent", width =  11, height = 10, units = "in", dpi = 600)