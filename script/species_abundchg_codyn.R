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

####### species abudances by site ########
#### abund change by spp, not path ##### 
comp.species.simper <- comp.species %>% dplyr::filter(Spcode %in% c("BOGR","BRJA","BRTE","CAEL","HECO","KOMA","PASM","VUOC",
                                                                    "ELEL","BOCU","ANGE","SCSC","SPAS","SONU","POPR",
                                                                    "CAEL"))
comp.species.simper$Trt <- factor(comp.species.simper$Trt, levels=c('con', 'chr', 'int'))
comp.species.simper$Site <- factor(comp.species.simper$Site, levels=c('SGS', 'CHY', 'HYS','KNZ'))


sites<-unique(comp.species.simper$Site)
abunchng_allspp<-data.frame() #Make dataframe called abunchng

#Making a dataset with treatment info
plotinfo<-comp.species.simper %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Site,Plot, Trt) %>% 
  unique()

#need to put in console i=1 to test if it works
for(i in 1:length(sites)){
  sub<-comp.species.simper%>% 
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
  abunchng_allspp<-rbind(abunchng_allspp,out)
}

#Merge dataframe with trt info
abunchng_allspp.trt<-abunchng_allspp %>% 
  left_join(plotinfo)

#Create a new column showing the year comparisons
abunchng_allspp.trt$Year.Year2<-str_c(abunchng_allspp.trt$Year,'-',abunchng_allspp.trt$Year2)

#Average and CI by time point comparison
#can remove trt so that n=>1 for all spp, but then it's not super useful
abunchng_allspp.ave<-abunchng_allspp.trt%>%
  dplyr::group_by(Site, Year.Year2,Trt,Spcode)%>%
  dplyr::summarize(mean_chng = mean(change,na.rm=T), n = n(),sd = sd(change,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_chng - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_chng + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

drought_ACave.spp <- abunchng_allspp.ave %>% dplyr::filter(Year.Year2 %in% c("2013-2014","2013-2015","2013-2016","2013-2017"))

recovery_ACave.spp <- abunchng_allspp.ave %>% dplyr::filter(Year.Year2 %in% c("2013-2018","2013-2019","2013-2020","2013-2021"))

drought_ACave.spp.SGS <- drought_ACave.spp %>% dplyr::filter(Site=='SGS')
sppchange.SGS.drt<-ggplot(subset(drought_ACave.spp.SGS, Spcode %in% c('BOGR',"VUOC","BRTE","ELEL")),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pretreatment")+
  xlab("Years of Treatment")+
  facet_wrap(~Spcode)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), legend.position = "top",legend.title = element_blank(),
                   strip.text = element_text(size = 14))
sppchange.SGS.drt
ggsave(filename = "sppchange.SGS.drt.jpeg", plot = sppchange.SGS.drt, bg = "transparent", width =  9, height = 6, units = "in", dpi = 600)



recovery_ACave.spp.SGS <- recovery_ACave.spp %>% dplyr::filter(Site=='SGS')
sppchange.SGS.rec<-ggplot(subset(recovery_ACave.spp.SGS, Spcode %in% c('BOGR',"VUOC","BRTE","ELEL")),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pretreatment")+
  xlab("Years of Recovery")+
  facet_wrap(~Spcode)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), legend.position = "top",legend.title = element_blank(),
                   strip.text = element_text(size = 14))
sppchange.SGS.rec

ggsave(filename = "sppchange.SGS.rec.jpeg", plot = sppchange.SGS.rec, bg = "transparent", width =  9, height = 6, units = "in", dpi = 600)


sppchange.CHY.drt<-ggplot(subset(drought_ACave.spp, Site=='CHY'),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  facet_wrap(~Spcode)+
  #scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
sppchange.CHY.drt

sppchange.HYS.drt<-ggplot(subset(drought_ACave.spp, Site=='HYS'),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  facet_wrap(~Spcode)+
  #scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"), legend.position = "top",legend.title = element_blank())
sppchange.HYS.drt
