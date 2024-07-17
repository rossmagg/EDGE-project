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

# C3/C4 Abundance change during drought, compared to 2013 
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

drought_ACave$Trt <- factor(drought_ACave$Trt, levels=c('con', 'chr', 'int'))

#fig grid 

all_path_yr.drt<-ggplot(drought_ACave,aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pre-treatment")+
  xlab("Years of Treatment")+
  geom_hline(yintercept = 0, color="grey")+
  facet_grid(Path~Site, scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=14),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "top",
                   legend.text = element_text(size=11),
                   legend.title = element_text(size=11),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
all_path_yr.drt

ggsave(filename = "all_path_yr.drt.pdf", plot = all_path_yr.drt, bg = "transparent", width =  11, height = 8, units = "in", dpi = 600)


## C3/C4 Abundance change during recovery, compared to 2017

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
abunchngtrt17$Trt <- factor(abunchngtrt17$Trt, levels=c('con', 'chr', 'int'))

#Average and CI by time point comparison

abunchng.ave17<-abunchngtrt17%>%
  dplyr::group_by(Site, Year.Year2,Trt,Path)%>%
  dplyr::summarize(mean_chng = mean(change,na.rm=T), n = n(),sd = sd(change,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_chng - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_chng + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()


rec_ACave17 <- abunchng.ave17 %>% dplyr::filter(Year.Year2 %in% c("2017-2018","2017-2019","2017-2020","2017-2021"))

rec_ACave17$Trt <- factor(rec_ACave17$Trt, levels=c('con', 'chr', 'int'))

#fig grid 

all_path_yr.rec<-ggplot(rec_ACave17,aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from end of treatment")+
  xlab("Years of Recovery")+
  geom_hline(yintercept = 0, color="grey")+
  facet_grid(Path~Site, scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=14),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "top",
                   legend.text = element_text(size=11),
                   legend.title = element_text(size=11),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
all_path_yr.rec

ggsave(filename = "all_path_yr.rec.pdf", plot = all_path_yr.rec, bg = "transparent", width =  11, height = 8, units = "in", dpi = 600)


### code can be deleted below ###

#figs by site so scales are free for each path by site 
SGS_path_drt <- ggplot(subset(drought_ACave,Site=="SGS"),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pre-treatment")+
  xlab("Years of Treatment")+
  geom_hline(yintercept = 0, color="grey")+
  facet_wrap(~Path, scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=12),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none")
SGS_path_drt 

CHY_path_drt <- ggplot(subset(drought_ACave,Site=="CHY"),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pre-treatment")+
  xlab("Years of Treatment")+
  geom_hline(yintercept = 0, color="grey")+
  facet_wrap(~Path, scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=12),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none")
CHY_path_drt 

HYS_path_drt <- ggplot(subset(drought_ACave,Site=="HYS"),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pre-treatment")+
  xlab("Years of Treatment")+
  geom_hline(yintercept = 0, color="grey")+
  facet_wrap(~Path, scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=12),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none")
HYS_path_drt 


KNZ_path_drt <- ggplot(subset(drought_ACave,Site=="KNZ"),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change from pre-treatment")+
  xlab("Years of Treatment")+
  geom_hline(yintercept = 0, color="grey")+
  facet_wrap(~Path, scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=12),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none")
KNZ_path_drt 

all_path_drt_yr <- ggarrange(SGS_path_drt,CHY_path_drt,HYS_path_drt,KNZ_path_drt)

#figs by site 

SGS_path_rec <- ggplot(subset(rec_ACave17,Site=="SGS"),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change since end of treatment")+
  xlab("Years of recovery")+
  geom_hline(yintercept = 0, color="grey")+
  facet_wrap(~Path, scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=12),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none")
SGS_path_rec

CHY_path_rec <- ggplot(subset(rec_ACave17,Site=="CHY"),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change since end of treatment")+
  xlab("Years of recovery")+
  geom_hline(yintercept = 0, color="grey")+
  facet_wrap(~Path, scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=12),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none")
CHY_path_rec

HYS_path_rec <- ggplot(subset(rec_ACave17,Site=="HYS"),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change since end of treatment")+
  xlab("Years of recovery")+
  geom_hline(yintercept = 0, color="grey")+
  facet_wrap(~Path, scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=12),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none")
HYS_path_rec

KNZ_path_rec <- ggplot(subset(rec_ACave17,Site=="KNZ"),aes(x=Year.Year2 , y=mean_chng,color=Trt))+
  ylab("Change since end of treatment")+
  xlab("Years of recovery")+
  geom_hline(yintercept = 0, color="grey")+
  facet_wrap(~Path, scales = "free")+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  scale_x_discrete(labels=c('1', '2', '3','4'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=12),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none")
KNZ_path_rec

all_path_rec_yr <- ggarrange(SGS_path_rec,CHY_path_rec,HYS_path_rec,KNZ_path_rec)

ggsave(filename = "all_path_rec_yr.pdf", plot = all_path_rec_yr, bg = "transparent", width =  11, height = 8, units = "in", dpi = 600)

