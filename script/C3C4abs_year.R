## C3/C4 avg cover over time ##
library(tidyverse)
library(ggarrange)


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

path_trt_ave<-comppath.plotsum%>%
  dplyr::group_by(Site,Year,Trt,Path)%>%
  dplyr::summarize(ave = mean(pathsum,na.rm=T), n = n(),sd = sd(pathsum,na.rm=T), se = sd/sqrt(n),
                   LowerCI = ave - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = ave + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

path_trt_ave$Trt <- factor(path_trt_ave$Trt, levels=c('con', 'chr', 'int'))
path_trt_ave$Site <- factor(path_trt_ave$Site, levels=c('SGS', 'CHY', 'HYS','KNZ'))

Drt_abs_ave <-path_trt_ave %>% dplyr::filter(Year %in% c("2014","2015","2016","2017"))
Rec_abs_ave <- path_trt_ave %>% dplyr::filter(Year %in% c("2018","2019","2020","2021"))

#fig grid drt 

Drt_path_abs_grid <- ggplot(Drt_abs_ave,aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year of Treatment")+
  facet_grid(Path~Site, scales = "free")+
  scale_x_discrete(labels=c('1', '2','3','4'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
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
Drt_path_abs_grid 
ggsave(filename = "Drt_path_abs_grid.pdf", plot = Drt_path_abs_grid, bg = "transparent", width =  8, height = 6, units = "in", dpi = 600)

#fig grid recovery 

Rec_path_abs_grid <- ggplot(Rec_abs_ave,aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year of Recovery")+
  facet_grid(Path~Site, scales = "free")+
  scale_x_discrete(labels=c('1', '2','3','4'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
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
Rec_path_abs_grid 
ggsave(filename = "Rec_path_abs_grid.pdf", plot = Rec_path_abs_grid, bg = "transparent", width =  8, height = 6, units = "in", dpi = 600)


#figs by site

#SGS
Drt_path_abs_SGS <- ggplot(subset(Drt_abs_ave, Site=="SGS"),aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year of Treatment")+
  facet_wrap(~Path, scales = "free")+
  scale_x_discrete(labels=c('1', '2','3','4'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=14),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
Drt_path_abs_SGS

Rec_path_abs_SGS <- ggplot(subset(Rec_abs_ave, Site=="SGS"),aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year of Recovery")+
  facet_wrap(~Path, scales = "free")+
  scale_x_discrete(labels=c('1', '2','3','4'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=14),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
Rec_path_abs_SGS

#CHY
Drt_path_abs_CHY <- ggplot(subset(Drt_abs_ave, Site=="CHY"),aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year of Treatment")+
  facet_wrap(~Path, scales = "free")+
  scale_x_discrete(labels=c('1', '2','3','4'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=14),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
Drt_path_abs_CHY


Rec_path_abs_CHY<- ggplot(subset(Rec_abs_ave, Site=="CHY"),aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year of Recovery")+
  facet_wrap(~Path, scales = "free")+
  scale_x_discrete(labels=c('1', '2','3','4'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=14),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
Rec_path_abs_CHY

#HYS
Drt_path_abs_HYS<- ggplot(subset(Drt_abs_ave, Site=="HYS"),aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year of Treatment")+
  facet_wrap(~Path, scales = "free")+
  scale_x_discrete(labels=c('1', '2','3','4'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=14),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
Drt_path_abs_HYS

Rec_path_abs_HYS<- ggplot(subset(Rec_abs_ave, Site=="HYS"),aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year of Recovery")+
  facet_wrap(~Path, scales = "free")+
  scale_x_discrete(labels=c('1', '2','3','4'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=14),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
Rec_path_abs_HYS

#KNZ
Drt_path_abs_KNZ<- ggplot(subset(Drt_abs_ave, Site=="KNZ"),aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year of Treatment")+
  facet_wrap(~Path, scales = "free")+
  scale_x_discrete(labels=c('1', '2','3','4'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=14),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
Drt_path_abs_KNZ

Rec_path_abs_KNZ<- ggplot(subset(Rec_abs_ave, Site=="KNZ"),aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year of Recovery")+
  facet_wrap(~Path, scales = "free")+
  scale_x_discrete(labels=c('1', '2','3','4'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=14),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
Rec_path_abs_KNZ


Drt_path_abs <- Drt_path_abs_SGS+Drt_path_abs_CHY+Drt_path_abs_HYS+Drt_path_abs_KNZ 
Drt_path_abs
ggsave(filename = "Drt_path_abs.pdf", plot = Drt_path_abs, bg = "transparent", width =  8, height = 6, units = "in", dpi = 600)


Rec_path_abs<- Rec_path_abs_SGS+Rec_path_abs_CHY+Rec_path_abs_HYS+Rec_path_abs_KNZ
Rec_path_abs
ggsave(filename = "Rec_path_abs.pdf", plot = Rec_path_abs, bg = "transparent", width =  8, height = 6, units = "in", dpi = 600)


#### Ave ABS cover by period 

#Drought period 
comppath.plotsum_Drt <- comppath.plotsum %>% dplyr::filter(Year %in% c("2014","2015","2016","2017"))

path_ave_Drt<-comppath.plotsum_Drt %>%
  dplyr::group_by(Site,Trt,Path)%>%
  dplyr::summarize(ave = mean(pathsum,na.rm=T), n = n(),sd = sd(pathsum,na.rm=T), se = sd/sqrt(n),
                   LowerCI = ave - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = ave + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

path_ave_Drt$Trt <- factor(path_ave_Drt$Trt, levels=c('con', 'chr', 'int'))
path_ave_Drt$Site <- factor(path_ave_Drt$Site, levels=c('SGS', 'CHY', 'HYS','KNZ'))

Drt_ave_abs<-ggplot(path_ave_Drt,aes(x=Site, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Site")+
  facet_wrap(~Path)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
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
Drt_ave_abs

ggsave(filename = "Drt_ave_abs.pdf", plot = Drt_ave_abs, bg = "transparent", width =  8, height = 6, units = "in", dpi = 600)

#Recovery period 
comppath.plotsum_Rec <- comppath.plotsum %>% dplyr::filter(Year %in% c("2018","2019","2020","2021"))

path_ave_Rec<-comppath.plotsum_Rec %>%
  dplyr::group_by(Site,Trt,Path)%>%
  dplyr::summarize(ave = mean(pathsum,na.rm=T), n = n(),sd = sd(pathsum,na.rm=T), se = sd/sqrt(n),
                   LowerCI = ave - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = ave + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

path_ave_Rec$Trt <- factor(path_ave_Rec$Trt, levels=c('con', 'chr', 'int'))
path_ave_Rec$Site <- factor(path_ave_Rec$Site, levels=c('SGS', 'CHY', 'HYS','KNZ'))

Rec_ave_abs<-ggplot(path_ave_Rec,aes(x=Site, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Site")+
  facet_wrap(~Path)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
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
Rec_ave_abs

ggsave(filename = "Rec_ave_abs.pdf", plot =Rec_ave_abs, bg = "transparent", width =  8, height = 6, units = "in", dpi = 600)

