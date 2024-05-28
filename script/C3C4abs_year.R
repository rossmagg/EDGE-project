## C3/C4 avg cover over time ##

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


#figs by site 

path_abs_SGS <- ggplot(subset(path_trt_ave, Site=="SGS"),aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year")+
  facet_wrap(~Path, scales = "free")+
  scale_x_discrete(labels=c('0', '1', '2','3','4','5','6','7','8'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=12),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none")
path_abs_SGS

path_abs_CHY <- ggplot(subset(path_trt_ave, Site=="CHY"),aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year")+
  facet_wrap(~Path, scales = "free")+
  scale_x_discrete(labels=c('0', '1', '2','3','4','5','6','7','8'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=12),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none")
path_abs_CHY

path_abs_CHY <- ggplot(subset(path_trt_ave, Site=="CHY"),aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year")+
  facet_wrap(~Path, scales = "free")+
  scale_x_discrete(labels=c('0', '1', '2','3','4','5','6','7','8'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=12),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none")
path_abs_CHY

path_abs_HYS <- ggplot(subset(path_trt_ave, Site=="HYS"),aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year")+
  facet_wrap(~Path, scales = "free")+
  scale_x_discrete(labels=c('0', '1', '2','3','4','5','6','7','8'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=12),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none")
path_abs_HYS

path_abs_KNZ <- ggplot(subset(path_trt_ave, Site=="KNZ"),aes(x=Year, y=ave,color=Trt))+
  ylab("Mean absolute cover")+
  xlab("Year")+
  facet_wrap(~Path, scales = "free")+
  scale_x_discrete(labels=c('0', '1', '2','3','4','5','6','7','8'))+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=3,position=position_dodge(.65))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.65), width=0, size=1, show.legend = TRUE)+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.text.y = element_text(size=12, color = "black"),
                   strip.text = element_text(size=12),
                   axis.title.x = element_text(size=14),
                   axis.title.y = element_text(size=14),
                   legend.position = "none")
path_abs_KNZ

C3C4abs_fig <- ggarrange(path_abs_SGS, path_abs_CHY,path_abs_HYS,path_abs_KNZ)

ggsave(filename = "C3C4abs_fig.jpeg", plot = C3C4abs_fig, bg = "transparent", width =  12, height = 8, units = "in", dpi = 600)