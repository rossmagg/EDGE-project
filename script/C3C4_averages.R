comp_zeros <- read.csv("C:\\Users\\Maggie Ross\\OneDrive - Colostate\\EDGE\\data\\spcompdata_final_revised_April2022\\spcomp_zeros_trt.csv")
comp_zeros$Species = NULL
colnames(comp_zeros)[colnames(comp_zeros) == "max.cover"] <- "cover" 
comp_zeros<- comp_zeros %>%
  dplyr::filter(Year %in% c(2013,2014,2015,2016,2017,2018,2019,2020,2021)) 

#add zeros for averaging 
comp.wide <- comp_zeros %>%
  pivot_wider(names_from = Year, values_from = cover)
comp.wide[is.na(comp.wide)] <- 0

#back to long for calculating 
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

trt <- read.csv("C:\\Users\\Maggie Ross\\OneDrive - Colostate\\EDGE\\data\\trt_info_2.csv")
trt<-trt %>% 
  dplyr::mutate_at(c('Site','Plot','Block','Trt'),as.factor)
comp.aveplot_trt<- merge(trt,comp.aveplot_zeros, fix.by=c("Site", "Plot"))

comp.species <- comp.aveplot_trt %>% dplyr::filter(Spcode %in% c("45","46","47","48","49","55","130","137","197","317","95","44", "9","7","89","221","235","234", "204","50","53","194","51","14","88","96","388","94","56","57","325","172"))
comp.species$Spcode <- factor(comp.species$Spcode, levels = c("45","46","47","48","49","55","130","137","197","317","95","44","9","7","89","221","235","234", "204","50","53","194","51","14","88","96","388","94","56","57","325","172"), 
                              labels = c("BODA","BOGR","BOHI","BRJA","BRTE","CAEL","HECO","KOPY","PASM","VUOC","ELEL","BOCU","ANTE","ANGE","DIOL","SCSC","SPAS","SONU","POPR","CABR","CAME","PAVI","CAHE","ARPU","DIAC","ERSP","PACA","ELSP","CAFI","CAGR","CABL","MUTO"))

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
    Spcode == "CABL" ~ "C3",
    Spcode == "MUTO" ~ "C4"))

comppath.plotsum <- comp.path %>%
  dplyr::group_by(Site,Year,Block,Trt,Plot,Path) %>% 
  dplyr::summarise(pathsum = sum(avg.cover))

comppath.plotsum$Trt <- factor(comppath.plotsum$Trt, levels=c('con', 'chr', 'int'))
comppath.plotsum$Site <- factor(comppath.plotsum$Site, levels=c('SGS', 'CHY', 'HYS','KNZ'))

path.ave.comb<-comppath.plotsum%>%
  dplyr::group_by(Site,Trt,Path)%>%
  dplyr::summarize(mean_cover = mean(pathsum,na.rm=T), n = n(),sd = sd(pathsum,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_cover - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_cover + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

#fig of con aves with all years 
con.AVES.FIG <- ggplot(subset(path.ave.comb, Trt=="con"), aes(x=Site, y=mean_cover))+ 
  facet_wrap(~Path,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  geom_bar(stat="identity", position = "dodge", width = .9)+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI, group=Trt),position=position_dodge(.9),width=.1)+ 
  #scale_y_continuous(expand=c(0,0))+
  theme_classic()+theme(legend.position="none",legend.title = element_blank())
con.AVES.FIG

ggsave(filename = "c3c4ave.allyrs.jpeg", plot = con.AVES.FIG, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)


#fig of con aves with just pretreatment year
comppath.plotsum.2013<- comppath.plotsum %>%
  dplyr::filter(Year == 2013)


path.ave.2013<-comppath.plotsum.2013%>%
  dplyr::group_by(Site,Trt,Path)%>%
  dplyr::summarize(mean_cover = mean(pathsum,na.rm=T), n = n(),sd = sd(pathsum,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_cover - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_cover + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

con.AVES.FIG.2013 <- ggplot(subset(path.ave.2013, Trt=="con"), aes(x=Site, y=mean_cover))+ 
  facet_wrap(~Path,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  geom_bar(stat="identity", position = "dodge", width = .9)+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI, group=Trt),position=position_dodge(.9),width=.1)+ 
  #scale_y_continuous(expand=c(0,0))+
  theme_classic()+theme(legend.position="none",legend.title = element_blank())
con.AVES.FIG.2013

ggsave(filename = "c3c4ave.2013.jpeg", plot = con.AVES.FIG.2013, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)

#c3 c4 on same bar
c3c4barcomb13 <- ggplot(subset(path.ave.2013, Trt=="con"), aes(x=Site, y=mean_cover, fill=Path))+ 
  #facet_wrap(~Path,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  scale_fill_manual(name="Path",values=c("#C3D7A4","#52854C"),labels = c("C3", "C4,"))+
  geom_bar(stat="identity", position = "stack")+
  #geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI, group=Trt),position=position_dodge(.9),width=.1)+ 
  scale_y_continuous(expand=c(0,0))+
  theme_classic()+theme(legend.position="top",legend.title = element_blank())
c3c4barcomb13 

ggsave(filename = "c3c4ave.2013stack.jpeg", plot = c3c4barcomb13, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)



#mean 2013 cover comb plots ave (not sep by trt)

path.ave.2013.notrt<-comppath.plotsum.2013%>%
  dplyr::group_by(Site,Path)%>%
  dplyr::summarize(mean_cover = mean(pathsum,na.rm=T), n = n(),sd = sd(pathsum,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_cover - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_cover + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

c3c4barcomb13.notrt <- ggplot(path.ave.2013.notrt, aes(x=Site, y=mean_cover, fill=Path))+ 
  #facet_wrap(~Path,strip.position="top",labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free")+
  geom_bar(stat="identity", position = "stack")+
  #geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI, group=Trt),position=position_dodge(.9),width=.1)+ 
  scale_y_continuous(expand=c(0,0))+
  theme_classic()+theme(legend.position="none",legend.title = element_blank())
c3c4barcomb13.notrt


# Panel fig showing C3/C4 cover for drought trts at 2013, 2017, 2021
# facetgrid? wrap by site and year. x=trt; y=ave color, fill=path 

comppath.plotsum_pan<- comppath.plotsum %>%
  dplyr::filter(Year %in% c(2013,2017,2021))

path.ave.pan<-comppath.plotsum_pan%>%
  dplyr::group_by(Site,Year,Trt,Path)%>%
  dplyr::summarize(mean_cover = mean(pathsum,na.rm=T), n = n(),sd = sd(pathsum,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_cover - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_cover + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

path.ave.pan.fig <- ggplot(path.ave.pan, aes(x=Trt, y=mean_cover, fill=Path))+ 
  facet_grid(vars(Site), vars(Year), scales = "free")+
  xlab("Treatment")+
  ylab("Mean cover")+
  geom_bar(stat="identity", position = "stack")+
  scale_fill_manual(name="Pathway",values=c("#C3D7A4","#52854C"),labels = c("C3", "C4"))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), legend.position = "top",legend.title = element_blank(),
                   strip.text = element_text(size = 14))
path.ave.pan.fig

ggsave(filename = "c3c4ave.panel.stack.jpeg", plot = path.ave.pan.fig, bg = "transparent", width =  10, height = 6, units = "in", dpi = 600)


#Ratio

#c3/c4 divided

c3c4_wide<-spread(comppath.plotsum,Path,pathsum)

c3c4_div <- c3c4_wide %>%
  dplyr::group_by(Site,Year, Block, Trt, Plot) %>%
  dplyr::summarise(ratio=sum((C3)/(C4))) %>%
  as_tibble()

c3c4_div.ave <- c3c4_div %>%
  dplyr::group_by(Site,Year,Trt) %>%
  dplyr::summarise(ave.div=mean(ratio)) %>%
  as_tibble()
