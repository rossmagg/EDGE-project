
### C3/C4 abundance change by period ###

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

#calculate aves

path_drt_ave<-abunchngtrt%>%
  dplyr::group_by(Site,Trt,Path)%>%
  dplyr::summarize(mean_change = mean(change,na.rm=T), n = n(),sd = sd(change,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_change - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_change + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

path_drt_ave$Trt <- factor(path_drt_ave$Trt, levels=c('con', 'chr', 'int'))

#figure
path_drt_fig<-ggplot(path_drt_ave, aes(x=Site, y=mean_change,color=Trt))+
  ylab("Mean change during drought")+
  xlab("Site")+
  geom_hline(yintercept = 0, color="grey")+
  facet_wrap(~Path, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.40))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.40), width=0, linewidth=1.25, show.legend = TRUE)+
  scale_x_discrete(labels=c('SGS', 'HPG', 'HYS','KNZ'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                        axis.text.y = element_text(size=12, color = "black"),
                        strip.text = element_text(size=12),
                        axis.title.x = element_text(size=14),
                        axis.title.y = element_text(size=14),
                        legend.position = "top",
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())
path_drt_fig
ggsave(filename = "path_drt_fig.pdf", plot = path_drt_fig, bg = "transparent", width =  8, height = 6, units = "in", dpi = 600)

#### C3C4 change during recovery, compared to 2017 ####

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

path_rec_ave<-abunchngtrt.rec%>%
  dplyr::group_by(Site,Trt,Path)%>%
  dplyr::summarize(mean_change = mean(change,na.rm=T), n = n(),sd = sd(change,na.rm=T), se = sd/sqrt(n),
                   LowerCI = mean_change - qt(1 - (0.05 / 2), n - 1) * se,
                   UpperCI = mean_change + qt(1 - (0.05 / 2), n - 1) * se)%>%
  as_tibble()

path_rec_ave$Trt <- factor(path_rec_ave$Trt, levels=c('con', 'chr', 'int'))

#figure
path_rec_fig<-ggplot(path_rec_ave, aes(x=Site, y=mean_change,color=Trt))+
  ylab("Mean change during recovery")+
  xlab("Site")+
  geom_hline(yintercept = 0, color="grey")+
  facet_wrap(~Path, scales = "free", ncol = 2)+
  scale_color_manual(name="Treatment",values=c("#56B4E9","#009E73","#E69F00"),labels = c("Control", "Chronic","Intense"))+
  geom_point(size=4,position=position_dodge(.40))+
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI),position=position_dodge(.40), width=0, linewidth=1.25, show.legend = TRUE)+
  scale_x_discrete(labels=c('SGS', 'HPG', 'HYS','KNZ'))+
  theme_bw()+theme(strip.background =element_rect(fill="lightgrey"),axis.text.x = element_text(size=12, color="black"),
                        axis.text.y = element_text(size=12, color = "black"),
                        strip.text = element_text(size=12),
                        axis.title.x = element_text(size=14),
                        axis.title.y = element_text(size=14),
                        legend.position = "top",
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())
path_rec_fig
ggsave(filename = "path_rec_fig.pdf", plot = path_rec_fig, bg = "transparent", width =  8, height = 6, units = "in", dpi = 600)

