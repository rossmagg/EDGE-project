
library(tidyverse)
library(ggplot2)

growpre <- read.csv("C:\\Users\\Maggie Ross\\OneDrive - Colostate\\EDGE\\data\\precipitation data\\final_dfs\\growingseason_precip_totals_allyears.csv")

str(growpre.sum)

growpre.sum<- growpre%>%
  dplyr::group_by(Site,Year) %>%
  dplyr::summarize(totalprecip = sum(ambient_precip),n=n()) %>%
  as_tibble()


## ADD 2013 DATA point from 30yr sums
#these have been updated to include april and sept 1-15

#sgs
growpre.sum <- rbind(growpre.sum, c("SGS", 2013,228.0))
growpre.sum <- rbind(growpre.sum, c("CHY", 2013,258.3))
growpre.sum <- rbind(growpre.sum, c("HYS", 2013,390.1))
growpre.sum <- rbind(growpre.sum, c("KNZ", 2013,574.8))                     

write.csv(growpre.sum,"C:\\Users\\Maggie Ross\\OneDrive - Colostate\\EDGE\\data\\precipitation data\\precip_growing_sums_MR.csv", row.names = FALSE)
sums_edited <- read.csv("C:\\Users\\Maggie Ross\\OneDrive - Colostate\\EDGE\\data\\precipitation data\\precip_growing_sums_MR.csv")

## use this df ## 
sums_edited_drought <- read.csv("C:\\Users\\Maggie Ross\\OneDrive - Colostate\\EDGE\\data\\precipitation data\\data\\manipulated_outputs_etc\\precip_growing_sums_droughtreduction.csv")
#switch to long
sums_edited_drought_long <- sums_edited_drought %>%
  pivot_longer(cols=totalprecip:drought_precip,names_to = "precip_amb_drought", values_to = "total")

#####################################

#30 yr ave bars 

#SGS 

#this just includes half of sept data.. need to go back and get daily totals 
SGS <- read.csv("C:\\Users\\Maggie Ross\\OneDrive - Colostate\\EDGE\\data\\precipitation data\\final_dfs\\SGS_30yr_growing.csv")

SGS.ave<- SGS%>%
  dplyr::summarize(aveprecip_30yr = mean(growing)) %>%
  as_tibble()

#sgs fig with bar

str(sums_edited_drought)

windowsFonts(A = windowsFont("Arial"))

SGSprecip_bar<-ggplot(subset(sums_edited_drought_long, Site=="SGS"), aes(x=Year, y=total, fill=precip_amb_drought))+
  ggtitle("a) SGS")+
  ylab("Growing season precipitation (mm)")+xlab("Year")+
  geom_bar(stat="identity",position = "identity")+
  scale_fill_brewer(palette="Paired")+
  scale_fill_manual(name="Treatment",values=c("#CCCFFF","#3333FF"),labels = c("Drought", "Ambient"))+
  scale_y_continuous(expand=c(0,0), breaks = c(0,50,100,150,200,250,300), limits = c(0,300))+
  theme_classic()+
  annotate("rect", xmin=2013.5,xmax=2017.5, ymin=-Inf,ymax = Inf, alpha=.1, fill="#999999")+
  annotate("rect", xmin=2017.501,xmax=2021.5, ymin=-Inf,ymax = Inf, alpha=.1, fill="#99CCFF")+
  scale_x_continuous(breaks = c(2013,2014,2015,2016,2017,2018,2019,2020,2021),expand=c(.01,.01))+
  geom_hline(yintercept = 262.9, color="black", linetype="dashed", size=1)+
  theme(legend.position="none",strip.background = element_blank(),legend.title = element_blank(),axis.text.x=element_blank(), 
        axis.text.y=element_text(size=11, color="black"),plot.title = element_text(size=20),
        text = element_text(family = "A"),axis.title.x = element_blank(),axis.title.y = element_blank())                                                                                                                                   
SGSprecip_bar
ggsave(filename = "SGSprecip_bar.jpeg", plot = SGSprecip_bar, bg = "transparent", width =  9, height = 6, units = "in", dpi = 600)


#CHY 

CHY <- read.csv("C:\\Users\\Maggie Ross\\OneDrive - Colostate\\EDGE\\data\\precipitation data\\final_dfs\\CHY_growing_30yrs.csv")

CHY.sum<- CHY%>%
  dplyr::group_by(Year) %>%
  dplyr::summarize(totalprecip = sum(PRCP)) %>%
  as_tibble()

CHY.ave <- CHY.sum %>%
  dplyr::summarize(aveprecip_30yr = mean(totalprecip)) %>%
  as_tibble()


CHYprecip_bar<-ggplot(subset(sums_edited_drought_long, Site=="CHY"), aes(x=Year, y=total, fill=precip_amb_drought))+
  ggtitle("b) HPG")+
  ylab("Growing season precipitation (mm)")+xlab("Year")+
  geom_bar(stat="identity",position = "identity")+
  scale_fill_manual(name="Treatment",values=c("#CCCFFF","#3333FF"),labels = c("Drought", "Ambient"))+
  scale_y_continuous(expand=c(0,0), breaks = c(0,50,100,150,200,250,300), limits = c(0,300))+
  theme_classic()+
  annotate("rect", xmin=2013.5,xmax=2017.5, ymin=-Inf,ymax = Inf, alpha=.1, fill="#999999")+
  annotate("rect", xmin=2017.501,xmax=2021.5, ymin=-Inf,ymax = Inf, alpha=.1, fill="#99CCFF")+
  scale_x_continuous(breaks = c(2013,2014,2015,2016,2017,2018,2019,2020,2021),expand=c(.01,.01))+
  geom_hline(yintercept = 262.9, color="black", linetype="dashed", size=1)+
  theme(legend.position="none",strip.background = element_blank(),legend.title = element_blank(),axis.text.x=element_blank(), 
        axis.text.y=element_text(size=11, color="black"),axis.title=element_text(size=12),plot.title = element_text(size=20),
        text = element_text(family = "A"),axis.title.x = element_blank(),axis.title.y = element_blank())                                                                                                                                          
CHYprecip_bar


#HYS

HYS <- read.csv("C:\\Users\\Maggie Ross\\OneDrive - Colostate\\EDGE\\data\\precipitation data\\final_dfs\\HYS_growing_30yrs.csv")

HYS.sum<- HYS%>%
  dplyr::group_by(Year) %>%
  dplyr::summarize(totalprecip = sum(PRCP)) %>%
  as_tibble()

HYS.ave <- HYS.sum %>%
  dplyr::summarize(aveprecip_30yr = mean(totalprecip)) %>%
  as_tibble()

#HYS fig with bar 

HYSprecip_bar<-ggplot(subset(sums_edited_drought_long, Site=="HYS"), aes(x=Year, y=total, fill=precip_amb_drought))+
  ggtitle("c) HYS")+
  ylab("Growing season precipitation (mm)")+xlab("Year")+
  geom_bar(stat="identity", position = "identity")+
  scale_y_continuous(expand=c(0,0), breaks = c(0,100,200,300,400,500,600,700), limits = c(0,700))+theme_classic()+
  scale_fill_manual(name="Treatment",values=c("#CCCFFF","#3333FF"),labels = c("Drought", "Ambient"))+
  annotate("rect", xmin=2013.5,xmax=2017.5, ymin=-Inf,ymax = Inf, alpha=.1, fill="#999999")+
  annotate("rect", xmin=2017.501,xmax=2021.5, ymin=-Inf,ymax = Inf, alpha=.1, fill="#99CCFF")+
  scale_x_continuous(breaks = c(2013,2014,2015,2016,2017,2018,2019,2020,2021),expand=c(.01,.01))+
  geom_hline(yintercept = 394.5, color="black", linetype="dashed", size=1)+
  theme(legend.position="none",strip.background = element_blank(),legend.title = element_blank(),axis.text.x=element_text(size=11, color="black"), 
        axis.text.y=element_text(size=11, color="black"),axis.title=element_text(size=12),,plot.title = element_text(size=20),
        text = element_text(family = "A"),axis.title.x = element_blank(),axis.title.y = element_blank())                                                                                                                                 
HYSprecip_bar


#KNZ

KNZ <- read.csv("C:\\Users\\Maggie Ross\\OneDrive - Colostate\\EDGE\\data\\precipitation data\\final_dfs\\KNZ_growing_30yrs.csv")

KNZ.sum<- KNZ%>%
  dplyr::group_by(Year) %>%
  dplyr::summarize(totalprecip = sum(ppt)) %>%
  as_tibble()

#remove row with na
KNZ.sum[is.na(KNZ.sum)] <- 0
KNZ.sum<-KNZ.sum[KNZ.sum$totalprecip != 0, ]

KNZ.ave <- KNZ.sum %>%
  dplyr::summarize(aveprecip_30yr = mean(totalprecip)) %>%
  as_tibble()

#KNZ fig with bar 

KNZprecip_bar<-ggplot(subset(sums_edited_drought_long, Site=="KNZ"), aes(x=Year, y=total, fill=precip_amb_drought))+
  ggtitle("d) KNZ")+
  ylab("Growing season precipitation (mm)")+xlab("Year")+
  geom_bar(stat="identity", position = "identity")+
  scale_fill_manual(name="Treatment",values=c("#CCCFFF","#3333FF"),labels = c("Drought", "Ambient"))+
  scale_y_continuous(expand=c(0,0), breaks = c(0,200,400,600,800,1000), limits = c(0,1000))+theme_classic()+
  annotate("rect", xmin=2013.5,xmax=2017.5, ymin=-Inf,ymax = Inf, alpha=.1, fill="#999999")+
  annotate("rect", xmin=2017.501,xmax=2021.5, ymin=-Inf,ymax = Inf, alpha=.1, fill="#99CCFF")+
  scale_x_continuous(breaks = c(2013,2014,2015,2016,2017,2018,2019,2020,2021),expand=c(.01,.01))+
  geom_hline(yintercept = 611.42, color="black", linetype="dashed", size=1)+
  theme(legend.position="none",strip.background = element_blank(),legend.title = element_blank(),axis.text.x=element_text(size=11, color="black"), 
        axis.text.y=element_text(size=11, color="black"),axis.title=element_text(size=11),plot.title = element_text(size=20),
        text = element_text(family = "A"),axis.title.x = element_blank(),axis.title.y = element_blank())                                                                                                                                         
KNZprecip_bar

library(patchwork)

precip_bar_30 <- SGSprecip_bar+CHYprecip_bar+HYSprecip_bar+KNZprecip_bar
precip_bar_30
ggsave(filename = "Precip_30yr.jpeg", plot = precip_bar_30, bg = "transparent", width =  9, height = 6, units = "in", dpi = 600)