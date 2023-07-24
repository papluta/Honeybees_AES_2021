library(dplyr)
library(ggplot2)
library(tidyverse)
library(readr)
library(ggpubr)

# load full raw honey bee colony development data
coldev <- read.csv('hb_colony_development.csv')

# load hive mortality
mort <- read.csv("wintermortality2021_binom.csv")
mort <- mort[-c(101:104),] # delete N1145 (had to be removed from the site per farmer's request)

# load land-use
load('Landuse.RData')

# COLONY DELEVOPMENT PLOTS
hb_pop <- subset(coldev, !is.na(pop_round)) # just the rounds where colony size was evaluated
bee <- ggplot(data = hb_pop, aes(x = factor(pop_round), y = no_bees)) + 
  stat_boxplot(geom = "errorbar", width = 0.5) + 
  geom_boxplot() + 
  xlab("round") + 
  ylab("absolute bee number") + 
  theme_bw(base_size = 14, base_family = "sans")
brood <- ggplot(data = hb_pop, aes(x = factor(pop_round), y = no_brood)) + 
  stat_boxplot(geom = "errorbar", width = 0.5) + 
  geom_boxplot() + 
  xlab("round") + 
  ylab("brood cell number") + 
  theme_bw(base_size = 14, base_family = "sans")
comb <- ggplot(data = hb_pop, aes(x = factor(pop_round), y = bees_comb)) + 
  stat_boxplot(geom = "errorbar", width = 0.5) + 
  geom_boxplot() + 
  xlab("round") + 
  ylab("bee-covered area of combs") + 
  theme_bw(base_size = 14, base_family = "sans")
weight <- ggplot(data = hb_pop, aes(x = factor(pop_round), y = weight)) + 
  stat_boxplot(geom = "errorbar", width = 0.5) + 
  geom_boxplot() + 
  xlab("round") + 
  ylab("weight") + 
  theme_bw(base_size = 14, base_family = "sans")
ggarrange(bee, weight, brood, comb)
ggsave('figures/SM_colonydev_ex.pdf', width = 10, height = 10)

#total acaricide mite fall
total_varroa <- subset(coldev, treated_previously=='1') #mite fall after treatment
total_varroa <- subset(total_varroa, !is.na(varroa_board))
total_varroa <- total_varroa %>% group_by(sample) %>% summarise(tvar=sum(varroa_board))
total_varroa$site <- as.factor(sapply(strsplit(as.character(total_varroa$sample), split='V', fixed=TRUE),function(x) (x[1])))
#dotchart(total_varroa$tvar, groups=total_varroa$site)

# GENERAL VARROA PLOTS
hb_var <- subset(coldev, !is.na(v_round))  # choosing only rounds where varroa fall was assessed
hb_var <- subset(hb_var, !is.na(daily_varroa))  # deleting missing values
ggplot(data = hb_var, aes(x = factor(v_round), y = log(varroa_board+1), fill = as.factor(treated_previously))) + 
  stat_boxplot(geom = "errorbar", width = 0.5) + 
  geom_boxplot() + 
  xlab("round") + 
  ylab("log mite fall") + 
  theme_bw(base_size = 12, base_family = "sans")+
  labs(fill = 'Treated with acaricide')+
  theme(legend.position = 'bottom',legend.direction = "vertical")
ggsave('figures/SM_varroa_ex.pdf', width = 10, height = 10)


#CHOOSING 4TH ROUND FOR POPULATION MODELS
dev4 <- subset(coldev, pop_round=='4') #at the time of pathogen sampling
dev4 <- dev4[-c(6)] #deleting 'site' column to avoid duplication
dev4 <- subset(dev4, !is.na(no_bees))
dev4 <- subset(dev4, !is.na(daily_varroa))

dev4 <- merge(dev4, total_varroa, by='sample', all.x = T)
dev4 <- merge(dev4, land, by='site', all.x = T)
dev4 <- merge(dev4, mort, by='sample', all.x = T)
dev4 <- subset(dev4, !site=='N1145') #deleting N1145 because it was removed from the site


# MORTALITY PLOTS
hist(mort$status)

C <- ggplot(dev4, aes(x=CropBV1_perc,y=status)) + 
  geom_point() + 
  geom_smooth(method="lm") +
  xlab('Organic crop %') +
  ylab('Survival') +
  theme_classic()
S <- ggplot(dev4, aes(x=semi_perc,y=status)) + 
  geom_point() + 
  geom_smooth(method="gam") +
  xlab('Semi-natural habitat %') +
  ylab('Survival') +
  theme_classic()
F <- ggplot(dev4, aes(x=Annual_flower,y=status)) + 
  geom_point() + 
  geom_smooth(method="lm") +
  xlab('Annual flower ha') +
  ylab('Survival') +
  theme_classic()
D <- ggplot(dev4, aes(x=no_bees,y=status)) + 
  geom_point() + 
  geom_smooth(method="gam") +
  xlab('Number of bees') +
  ylab('Survival') +
  theme_classic()
ggarrange(C,S,F,D)

# mortality and varroa
dev4$tvar100 <- dev4$tvar / dev4$no_bees * 100 #standardizing treatment fall per 100 bees
ggplot(dev4, aes(x=tvar100,y=status)) + 
  geom_point() + 
  geom_smooth(method="lm") +
  xlab('Treatment fall per 100 bees') +
  ylab('Survival') +
  theme_classic()

# adding parasites - subset ONLY

parasite <- read.table("HBULK.csv", 
                       header=TRUE, stringsAsFactors=TRUE, sep=",", na.strings="NA", dec=".", 
                       strip.white=TRUE, fill=T)
parasite <- merge(parasite, land, by='site')
parasite <- parasite[-c(16:34)]
parasite$Z.CropBV1_perc <- scale(parasite$CropBV1_perc)
parasite$Z.Annual_flower <- scale(parasite$Annual_flower)
parasite$Z.semi_perc <- scale(parasite$semi_perc)
parasite$Z.p_rich <- scale(parasite$p_rich)
parasite1 <- subset(parasite, date == 'first') #separating 1st and 2nd sampling
parasite2 <- subset(parasite, date == 'second')
colnames(parasite2)[(14)] <- c('p_rich2') 
colnames(parasite2)[(23)] <- c('Z.p_rich2')
parasite2$p_rich1 <- parasite1$p_rich #adding 1st p_rich to the df for the models
parasite2$Z.p_rich1 <- parasite1$Z.p_rich
parasite2 <- parasite2[-c(3:4)]
pardev <- merge(parasite2, dev4, by='sample', all.x=T)

ggplot(pardev, aes(x=p_rich2,y=status)) + 
  geom_point() + 
  geom_smooth(method="lm") +
  xlab('Parasite species richness') +
  ylab('Survival') +
  theme_classic()

ggplot(pardev, aes(x=p_rich2,y=no_bees)) + 
  geom_point() + 
  geom_smooth(method="lm") +
  xlab('Parasite species richness') +
  ylab('Number of bees') +
  theme_classic()

Cp <- ggplot(pardev, aes(x=CropBV1_perc,y=p_rich2)) + 
  geom_point() + 
  geom_smooth(method="lm") +
  xlab('Organic crop %') +
  ylab('Parasite species richness') +
  theme_classic()
Sp <- ggplot(pardev, aes(x=semi_perc,y=p_rich2)) + 
  geom_point() + 
  geom_smooth(method="lm") +
  xlab('Semi-natural habitat %') +
  ylab('Parasite species richness') +
  theme_classic()
Fp <- ggplot(pardev, aes(x=Annual_flower,y=p_rich2)) + 
  geom_point() + 
  geom_smooth(method="lm") +
  xlab('Annual flower ha') +
  ylab('Parasite species richness') +
  theme_classic()
ggarrange(Cp,Sp,Fp)

ggplot(parasite, aes(x=site, y=p_rich, fill=date))+
  geom_bar(stat='identity')+
  facet_wrap(~date)

HBULK_prev <- read.table('C:/Users/patry/OneDrive/PhD/tralala/combee/dfs/HBULK_prev.csv', header = T, strip.white = T, stringsAsFactors = T, dec='.', sep = ',')
HBULK_prev$noc1 <- rep(0, times=16) #there was no Nosema ceranae and bombi so I did not include them in the final CSV, that's why now I have to add the zeros
HBULK_prev$nob1 <- rep(0, times=16)
HBULK_prev$noc2 <- HBULK_prev$noc1
HBULK_prev$nob2 <- HBULK_prev$nob1
dwva1 <- sum(HBULK_prev$dwva1)/16
HBULK_gprev <- as.data.frame(c(1))
HBULK_gprev$dwva1 <- sum(HBULK_prev$dwva1)/16
HBULK_gprev$dwva2 <- sum(HBULK_prev$dwva2)/16
HBULK_gprev$dwvb1 <- sum(HBULK_prev$dwvb1)/16
HBULK_gprev$dwvb2 <- sum(HBULK_prev$dwvb2)/16
HBULK_gprev$bqcv1 <- sum(HBULK_prev$bqcv1)/16
HBULK_gprev$bqcv2 <- sum(HBULK_prev$bqcv2)/16
HBULK_gprev$sbv1 <- sum(HBULK_prev$sbv1)/16
HBULK_gprev$sbv2 <- sum(HBULK_prev$sbv2)/16
HBULK_gprev$abpv1 <- sum(HBULK_prev$abpv1)/16
HBULK_gprev$abpv2 <- sum(HBULK_prev$abpv2)/16
HBULK_gprev$tryp1 <- sum(HBULK_prev$tryp1)/16
HBULK_gprev$tryp2 <- sum(HBULK_prev$tryp2)/16
HBULK_gprev$neo1 <- sum(HBULK_prev$neo1)/16
HBULK_gprev$neo2 <- sum(HBULK_prev$neo2)/16
HBULK_gprev$noc1 <- sum(HBULK_prev$noc1)/16
HBULK_gprev$noc2 <- sum(HBULK_prev$noc2)/16
HBULK_gprev$nob1 <- sum(HBULK_prev$nob1)/16
HBULK_gprev$nob2 <- sum(HBULK_prev$nob2)/16
HBULK_gprev <- HBULK_gprev[-c(1)]
HBULK_gprev <- unlist(HBULK_gprev)
HBULK_gprev <- as.data.frame(HBULK_gprev)
HBULK_gprev$parasite <- rownames(HBULK_gprev)
rownames(HBULK_gprev) <- NULL
colnames(HBULK_gprev) <- c('Prevalence', 'Parasite')

type <- c('dwva', 'dwvb', 'bqcv', 'sbv', 'abpv', 'tryp', 'neo', 'noc', 'nob')
HBULK_gprev$type <- rep(type, times=1, each=2)
date <- c('first', 'second')
HBULK_gprev$date <- rep(date, times=9)
HBULK_gprev$date <- as.factor(HBULK_gprev$date)
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#316288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", 'yellow','azure4')
library(ggpattern)
library(viridis)
library(ggplot2)
lvl_order <- c('dwva1','dwva2','dwvb1','dwvb2','bqcv1','bqcv2','sbv1','sbv2','abpv1','abpv2','tryp1','tryp2','neo1','neo2', 'noc1','noc2','nob1','nob2')
custom <- c('dwva1'='DWV-A', 'dwva2'='', 'dwvb1'='DWV-B','dwvb2'='','bqcv1'='BQCV','bqcv2'='','sbv1'='SBV','sbv2'='','abpv1'='ABPV','abpv2'='','tryp1'='Tryp.','tryp2'='','neo1'='Neog.','neo2'='', 'noc1'='N.ceranae', 'noc2'='', 'nob1'='N.bombi', 'nob2'='')
fill <- c("#DDCC77","#DDCC77","#117733","#117733","#316288","#316288","#AA4499","#AA4499","#44AA99","#44AA99","#999933","#999933","#661100","#661100","#88CCEE","#88CCEE","#CC6677","#CC6677")
pal <- c('#BDD9BF', '#929084', '#FFC857', '#A997DF', '#E5323B', '#2E4052','#058ED9', '#E1DAAE','#E1DAAE')
pal2 <- rep(pal, each=2)
library(RColorBrewer)
library(ggtext)

#pattern
ggplot(HBULK_gprev, aes(x=factor(Parasite, level=lvl_order), y=Prevalence)) +
  geom_bar_pattern(stat = "identity", aes(pattern = date), position = 'dodge',
                   pattern_fill = 'white', colour  = 'black', fill=pal2) +
  scale_pattern_manual(values=c("none", "stripe")) +
  #scale_fill_manual(values= '#88CCEE')
  theme_classic() +
  #scale_pattern_fill_manual(values = c('#88CCEE', "#88CCEE", '#CC6677', '#CC6677', '#DDCC77', '#DDCC77', '#117733', '#117733', '#316288', '#316288', '#AA4499', '#AA4499', '#44AA99', '#44AA99')) +
  #labs(title= 'Honey bees') +
  xlab(NULL) +
  ylab("Parasite prevalence") +
  guides(pattern=guide_legend(title= 'Sampling time', title.theme = element_text(size=18)), fill='none')+
  theme(axis.ticks.x = element_blank(),legend.text = element_text(size=18), legend.title = element_text(size=22),title = element_text(size=24), axis.text=element_text(size=20, colour = 'black', vjust = 0.1), axis.title=NULL, axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)) + #face="bold"
  #ylim(0:1) +
  scale_x_discrete(labels=custom) +
  scale_y_continuous(expand = c(0, 0.005), limits = c(0, 1))

### alpha
HBULK_gprev2 <- subset(HBULK_gprev, date == 'second')
ggplot(HBULK_gprev, aes(x=factor(Parasite, level=lvl_order), y=Prevalence, alpha=date)) +
  geom_bar(stat = "identity", position = 'dodge', colour  = 'black', fill=pal2) +
  theme_bw() +
  scale_alpha_discrete(range = c(0.5, 1)) +
  xlab(NULL) +
  ggtitle(' ')+
  ylab("Parasite prevalence") +
  guides(pattern=guide_legend(title= 'Sampling time', title.theme = element_text(size=18)), fill='none')+
  theme(axis.ticks.x = element_blank(),legend.text = element_text(size=18), legend.title = element_text(size=22),title = element_text(size=24), axis.text=element_text(size=20, colour = 'black', vjust = 0.1), axis.title=NULL, axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)) + #face="bold"
  #ylim(0:1) +
  scale_x_discrete(labels=custom) +
  scale_y_continuous(expand = c(0, 0.005), limits = c(0, 1))+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(color='black'), axis.text.x = element_text(hjust = 0.1))
