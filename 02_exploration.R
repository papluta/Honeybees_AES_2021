library(dplyr)
library(ggplot2)
library(tidyverse)
library(readr)
library(ggpubr)

# load full raw honey bee colony development data
load('Data/Data_bees_raw.RData')

# load hive mortality
mort <- read.csv("Data/wintermortality2021_binom.csv")
mort <- mort[-c(101:104),] # delete N1145 (had to be removed from the site per farmer's request)

# load land-use
load('Data/data_all_2km.RData')

# COLONY DELEVOPMENT PLOTS
hb_pop <- subset(data.bees.raw, !is.na(pop_round)) # just the rounds where colony size was evaluated
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
total_varroa <- data.bees.raw %>% drop_na(varroa_board) %>% group_by(sample) %>% summarise(tvar=sum(varroa_board))
total_varroa$site <- as.factor(sapply(strsplit(as.character(total_varroa$sample), split='V', fixed=TRUE),function(x) (x[1])))

# GENERAL VARROA PLOTS
hb_var <- data.bees.raw %>% drop_na(v_round) %>% drop_na(daily_varroa)  # deleting missing values
ggplot(data = hb_var, aes(x = factor(v_round), y = log(varroa_board+1), fill = as.factor(treated_previously))) + 
  stat_boxplot(geom = "errorbar", width = 0.5) + 
  geom_boxplot() + 
  xlab("round") + 
  ylab("log mite fall") + 
  theme_bw(base_size = 12, base_family = "sans")+
  labs(fill = 'Treated with acaricide')+
  theme(legend.position = 'bottom',legend.direction = "vertical")
#ggsave('figures/SM_varroa_ex.pdf', width = 10, height = 10)


#### PREVALENCE PLOT ####
library(ggpattern)


load('Data/Data_final.RData')
parplot <- data.final.august %>% bind_rows(data.final.july) %>% select(Site, month, dwva:nob) %>% 
  pivot_longer(cols = dwva:nob, names_to = 'Parasite', values_to = 'Status') %>% group_by(Parasite, month) %>%
  summarise(Prevalence = mean(Status))
  


type <- c('dwva', 'dwvb', 'bqcv', 'sbv', 'abpv', 'cbpv', 'sbpv','tryp', 'neo', 'noc', 'nob')

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#316288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", 'yellow','azure4')

custom <- c('dwva'='DWV-A', 'dwvb'='DWV-B','bqcv'='BQCV','sbv'='SBV','abpv'='ABPV', 'cbpv' = 'CBPV', 'sbpv' = 'SBPV','tryp'='Tryp.', 'neo'='Neog.', 'noc'='V.ceranae', 'nob'='V.bombi')
#fill <- c("#DDCC77","#DDCC77","#117733","#117733","#316288","#316288","#AA4499","#AA4499","#44AA99","#44AA99","#999933","#999933","#661100","#661100","#88CCEE","#88CCEE","#CC6677","#CC6677")
pal <- c('#BDD9BF', '#929084', '#FFC857', '#A997DF', '#E5323B', '#77b4dd','#DDCC77', '#E1DAAE','#117733', '#316288', '#FFC857')
pal2 <- rep(pal, each=2)

### alpha
ggplot(parplot, aes(x=factor(Parasite, level=type), y=Prevalence, alpha=factor(month, levels = c('July', 'August')))) +
  geom_bar(stat = "identity", position = 'dodge', colour  = 'black', fill = pal2) +
  theme_bw() +
  scale_alpha_discrete(range = c(0.5, 1)) +
  xlab(NULL) +
  ggtitle(' ')+
  ylab("Parasite prevalence") +
  theme(axis.ticks.x = element_blank(),legend.text = element_text(size=18), legend.title = element_text(size=22),title = element_text(size=24), axis.text=element_text(size=20, colour = 'black'), axis.title=NULL, axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)) + 
  scale_x_discrete(labels=custom) +
  scale_y_continuous(expand = c(0, 0.005), limits = c(0, 1))+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(color='black'))

#ggsave('Figures/Prevalence.pdf', height = 18, width = 33, units = 'cm')

#pattern
ggplot(parplot, aes(x=factor(Parasite, level=type), y=Prevalence)) +
  geom_bar_pattern(stat = "identity", aes(pattern = factor(month, levels = c('July', 'August'))), position = 'dodge',
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

