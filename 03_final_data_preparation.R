library(readr)
library(dplyr)
library(tidyverse)


### bee data
load('data_allmonths.RData')

data.cut <- data.merged %>% filter(round %in% c(3,4))

### parasite data
par <- read.csv('HBULK.csv')

data.p <- data.cut  %>% inner_join(par, by = c('month' = 'date', 'sample' = 'sample', 'site' = 'site'))

data.august <- data.p %>% filter(month == 'August') %>%
  mutate(Varroa.p = varroa_board/no_bees*100, no_bees = as.integer(no_bees)) %>% # calculating varroa per 100 bees
  drop_na(Varroa.p) %>%
  mutate(Varroa.p.s = scale(Varroa.p)) %>% rename(Site = site)

data.july <- data.p %>% filter(month == 'July') %>%
  mutate(Varroa.p = varroa_board/no_bees*100, no_bees = as.integer(no_bees)) %>% # calculating varroa per 100 bees
  drop_na(Varroa.p) %>%
  mutate(Varroa.p.s = scale(Varroa.p)) %>% rename(Site = site)

### landuse
load('radii.RData')
radii2 <- data.frame(Site = radii[['r500']]$Site, Ann.fl500 = radii[['r500']]$Ann.fl, OSR500 = radii[['r500']]$OSR, 
                   SNH750 = radii[['r750']]$SNH, Org.farm2000 = radii[['r2000']]$Org.farm)

land <- radii2 %>% mutate(Org.farm.p = Org.farm2000*100/(pi*2000^2/10000), # calculating the % of each variable in the radius
                         SNH.p = SNH750*100/(pi*750^2/10000),
                         Ann.fl.p = Ann.fl500*100/(pi*500^2/10000),
                         OSR.p = OSR500*100/(pi*500^2/10000)) %>%
  mutate(across(Org.farm.p:OSR.p, ~c(scale(.)), .names = '{col}.s')) %>% mutate(Site = as.factor(Site))


### merge
data.final.august <- left_join(data.august, land, by = c('Site' = 'Site'))
data.final.july <- left_join(data.july, land, by = c('Site' = 'Site'))

#save(data.final.august, data.final.july, file = 'Data_final.RData')
