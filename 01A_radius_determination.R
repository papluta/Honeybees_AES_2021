library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

# list of radii in 250 m increments
load('Data/Radii_raw.RData')

# calculating the 4 landscape variables for each radius
land_fun <- function(x){ x %>% 
    mutate(SNH = semi_natur + Fallow + Other_AUM + Grassy_str + Flower_fieBS2 + Flower_fieBS12/2 + Flower_fie,
           Ann.fl = Flower_fieBS11 + Flower_fieBS12/2) %>%
    rename(Org.farm = CropBV1) %>% 
    dplyr::select(Site, Ann.fl, SNH, Org.farm, OSR)
}

radii <- radii.raw %>% lapply(land_fun)

#save(radii, file = 'radii.RData')

# loading the abs. number of bees
load('Data/Data_final.RData')
par <- data.final.august %>% dplyr::select(Site, no_bees)

# putting all radii of each variable together

AF <- radii %>% lapply(ungroup) %>% lapply(function(x) x %>% dplyr::select(Ann.fl)) %>% bind_cols(.name_repair = 'universal') %>% mutate(Site = radii[[1]]$Site) %>% right_join(par, by = join_by(Site)) 

Org_crop <- radii %>% lapply(ungroup) %>% lapply(function(x) x %>% dplyr::select(Org.farm)) %>% bind_cols(.name_repair = 'universal') %>% mutate(Site = radii[[1]]$Site) %>% right_join(par, by = join_by(Site))

SNH <- radii %>% lapply(ungroup) %>% lapply(function(x) x %>% dplyr::select(SNH)) %>% bind_cols(.name_repair = 'universal') %>% mutate(Site = radii[[1]]$Site) %>% right_join(par, by = join_by(Site))

OSR <- radii %>% lapply(ungroup) %>% lapply(function(x) x %>% dplyr::select(OSR)) %>% bind_cols(.name_repair = 'universal') %>% mutate(Site = radii[[1]]$Site) %>% right_join(par, by = join_by(Site))


#############  PLOTS OF PEARSON'S CORRELATIONS  ##############

par(mfrow=c(2,2))
#Annual_fl
plot(cor(AF[,c(1:7)], AF$no_bees, use="complete"), ylab = 'Annual flower')+
  abline(h=0)


#Org_crop
plot(cor(Org_crop[,c(1:7)], Org_crop$no_bees, use="complete"), ylab = 'Organic crop')+
  abline(h=0)

#semi
plot(cor(SNH[,c(1:7)], SNH$no_bees, use="complete"), ylab = 'SNH')+
  abline(h=0)

#OSR
plot(cor(OSR[,c(1:7)], OSR$no_bees, use="complete"), ylab = 'OSR')+
  abline(h=0)

