# Code for the re-evaluation in JAE

library(readr)
library(dplyr)
library(tidyverse)
library(lme4)
library(car)
library(MuMIn)
library(DHARMa)
options(na.action='na.fail')


load('data_allmonths.RData')

data.cut <- data.merged %>% filter(round %in% c(3,4))

par <- read.csv('HBULK.csv')

data.p <- data.cut  %>% inner_join(par, by = c('month' = 'date', 'sample' = 'sample', 'site' = 'site'))

data.august <- data.p %>% filter(month == 'August') %>%
  mutate(Varroa.p = varroa_board/no_bees*100, no_bees = as.integer(no_bees)) %>% # calculating varroa per 100 bees
  drop_na(Varroa.p) %>%
  mutate(Varroa.p.s = scale(Varroa.p))

##### Reviewer 1 JAE: checking if the results hold if we use different radii
#### In the initial submission we use pearson's correlation coefficient of landscape variables at different radii and abs. number of bees
#### to determine a best radius for each landscape variable (details in the ms; 500m for Annual flowers and OSR, 750m for SNH and 2000m for Organic farming).
#### The reviewer doubts regarded if the results will hold if we use different radii. We will reanalyse the data for the three radii (500, 750, 20000):

load('radii.RData')

data.allradii <- lapply(radii, function(x) merge(x, data.august, by.x = 'Site', by.y ='site'))

d500 <- data.allradii[[1]] %>% mutate(Org.farm.p = Org.farm*100/(pi*500^2/10000), # calculating the % of each variable in 500m radius
                                      SNH.p = SNH*100/(pi*500^2/10000),
                                      Ann.fl.p = Ann.fl*100/(pi*500^2/10000),
                                      OSR.p = OSR*100/(pi*500^2/10000)) %>%
  mutate(across(Org.farm.p:OSR.p, ~c(scale(.)), .names = '{col}.s')) %>% mutate(Site = as.factor(Site))

d750 <- data.allradii[[2]] %>% mutate(Org.farm.p = Org.farm*100/(pi*750^2/10000), # calculating the % of each variable in 750m radius
                                      SNH.p = SNH*100/(pi*750^2/10000),
                                      Ann.fl.p = Ann.fl*100/(pi*750^2/10000),
                                      OSR.p = OSR*100/(pi*2000^2/10000)) %>%
  mutate(across(Org.farm.p:OSR.p, ~c(scale(.)), .names = '{col}.s')) %>% mutate(Site = as.factor(Site))

d2000 <- data.allradii[[7]] %>% mutate(Org.farm.p = Org.farm*100/(pi*500^2/10000), # calculating the % of each variable in 2km radius
                                       SNH.p = SNH*100/(pi*2000^2/10000),
                                       Ann.fl.p = Ann.fl*100/(pi*500^2/10000),
                                       OSR.p = OSR*100/(pi*2000^2/10000)) %>%
  mutate(across(Org.farm.p:OSR.p, ~c(scale(.)), .names = '{col}.s')) %>% mutate(Site = as.factor(Site))



radiuslist <- list(d500 = d500, d750 = d750,d2000 = d2000) #all data frames in a list


##### PARASITE RICHNESS MODELS COMPARISON ######

# the three radii models + dredge:
models.selection <- lapply(radiuslist, function(x) dredge(lmer(log(p_rich) ~ Varroa.p.s + Org.farm.p.s + Ann.fl.p.s + SNH.p.s + OSR.p.s + (1|Site), REML = F, data = x), m.lim = c(0,5), rank = 'AICc'))

#### SAME RESULT OF MODEL SELECTION AMONG THE THREE RADII
models.selection[['d500']] %>% filter(delta < 2) # delta < 2, because within delta = 2 models considered equally good (predictive-wise)
models.selection[['d750']] %>% filter(delta < 2)
models.selection[['d2000']] %>% filter(delta < 2)

best.models.p <- lapply(radiuslist, function(x) lmer(log(p_rich) ~ Varroa.p.s + Org.farm.p.s + (1|Site), REML = T, data = x))
lapply(best.models.p, function(x) coef(summary(x))[ , c("Estimate", 'Std. Error')]) # estimates are the same
lapply(best.models.p, summary) # full summary
lapply(best.models.p, function(x) Anova(x, type = 3))

## visualisation of coefs
best.models.p %>% lapply(function(x) coef(summary(x))[ , "Estimate"]) %>% bind_rows(.id = 'id') %>% 
  pivot_longer(cols = where(is.numeric), names_to = 'Predictor', values_to = 'Estimate') %>% 
  ggplot(aes(Predictor, Estimate, col = factor(id))) +
  geom_point(size = 3, shape = 19, position = position_dodge(width = 0.2))+
  theme_minimal(base_size = 16)+
  geom_hline(yintercept=0)+
  scale_color_manual(values = c('red', 'blue', 'green'))+
  labs(title = 'Parasite richness model estimates')
#ggsave('coefs_parasiterichness.pdf', device = 'pdf')

### model assumptions validation:
### change from [[1]] to [[3]] to see different radii
testDispersion(best.models.p[[1]])
plot(so <- simulateResiduals(fittedModel = best.models.p[[1]], plot = F))

##### VARROA MODELS COMPARISON ######

# the three radii models + dredge:
models.selection <- lapply(radiuslist, function(x) dredge(lmer(Varroa.p ~ Org.farm.p.s + Ann.fl.p.s + SNH.p.s + OSR.p.s + (1|Site), family = Gamma(link='log'), data = x), m.lim = c(0,5), rank = 'AICc'))

#### MODEL SELECTION AMONG THE THREE RADII
models.selection[['d500']] %>% filter(delta < 2)
models.selection[['d750']] %>% filter(delta < 2)
models.selection[['d2000']] %>% filter(delta < 2) # model selection also wants to include OSR!

# so first with Annual flowers and SNH
best.models.v <- lapply(radiuslist, function(x) glmer(Varroa.p ~ Ann.fl.p.s + SNH.p.s + (1|Site), family = Gamma(link='log'), data = x))
lapply(best.models.v, function(x) coef(summary(x))[ , c("Estimate", 'Std. Error')]) # estimates are the same
lapply(best.models.v, summary) # full summary
lapply(best.models.v, Anova)

# then with OSR and SNH
best.models.v <- lapply(radiuslist, function(x) glmer(Varroa.p ~ OSR.p.s + SNH.p.s + (1|Site), family = Gamma(link='log'), data = x))
lapply(best.models.v, function(x) coef(summary(x))[ , c("Estimate", 'Std. Error')]) # estimates are the same
lapply(best.models.v, summary) # full summary
lapply(best.models.v, Anova)

# then with Annual flower and OSR
best.models.v <- lapply(radiuslist, function(x) glmer(Varroa.p ~ OSR.p.s + Ann.fl.p.s + (1|Site), family = Gamma(link='log'), data = x))
lapply(best.models.v, function(x) coef(summary(x))[ , c("Estimate", 'Std. Error')]) # estimates are the same
lapply(best.models.v, summary) # full summary
lapply(best.models.v, Anova)

## visualisation of coefs
best.models.v %>% lapply(function(x) coef(summary(x))[ , "Estimate"]) %>% bind_rows(.id = 'id') %>% 
  pivot_longer(cols = where(is.numeric), names_to = 'Predictor', values_to = 'Estimate') %>% 
  ggplot(aes(Predictor, Estimate, col = factor(id))) +
  geom_point(size = 3, shape = 19, position = position_dodge(width = 0.2))+
  theme_minimal(base_size = 16)+
  geom_hline(yintercept=0)+
  scale_color_manual(values = c('red', 'blue', 'green'))+
  labs(title = 'Varroa model estimates')
ggsave('coefs_varroa.pdf', device = 'pdf')

### model assumptions validation:
### change from [[1]] to [[3]] to see different radii
testDispersion(best.models.v[[1]])
plot(so <- simulateResiduals(fittedModel = best.models.v[[1]], plot = F))


##### COLONY GROWTH MODELS COMPARISON ######

# the three radii models + dredge:
models.selection <- lapply(radiuslist, function(x) dredge(glmer.nb(no_bees ~ log(p_rich) + Varroa.p.s + Org.farm.p.s + Ann.fl.p.s + SNH.p.s + OSR.p.s + (1|Site), family='nbinom2', data = x), m.lim = c(0,5), rank = 'AICc'))

#### SAME RESULT OF MODEL SELECTION AMONG THE THREE RADII
models.selection[['d500']] %>% filter(delta < 2) # delta < 2, because within delta = 2 models considered equally good (predictive-wise)
models.selection[['d750']] %>% filter(delta < 2)
models.selection[['d2000']] %>% filter(delta < 2)

best.models.c <- lapply(radiuslist, function(x) glmer.nb(no_bees ~ log(p_rich) + Org.farm.p.s + (1|Site), data = x))
lapply(best.models.c, function(x) coef(summary(x))[ , c("Estimate", 'Std. Error')]) # estimates are the same
lapply(best.models.c, summary) # full summary
lapply(best.models.c, Anova)

## visualisation of coefs
best.models.c %>% lapply(function(x) coef(summary(x))[ , "Estimate"]) %>% bind_rows(.id = 'id') %>% 
  pivot_longer(cols = where(is.numeric), names_to = 'Predictor', values_to = 'Estimate') %>% 
  ggplot(aes(Predictor, Estimate, col = factor(id))) +
  geom_point(size = 3, shape = 19, position = position_dodge(width = 0.2))+
  theme_minimal(base_size = 16)+
  geom_hline(yintercept=0)+
  scale_color_manual(values = c('red', 'blue', 'green'))+
  labs(title = 'Colony growth model estimates')
ggsave('coefs_colonygrowth.pdf', device = 'pdf')

### model assumptions validation:
### change from [[1]] to [[3]] to see different radii
testDispersion(best.models.c[[1]])
plot(so <- simulateResiduals(fittedModel = best.models.c[[1]], plot = F))

###### SEM COMPARISON ######
library(piecewiseSEM)

# 500m radius
m1 <- glmer(Varroa.p ~ Ann.fl.p.s + SNH.p.s + (1|Site), family = Gamma(link='log'), d500) 
m2 <- lmer(log(p_rich) ~ Varroa.p + Org.farm.p.s + (1|Site), d500) 
m3 <- glmer.nb(no_bees ~ log(p_rich) + SNH.p.s + Org.farm.p.s + (1|Site), family='nbinom2', d500) 
summary(psem(m1, m2, m3))


# 750 radius
m1 <- glmer(Varroa.p ~ Ann.fl.p.s + SNH.p.s + (1|Site), family = Gamma(link='log'), d750) 
m2 <- lmer(log(p_rich) ~ Varroa.p + Org.farm.p.s + (1|Site), d750) 
m3 <- glmer.nb(no_bees ~ log(p_rich) + SNH.p.s + Org.farm.p.s + (1|Site), family='nbinom2', d750) 
summary(psem(m1, m2, m3))

# 750 radius
m1 <- glmer(Varroa.p ~ Ann.fl.p.s + SNH.p.s + (1|Site), family = Gamma(link='log'), d2000) 
m2 <- lmer(log(p_rich) ~ Varroa.p + Org.farm.p.s + (1|Site), d2000) 
m3 <- glmer.nb(no_bees ~ log(p_rich) + SNH.p.s + Org.farm.p.s + (1|Site), family='nbinom2', d2000) 
summary(psem(m1, m2, m3))


### trying varroa log
#500m
m1 <- lmer(log(Varroa.p) ~ Ann.fl.p.s + SNH.p.s + (1|Site), d500) 
m2 <- lmer(log(p_rich) ~ log(Varroa.p) + Org.farm.p.s + (1|Site), d500) 
m3 <- glmer.nb(no_bees ~ log(p_rich) + SNH.p.s + Org.farm.p.s + (1|Site), family='nbinom2', d500) 
summary(psem(m1, m2, m3))

#750
m1 <- lmer(log(Varroa.p) ~ Ann.fl.p.s + SNH.p.s + (1|Site), d750) 
m2 <- lmer(log(p_rich) ~ log(Varroa.p) + Org.farm.p.s + (1|Site), d750) 
m3 <- glmer.nb(no_bees ~ log(p_rich) + SNH.p.s + Org.farm.p.s + (1|Site), family='nbinom2', d750) 
summary(psem(m1, m2, m3))

#2000
m1 <- lmer(log(Varroa.p) ~ Ann.fl.p.s + SNH.p.s + (1|Site), d2000) 
m2 <- lmer(log(p_rich) ~ log(Varroa.p) + Org.farm.p.s + (1|Site), d2000) 
m3 <- glmer.nb(no_bees ~ log(p_rich) + SNH.p.s + Org.farm.p.s + (1|Site), family='nbinom2', d2000) 
summary(psem(m1, m2, m3))
