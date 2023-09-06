library(readr)
library(dplyr)
library(tidyverse)
library(lme4)
library(car)
library(MuMIn)
library(DHARMa)
library(effects)
options(na.action='na.fail')


load(file = 'Data/Data_final.RData')
data.final.august$Site <- as.factor(data.final.august$Site)
load(file = 'Data/Flower_cover.RData')
data.final.august <- data.final.august %>% left_join(fl.cv[1:2], by = c('Site' = 'Site')) %>% mutate(m.cover.p.s = scale(m.cover.p)[,1])


###################################
### 1. PARASite  RICHNESS MODEL ###
###################################


# trying to fit parasite richness (p_rich) to different distributions
#prich.M <- glmer(cbind(p_rich, 11-p_rich2) ~ Varroa.p.s + Org.farm.p.s + (1|Site), family = 'binomial', data.final.august) ### doesn't work well, assumes all paraSites have equal probability of being contracted - not true
#prich.M <- glmer(p_rich ~ Varroa.p.s + Org.farm.p.s + (1|Site), family = 'poisson', data.final.august) ### not good either - very underdispersed
#prich.M <- glmer(p_rich ~ Varroa.p.s + Org.farm.p.s + (1|Site), family = 'nbinom1', data.final.august) ### not good either - very underdispersed
#library(glmmTMB)
#data$p_rich_r <- data$p_rich/11
#prich.M <- glmmTMB(p_rich_r ~ Varroa.p.s + Org.farm.p.s + (1|Site), family = beta_family(), data) ### works well but cannot be used in piecewiseSEM


#using logarithm of parasite richness


prich.M<- lmer(log(p_rich+1) ~ Varroa.p.s + (↕ + Ann.fl.p.s + SNH.p.s + OSR.p.s)^2 + (1|Site), REML=F, data = data.final.august)

# model selection
b<-dredge(prich.M, m.lim = c(0, 5), rank="AICc") 
subset(b, delta < 2)

# best model
prich.M<- lmer(log(p_rich+1) ~ Varroa.p.s + Org.farm.p.s + (1|Site), REML=T, data = data.final.august)

# checking if the results hold w/o the influential point (Nor918V1)
data.out <- data.final.august %>% filter(sample != 'Nor918V1') 
prich.M <- lmer(log(p_rich+1) ~ Varroa.p.s + Org.farm.p.s + (1|Site), REML=T, data = data.out)

vif(prich.M)
summary(prich.M)
car::Anova(prich.M, type=3)
plot(allEffects(prich.M))
testDispersion(prich.M)
plot(so <- simulateResiduals(fittedModel = prich.M, plot = F))
r.squaredGLMM(prich.M)


#######################
### 2. VARROA MODEL ###
#######################

#trying to fit Varroa
#data$Varroa_L <- log(data$Varroa.p)
#var.M <- lmer(Varroa_L ~ SNH.p.s + Ann.fl.p.s + (1|Site), REML=T, data.final.august) ### works alright, but recommended only as a last resort
#var.M <- glmer(Varroa.p ~ SNH.p.s + Ann.fl.p.s + (1|Site), family = Gamma(link='log'), data.final.august) #same result as logarithm, but a bit underdispersed (though not significantly)
#var.M <- glmmTMB(Varroa.p ~ SNH.p.s + Ann.fl.p.s + (1|Site), family = Gamma(link='log'), REML=T, data.final.august) #same result, better dispersion but cannot be used in piecewiseSEM
#var.M <- glmer(varroa_board ~  Ann.fl.p.s + SNH.p.s + (1|Site), family = 'poisson', weights=no_bees, data.final.august) #poisson with weights as no_bees (each bee had an equal prob. of aquiring a varroa mite)


var.M <- glmer(Varroa.p ~ (Org.farm.p.s + Ann.fl.p.s + SNH.p.s + OSR.p.s)^2 + (1|Site), family = Gamma(link='log'), data.final.august)
vif(var.M)
# model selection
b<-dredge(var.M, m.lim = c(0, 4), rank="AICc") 
subset(b, delta < 2)

# best model with interaction, but VIF is inflated - removed interaction
var.M <- glmer(Varroa.p ~ Ann.fl.p.s + SNH.p.s + (1|Site), family = Gamma(link='log'), data.final.august)

vif(var.M)
summary(var.M)
car::Anova(var.M, type=3)
plot(allEffects(var.M))
testDispersion(var.M)
plot(so <- simulateResiduals(fittedModel = var.M, plot=T))
r.squaredGLMM(var.M)


##############################
### 3. COLONY GROWTH MODEL ###
##############################

col.M <- glmer.nb(no_bees ~ log(p_rich) + Varroa.p.s + (Org.farm.p.s + Ann.fl.p.s + SNH.p.s + OSR.p.s)^2 + (1|Site), family='nbinom2', data.final.august) 

# model selection
b<-dredge(col.M, m.lim = c(0, 3), rank="AICc") 
subset(b, delta < 2)

# best model 
col.M <- glmer.nb(no_bees ~ log(p_rich) + Org.farm.p.s + (1|Site), family='nbinom2', data.final.august) # second best because of residuals

vif(col.M)
summary(col.M)

car::Anova(col.M, type=3)
plot(allEffects(col.M))
testDispersion(col.M)
plot(so <- simulateResiduals(fittedModel = col.M, plot = F))
r.squaredGLMM(col.M)



##### individual parasites instead of parasite richness, model 3b:

#### NOT WORKING with BQCV - too many positives
col.M <- glmer.nb(no_bees ~ dwvb + tryp + bqcv + Varroa.p.s + Org.farm.p.s + Ann.fl.p.s + SNH.p.s + OSR.p.s + (1|Site), data.final.august) 

#### without BQCV
col.M <- glmer.nb(no_bees ~ dwvb + tryp + Varroa.p.s + Org.farm.p.s + Ann.fl.p.s + SNH.p.s + OSR.p.s + (1|Site), data.final.august) 
col.M <- glmer.nb(no_bees ~ bqcv + (1|Site), data.final.august) 

b<-dredge(col.M, m.lim = c(0, 3), rank="AICc") 
subset(b, delta < 2)

col.M <- glmer.nb(no_bees ~ tryp  + Org.farm.p.s + SNH.p.s +  (1|Site), control = glmerControl(optimizer = 'bobyqa'), family='nbinom2',data.final.august)
vif(col.M)
summary(col.M)

car::Anova(col.M, type=3)
plot(allEffects(col.M))
testDispersion(col.M)
plot(so <- simulateResiduals(fittedModel = col.M, plot = F))
r.squaredGLMM(col.M)


#################################
### 4. COLONY SURVIVAL MODEL  ###
#################################

survival <- read.csv('wintermortality2021_binom.csv')

data_S <- data.final.august %>% left_join(survival, by=c('sample' = 'sample')) %>%
                                            filter(Site != 'Nor1070') #merging with survival - 2 colonies less, because they had to been removed from the Site in the fall
data_S$no_bees.s <- scale(data_S$no_bees)

# not going through as one model - splitting into two
surv.M <- glmer(status ~ Org.farm.p.s + Ann.fl.p.s + SNH.p.s + OSR.p.s + (1|Site), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), family=binomial, data_S)
surv.M <- glmer(status ~ log(p_rich+1) + Varroa.p.s + no_bees.s + (1|Site), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), family=binomial, data_S)

# model selection - first the landscape one, then the parasite one
b<-dredge(surv.M, m.lim = c(0, 3), rank="AICc") 
subset(b, delta < 2)

# best landscape model = NULL
# best parasite model
surv.M <- glmer(status ~  no_bees.s + (1|Site), family=binomial, data_S)

summary(surv.M)
car::Anova(surv.M, type=3)
plot(allEffects(surv.M))
testDispersion(surv.M)
plot(so <- simulateResiduals(fittedModel = surv.M, plot = F))
r.squaredGLMM(surv.M)



###############
### 5. SEM  ### 
###############

library(piecewiseSEM)

# best models - skipped colony survival
m1 <- glmer(Varroa.p ~ Ann.fl.p.s + SNH.p.s + (1|Site), family = Gamma(link='log'), data.final.august) 
m2 <- lmer(log(p_rich) ~ Varroa.p + Org.farm.p.s + (1|Site), data.final.august) 
m3 <- glmer.nb(no_bees ~ log(p_rich) + SNH.p.s + Org.farm.p.s + (1|Site), family='nbinom2', data.final.august) 

sem <- psem(m1, m2, m3)

summary(sem)

plot(sem) #not working with missing std. estimates



########### Standarization of the estimates for Gamma  distribution 
########### (not implemented in piecewiseSEM) 


## m1
R2 <- cor(data.final.august$Varroa.p, predict(m1, type = "response"))^2 # non-linear predictions

sd.yhat <- sqrt(var(predict(m1, type = "link"))/R2)

summary(m1)$coefficients[2, 1] * sd(data.final.august$Ann.fl.p.s)/sd.yhat #-0.35
summary(m1)$coefficients[3, 1] * sd(data.final.august$SNH.p.s)/sd.yhat #0.27



### m3 - nibinom, checking if it's working

R2 <- cor(data$no_bees, predict(m3, type = "response"))^2 # non-linear predictions

sd.yhat <- sqrt(var(predict(m3, type = "link"))/R2)

summary(m3)$coefficients[4, 1] * sd(data.final.august$Org.farm.p.s)/sd.yhat 
summary(m3)$coefficients[3, 1] * sd(data.final.august$SNH.p.s)/sd.yhat 
summary(m3)$coefficients[2, 1] * sd(data.final.august$log(p_rich))/sd.yhat 


######### total effects in SEM ##########

# org farm
-0.33 * -0.4 + 0.38 

# snh
0.27 * 0.58 * -0.4

#annfl
-0.35 * 0.58 * -0.4

#varroa
0.58 * -0.4


#############################################
### Recalculating coefficients from GLMMs ###
#############################################


## for Gamma and negative binomial (log link) - exponentiating:
#exp(X)

## colony growth model
## Org_crop = 1.17
## SNH = 0.9
## p_rich = 0.82

## for binomial models (logit link) - reversing log-odds:
#exp(X)/(1+exp(X))


#############################################
### INDIVIDUAL PARASITE PREVALENCE MODELS ###
#############################################

data.both.wide <- inner_join(data.final.august, data.final.july[1:23], by = join_by('sample' == 'sample', 'Site' == 'Site'), suffix = c('2', '1'))

# DWV-A <10% positives

# DWV-B
dwvb.M <- glmer(dwvb2 ~ dwvb1 + Varroa.p.s2 + (Org.farm.p.s + Ann.fl.p.s + SNH.p.s + OSR.p.s)^2 + (1|Site), family='binomial', data=data.both.wide) #excluding 2-way interactions becasue of VIF

b<-dredge(dwvb.M, m.lim = c(0, 4), rank="AICc") 
subset(b, delta < 2)

dwvb.M <- glmer(dwvb2 ~  (1 | Site), family='binomial', data=data.both.wide)
vif(dwvb.M)
summary(dwvb.M)
plot(allEffects(dwvb.M))
testDispersion(dwvb.M)
plot(so <- simulateResiduals(fittedModel = dwvb.M, plot = F))

#BQCV - no interaction term because of the error
bqcv.M <- glmer(bqcv2 ~ bqcv1 + Org.farm.p.s + Ann.fl.p.s + SNH.p.s + OSR.p.s + (1 | Site), family='binomial', control=glmerControl(nAGQ0initStep=F, tolPwrss=1e-12), data=data.both.wide)

b<-dredge(bqcv.M, m.lim = c(0, 4), rank="AICc") 
subset(b, delta < 2)

bqcv.M <- glmer(bqcv2 ~  Org.farm.p.s + SNH.p.s + OSR.p.s + (1 | Site), family='binomial', data=data.both.wide) #third best model because of convergence problems
summary(bqcv.M)
plot(allEffects(bqcv.M))
testDispersion(bqcv.M)
plot(so <- simulateResiduals(fittedModel = bqcv.M, plot = F))


# SBV <10% positives
# ABPV <10% positives
# CBPV <10% positives
# SBPV <10% positives

# Trypanosomes
tryp.M <- glmer(tryp2 ~ tryp1 + (Org.farm.p.s + Ann.fl.p.s + SNH.p.s + OSR.p.s)^2 + (1 | Site), family='binomial', control=glmerControl(nAGQ0initStep=F, tolPwrss=1e-12), data=data.both.wide)

b<-dredge(tryp.M, m.lim = c(0, 5), rank="AICc") 
subset(b, delta < 2)

tryp.M <- glmer(tryp2 ~ tryp1  + (1 | Site), family='binomial', data=data.both.wide)
summary(tryp.M)
Anova(tryp.M, type = 3)
plot(allEffects(tryp.M))
testDispersion(tryp.M)
plot(so <- simulateResiduals(fittedModel = tryp.M, plot = F))
r.squaredGLMM(tryp.M)

# Neogregarines <10% positives
# Nosema bombi <10% positives
# Nosema ceranae <10% positives


############################################
### MORAN'S I - SPATIAL AUTOCORRELATION  ###
############################################


library(spdep)
coord <- read.csv("Data/coordinates.csv")
coord <- coord %>% distinct()

# paraSite  richness
plot(so <- simulateResiduals(fittedModel = m2, plot = F))
so2 <- recalculateResiduals(so, group = data$Site) #because of 2 data points from each location
testSpatialAutocorrelation(so2, coord$long, coord$lat, plot = T)
# colony growth
plot(so <- simulateResiduals(fittedModel = m3, plot = F))
so2 <- recalculateResiduals(so, group = data$Site) #because of 2 data points from each location
testSpatialAutocorrelation(so2, coord$long, coord$lat, plot = T)
# varroa infestation
plot(so <- simulateResiduals(fittedModel = m1, plot = F))
so2 <- recalculateResiduals(so, group = data$Site) #because of 2 data points from each location
testSpatialAutocorrelation(so2, coord$long, coord$lat, plot = T)
# colony survival
plot(so <- simulateResiduals(fittedModel = surv.M, plot = F))
so2 <- recalculateResiduals(so, group = data_S$Site) #because of 2 data points from each location
testSpatialAutocorrelation(so2, coord$long, coord$lat, plot = T)