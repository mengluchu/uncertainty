#===========================================================
#                   Preliminar models
#===========================================================


## Libraries
library(tidyverse)
library(boot)
library(MASS)
library(fitdistrplus)

#setwd("C:/Users/Usuario/Desktop/Projects/2021/KAUST/INLA")

#dat_fit = read.csv("dat_fit.csv", header=T)
dat_fit  = read.csv("https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/DENL17_uc.csv")


dim(dat_fit)

head(dat_fit)
glimpse(dat_fit)

options(scipen=999)

# Response variable
dat_fit$mean_value = as.numeric(dat_fit$mean_value)
hist(dat_fit$mean_value)

dat_fit = dat_fit[, c("nightlight_450", "population_1000",  "population_3000", 
        "road_class_1_5000", "road_class_2_100",  "road_class_3_300", "trop_mean_filt", 
        "road_class_3_3000", "road_class_1_100", "road_class_3_100", "road_class_3_5000",
        "Countrycode", "MeasurementType", "AirQualityStationType", "AirQualityStationArea",
        "mean_value")]

glimpse(dat_fit)
summary(dat_fit$mean_value)


par(mfrow=c(1,1))
hist(dat_fit$mean_value)
descdist(dat_fit$mean_value, boot = 2000)


# Selection of variables
full = glm(mean_value~.,data = dat_fit, family  = Gamma(link = "identity"))
step = stepAIC(full, trace = TRUE)
step$anova

backward  = stepAIC(glm(mean_value~.,data = dat_fit, family = Gamma(link = "identity")),direction="backward")
forward   = stepAIC(glm(mean_value~.,data = dat_fit, family = Gamma(link = "identity")),direction="forward")
both      = stepAIC(glm(mean_value~.,data = dat_fit, family = Gamma(link = "identity")),direction="both")


backward$anova
forward$anova
both$anova

# Compare deviances
step$deviance
backward$deviance
forward$deviance
both$deviance   

# Compare AIC
step$aic
backward$aic
forward$aic
both$aic

# Formula for the best model (all the models except 'forward')
formula(step)
formula(backward)
formula(both)

#===============================================================================================================
#                                         Comparative Models
#===============================================================================================================


# First case: The same variables that appear in the document. 

# Second case: Same variables that appear in the document + 
#              Countrycode + 
#              MeasurementType + 
#              AirQualityStationType + 
#              AirQualityStationArea


## Gamma GLM
m1 = glm(mean_value ~ nightlight_450 + population_3000 + road_class_1_5000 + 
         road_class_2_100 + road_class_3_300 + trop_mean_filt + road_class_3_3000 + 
         road_class_1_100 + road_class_3_100 + road_class_3_5000, 
         family = Gamma(link = "identity"), data = dat_fit, na.action = na.omit)
m1$aic


m2 = glm(mean_value ~ nightlight_450 + population_3000 + road_class_1_5000 + 
         road_class_2_100 + road_class_3_300 + trop_mean_filt + road_class_3_3000 + 
         road_class_1_100 + road_class_3_100 + road_class_3_5000 +
         Countrycode + MeasurementType + AirQualityStationType + AirQualityStationArea, 
         family = Gamma(link = "identity"), data = dat_fit, na.action = na.omit)

m2$aic



#=======================================================================================
#                           Lognormal case
#=======================================================================================
m3 = glm(mean_value ~ nightlight_450 + population_3000 + road_class_1_5000 + 
         road_class_2_100 + road_class_3_300 + trop_mean_filt + road_class_3_3000 + 
         road_class_1_100 + road_class_3_100 + road_class_3_5000, 
         family = gaussian(link = "log"), data = dat_fit, na.action = na.omit)
m3$aic


m4 = glm(mean_value ~ nightlight_450 + population_3000 + road_class_1_5000 + 
         road_class_2_100 + road_class_3_300 + trop_mean_filt + road_class_3_3000 + 
         road_class_1_100 + road_class_3_100 + road_class_3_5000 +
         Countrycode + MeasurementType + AirQualityStationType + AirQualityStationArea, 
         family = gaussian(link = "log"), data = dat_fit, na.action = na.omit)

m4$aic




#=======================================================================================
#                                  Gaussian
#=======================================================================================
m5 = glm(mean_value ~ nightlight_450 + population_3000 + road_class_1_5000 + 
           road_class_2_100 + road_class_3_300 + trop_mean_filt + road_class_3_3000 + 
           road_class_1_100 + road_class_3_100 + road_class_3_5000, 
         family = gaussian(link = "identity"), data = dat_fit, na.action = na.omit)
m5$aic


m6 = glm(mean_value ~ nightlight_450 + population_3000 + road_class_1_5000 + 
           road_class_2_100 + road_class_3_300 + trop_mean_filt + road_class_3_3000 + 
           road_class_1_100 + road_class_3_100 + road_class_3_5000 +
           Countrycode + MeasurementType + AirQualityStationType + AirQualityStationArea, 
         family = gaussian(link = "identity"), data = dat_fit, na.action = na.omit)

m6$aic






# Compare AIC
AIC(m1, m2, m3, m4, m5, m6) # Best models considering AIC are 'm1' and 'm2' (gamma distribution)
anova(m1, m2, m3, m4, m5, m6)

# Comparative graph for 'm1' and 'm2'

termplot(m1)
termplot(m2)


# loglikelihood comparison
library(lmtest)
lrtest(m1, m2, m3 , m4, m5, m6)  



#===========================================
#         Mejor predictor lineal
#===========================================
par(mfrow=c(3,2))
scatter.smooth(1:482, rstandard(m1, type='deviance'), col='gray', main = "m1 (model without factors)")
scatter.smooth(1:482, rstandard(m2, type='deviance'), col='gray', main = "m2 (model with factors)")
scatter.smooth(1:482, rstandard(m3, type='deviance'), col='gray', main = "m3 (model without factors)")
scatter.smooth(1:482, rstandard(m4, type='deviance'), col='gray', main = "m4 (model with factors)")
scatter.smooth(1:482, rstandard(m5, type='deviance'), col='gray', main = "m5 (model without factors)")
scatter.smooth(1:482, rstandard(m6, type='deviance'), col='gray', main = "m6 (model with factors)")



# 
# 
par(mfrow=c(3,2))
qqnorm(statmod::qresid(m1)); qqline(statmod::qresid(m1))
qqnorm(statmod::qresid(m2)); qqline(statmod::qresid(m2))
qqnorm(statmod::qresid(m3)); qqline(statmod::qresid(m3))
qqnorm(statmod::qresid(m4)); qqline(statmod::qresid(m4))
qqnorm(statmod::qresid(m5)); qqline(statmod::qresid(m5))
qqnorm(statmod::qresid(m6)); qqline(statmod::qresid(m6))




#===========================================================================
#                              Bayesian GLM
#===========================================================================
# INLA constrols
library(INLA)
library(INLAutils)
cres = list(return.marginals.predictor = TRUE, return.marginals.random = TRUE)
cinla <- list(strategy = 'adaptive', int.strategy = 'eb')  # Estrategias de estimaci?n

# Initial values of parameters
ini.zb <- c(1.834)

#=========================================================================
#                     Gamma model (without factors)
#=========================================================================
formula1 = mean_value ~ nightlight_450 + population_3000 + road_class_1_5000 + 
           road_class_2_100 + road_class_3_300 + trop_mean_filt + road_class_3_3000 + 
           road_class_1_100 + road_class_3_100 + road_class_3_5000


inla1 = inla(formula1,
             family = 'Gamma',
             data = dat_fit, 
             control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
             control.predictor=list(compute=TRUE, link = 1),
             control.results = cres, control.inla = cinla,
             control.mode = list(theta = ini.zb, restart = TRUE),
             verbose=TRUE)

summary(inla1)
inla1$waic$waic


#=========================================================================
#                   Gamma model (with factors)
#=========================================================================
formula2 = mean_value ~ nightlight_450 + population_3000 + road_class_1_5000 + 
           road_class_2_100 + road_class_3_300 + trop_mean_filt + road_class_3_3000 + 
           road_class_1_100 + road_class_3_100 + road_class_3_5000 +
           Countrycode + MeasurementType + AirQualityStationType + AirQualityStationArea
  
inla2 = inla(formula2,
            family = 'Gamma',
            data = dat_fit, 
            control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
            control.predictor=list(compute=TRUE, link = 1),
            control.results = cres, control.inla = cinla,
            control.mode = list(theta = ini.zb, restart = TRUE),
            verbose=TRUE)

summary(inla2)
inla2$waic$waic



#=========================================================================
#                   Lognormal model ( without factors)
#=========================================================================

formula3 = mean_value ~ nightlight_450 + population_3000 + road_class_1_5000 + 
           road_class_2_100 + road_class_3_300 + trop_mean_filt + road_class_3_3000 + 
           road_class_1_100 + road_class_3_100 + road_class_3_5000

inla3 = inla(formula3,
             family = 'lognormal',
             data = dat_fit, 
             control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
             control.predictor=list(compute=TRUE, link = 1),
             control.results = cres, control.inla = cinla,
             control.mode = list(theta = ini.zb, restart = TRUE),
             verbose=TRUE)

summary(inla3)





#=========================================================================
#                     Lognormal model ( with factors)
#=========================================================================
formula4 = mean_value ~ nightlight_450 + population_3000 + road_class_1_5000 + 
           road_class_2_100 + road_class_3_300 + trop_mean_filt + road_class_3_3000 + 
           road_class_1_100 + road_class_3_100 + road_class_3_5000 +
           Countrycode + MeasurementType + AirQualityStationType + AirQualityStationArea

inla4 = inla(formula4,
             family = 'lognormal',
             data = dat_fit, 
             control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
             control.predictor=list(compute=TRUE, link = 1),
             control.results = cres, control.inla = cinla,
             control.mode = list(theta = ini.zb, restart = TRUE),
             verbose=TRUE)


#==========================
# Summary
#==========================

summary(inla1)
summary(inla2)
summary(inla3) 
summary(inla4)


# WAIC
inla1$waic$waic
inla2$waic$waic
inla3$waic$waic
inla4$waic$waic

inla1$dic$dic
inla2$dic$dic
inla3$dic$dic
inla4$dic$dic


slcpo <- function(m, na.rm = TRUE) {
  - sum(log(m$cpo$cpo), na.rm = na.rm)
}


c("m1" = slcpo(inla1), 
  "m2" = slcpo(inla2), 
  "m3" = slcpo(inla3), 
  "m4" = slcpo(inla4))




library(brinla)
bri.fixed.plot(inla2)


b0 =  inla2$marginals.fixed[[1]]
plot_b0 = ggplot(data.frame(inla.smarginal(b0)), aes(x, y)) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Posterior distribution fo the intercept", y = "Intercept", x = "") + 
  theme_bw()

plot_b0

# compare models with cross validation (lowest is best)
cval <- function(result=inla1){
  test <- sum(result$cpo$failure)
  if(test != 0) print("some cpo estimates failed the quality test. do not use.")
  return(-2 * sum(log(result$cpo$cpo), na.rm = T))
}
cval(inla1); cval(inla2); cval(inla3); cval(inla4)


# look at model fit to data in both cases
par(mfrow=c(2,2))
hist(inla1$cpo$pit, main = "first case with INLA (gamma)") # should be uniformish
plot(inla1$summary.fitted.values$mean, dat_fit$mean_value, main = "first case with INLA (gamma)"); abline(0,1)
hist(inla2$cpo$pit, main = "second case with INLA (gamma)") # should be uniformish
plot(inla2$summary.fitted.values$mean, dat_fit$mean_value, main = "second case with INLA (gamma)"); abline(0,1)


