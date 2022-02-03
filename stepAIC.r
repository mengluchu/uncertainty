#===========================================================
#                   Feature selection
#===========================================================
## Libraries
library(tidyverse)
library(boot)
library(MASS)

setwd("C:/Users/Usuario/Desktop/Projects/2021/KAUST/INLA")

dat_fit = read.csv("dat2.csv", header=T)
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
                      "urbantype", "mean_value")]

glimpse(dat_fit)
summary(dat_fit$mean_value)


dat_fit$Countrycode  = as.factor(dat_fit$Countrycode)
dat_fit$MeasurementType  = as.factor(dat_fit$MeasurementType)
dat_fit$AirQualityStationType = as.factor(dat_fit$AirQualityStationType)
dat_fit$AirQualityStationArea = as.factor(dat_fit$AirQualityStationArea)
dat_fit$urbantype = as.factor(dat_fit$urbantype)


# Selection of variables
full = glm(sqrt(mean_value)~.,data = dat_fit, family  = gaussian(link = "identity"))
step = stepAIC(full, trace = TRUE)
step$anova

backward  = stepAIC(glm(sqrt(mean_value)~.,data = dat_fit, family = gaussian(link = "identity")),direction="backward")
forward   = stepAIC(glm(sqrt(mean_value)~.,data = dat_fit, family = gaussian(link = "identity")),direction="forward")
both      = stepAIC(glm(sqrt(mean_value)~.,data = dat_fit, family = gaussian(link = "identity")),direction="both")


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


# The best set of covariates is:

  nightlight_450 + population_1000 + population_3000 + 
  road_class_1_5000 + road_class_2_100 + road_class_3_300 + 
  trop_mean_filt + road_class_1_100 + Countrycode + MeasurementType + 
  AirQualityStationType + AirQualityStationArea + urbantype

