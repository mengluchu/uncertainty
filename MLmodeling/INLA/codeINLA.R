library(INLA)
library(APMtools)
library(dplyr)
library(dismo)
source("https://raw.githubusercontent.com/mengluchu/uncertainty/master/MLmodeling//INLA/INLA_util.R")
#################################################
d= read.csv("https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/DENL17_uc.csv")
# Covariates
covnames0 <- c("nightlight_450", "population_1000", "population_3000",
               "road_class_1_5000", "road_class_2_100", "road_class_3_300", "trop_mean_filt",
               "road_class_3_3000", "road_class_1_100", "road_class_3_100"
               , "road_class_3_5000", "road_class_1_300", "road_class_1_500",
                "road_class_2_1000", "nightlight_3150", "road_class_2_300", "road_class_3_1000", "temperature_2m_7")

d <- d[, c("mean_value", "Longitude", "Latitude", covnames0)]

# covnames0 <- NULL
covnames <- c("b0", covnames0)  # covnames is intercept and covnames0

# Data for estimation. Create variables y with the response, coox and cooy with the coordinates, and b0 with the intercept (vector of 1s)
d$y <- d$mean_value # response
d$coox <- d$Longitude
d$cooy <- d$Latitude
d$b0 <- 1 # intercept
d$real <- d$y

# Cross-validation
VLA <- lapply(1:20, FUN = INLA_crossvali, d = d, covnames = covnames)
(VLA <- data.frame(LA = rowMeans(data.frame(VLA))))

#############################################################################
### INLA example, using all data, i.e. training the same as testing ####
#############################################################################
dp <- d
# Call inla()
lres <- fnFitModelINLA(d, dp = NULL, covnames, TFPOSTERIORSAMPLES = FALSE, formulanew = NULL)
res <- lres[[1]]
stk.full <- lres[[2]]
mesh <- lres[[3]]
# Get predictions
dres <- fnGetPredictions(res, stk.full, mesh, d, dp, covnames, NUMPOSTSAMPLES = -1, cutoff_exceedanceprob = 30)
# Goodness of fit
APMtools::error_matrix(validation = dres$real, prediction = dres$pred_mean)
cor(dres$real, dres$pred_mean)
mean(dres$pred_ll <= dres$real &  dres$real <= dres$pred_ul)


 
# RMSE          7.5002554
# RRMSE         0.3136362
# IQR           7.4369859
# rIQR          0.3392630
# MAE           5.5197591
# rMAE          0.2308669
# rsq           0.6593912
# explained_var 0.6602920
# cor           0.8173491
# covprob       0.5226804

# with less variables
#RMSE          7.2588620
#RRMSE         0.3035446
#IQR           7.4579043
#MAE           5.3791243
##rIQR          0.3405347
#rMAE          0.2249874
#rsq           0.6807286
#explained_var 0.6818688
#cor           0.8296228
#covprob       0.4737113

# with stacked-learners using INLA
#RMSE          7.1359408
#RRMSE         0.2984085
#IQR           7.1528871
#rIQR          0.3269740
#MAE           5.2621438
#rMAE          0.2200768
#rsq           0.6902379
#explained_var 0.6910040
#cor           0.8343859
#covprob95     0.3582474
#covprob90     0.3144330
#covprob50     0.1484536