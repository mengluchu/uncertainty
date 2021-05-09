#devtools::install_github("mengluchu/APMtools") 
library(INLA)
library(APMtools)
library(dplyr)
library(h2o)
library(dismo)
#h2o.shutdown(prompt = F)
h2o.init(nthreads=-1)
source("https://raw.githubusercontent.com/mengluchu/uncertainty/master/MLmodeling//INLA/INLA_util.R")

############## run for 100 m resolution
resolution = 100 

d <- read.csv("https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/DENL17_uc.csv")

if (resolution ==100)
{
  d = d%>%dplyr::select(-c(industry_25,industry_50,road_class_1_25,road_class_1_50,road_class_2_25,road_class_2_50,   road_class_3_25,road_class_3_50))
} 

# Data for estimation. Create variables y with the response, coox and cooy with the coordinates, and b0 with the intercept (vector of 1s)
d$y <- d$mean_value # response
d$coox <- d$Longitude
d$cooy <- d$Latitude
d$b0 <- 1 # intercept
d$real <- d$y
dp <- d
 
# stacked model with INLA

covnames <- c("lasso", "rf", "xgb")
formula <- as.formula("y ~ 0 + f(lasso, model = 'clinear', range = c(0, 1), initial = 1/3) +
                         f(rf, model = 'clinear', range = c(0, 1), initial = 1/3) +
                         f(xgb, model = 'clinear', range = c(0, 1), initial = 1/3) + f(s, model = spde)") 
 
# Cross-validation
VLA <- lapply(1:20, FUN = INLA_stack_crossvali, d = d, formula = formula, covnames = covnames, family="gaussian")
(VLA <- data.frame(LA = rowMeans(data.frame(VLA))))

# stacked model without INLA
ensembleh2o <- lapply(1:2, FUN =ensemble, d = d )

save(ensembleh2o, file = "ensembleh2o.rda")
save(VLA, file = "VLA.rda")

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

'almost  the same if i use 3000 or 300 xgboos trees! lower learning rate + less trees = similar spatial pattern but slightly less RMSE. Maybe just use 300 trees.  
LA
RMSE          7.1859116
RRMSE         0.3005263
IQR           7.3031196
rIQR          0.3335296
MAE           5.3181910
rMAE          0.2224en85
rsq           0.6859345
explained_var 0.6869203
cor           0.8322331
covprob95     0.3623711
covprob90     0.3211340
covprob50     0.1453608'


