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
ensembleh2o <- lapply(1:20, FUN =ensemble, d = d )

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

ensemble= function(d, n, y_var= "y", prestring =  "road|nightlight|population|temp|wind|trop|indu|elev|radi" ){
  print(n)
  # Split data
  smp_size <- floor(0.8 * nrow(d)) 
  set.seed(n)
  training <- sample(seq_len(nrow(d)), size = smp_size)
  test <- seq_len(nrow(d))[-training] 
  # Fit model
  
  dptest <- d[test, ]
  ytest = dptest[,y_var]
  x_var = d%>%dplyr::select(matches(prestring))%>%names() 
  td = d%>%dplyr::select(c(x_var, y_var))
  td = td[training, ]
  
  merged_hex = as.h2o(td)
  dptest= as.h2o(dptest)
  rf = h2o.randomForest(y = y_var , x= x_var, training_frame = merged_hex, ntrees =2,  nfolds =2, keep_cross_validation_predictions = T, seed =1 )
  xgb = h2o.xgboost(y = y_var , x= x_var, training_frame = merged_hex, nfolds =2, ntrees =3, eta =  0.007, subsample = 0.7, max_depth = 6, gamma = 5, reg_lambda =2, reg_alpha = 0, keep_cross_validation_predictions = T, seed =1 )
  las = h2o.glm(y = y_var , x= x_var, training_frame = merged_hex, alpha = 1,  nfolds =2, keep_cross_validation_predictions = T, seed =1 )
  
  ens = h2o.stackedEnsemble(y = y_var, x= x_var, training_frame = merged_hex, base_models = c(rf, xgb, las))
  dtest <- fnMLPredictionsAll(d=d, training = training, test = test)
  dtesth2o = as.h2o(dtest)
  print(h2o.performance(ens, newdata = dtesth2o))
  print(h2o.predict(ens, newdata = dtesth2o))
  pred = as.data.frame(h2o.predict(object =ens, newdata = dtesth2o))$predict 
  
  val = APMtools::error_matrix(validation = ytest, prediction = pred)
  
  h2o.removeAll()
  return(val)
}
