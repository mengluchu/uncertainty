source("INLA_util.R")

ipak <- function(pkg){
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE , repos='http://cran.muenster.r-project.org')
  sapply(pkg, require, character.only = TRUE)
}
packages <- c( "devtools", "dplyr","data.table" , "ggplot2" , "RColorBrewer", "xgboost",  "glmnet", "ranger", "randomForest","tidyr" ,"tibble","stargazer", "sf", "CAST", "caret", "quantregForest", "pdp", "h2o", "dismo")
ipak(packages)
install_github("mengluchu/APMtools") 
library(APMtools)
ls("package:APMtools") 
 

resolution =100 # resolution of grid
nboot = 12  # number of bootstraps
y_var = "mean_value"
prestring =  "road|nightlight|population|temp|wind|trop|indu|elev|radi"
varstring = paste(prestring,y_var,sep="|")
 
mergedall = read.csv("https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/DENL17_uc.csv")
if (resolution ==100)
{
  mergedall = mergedall%>%dplyr::select(-c(industry_25,industry_50,road_class_1_25,road_class_1_50,road_class_2_25,road_class_2_50,   road_class_3_25,road_class_3_50))
}  
for (n in 1:10){
  covnames0 <- c("nightlight_450", "population_1000", "population_3000",
                 "road_class_1_5000", "road_class_2_100", "road_class_3_300", "trop_mean_filt",
                 "road_class_3_3000", "road_class_1_100", "road_class_3_100"
                 , "road_class_3_5000", "road_class_1_300", "road_class_1_500",
                 "road_class_2_1000", "nightlight_3150", "road_class_2_300", "road_class_3_1000", "temperature_2m_7")
  
  
  df = mergedall
  set.seed(n)
  smp_size <- floor(0.8 * nrow(df)) 
  
  training <- sample(seq_len(nrow(df)), size = smp_size)
  test = seq_len(nrow(df))[-training] 
  
  y_denl = df[,y_var]
  y_denl_test = y_denl[test] 
  x_p = df%>%dplyr::select(matches(varstring))%>%dplyr::select(-y_var)
  
  quantRF <- ranger(x = x_p[training,],
                    y = y_denl[training], mtry = NULL, num.trees = 1000,
                    quantreg = T) 
  # compute predictions (mean) for each validation site
  pred <- predict(quantRF, data = x_p[test,], what = mean)
  
  
  ## ----investigate-single-point,echo=FALSE,fig.pos='!h',fig.height=5,fig.width=4,fig.align='center', out.width='0.4\\textwidth',fig.cap= "Histogram of predictive distribution for one single prediction point (dotted lines: 90 \\% prediction interval, dashed line: mean prediction)."----
  
  ## predict 0.01, 0.02,..., 0.99 quantiles for validation data
  pred.distribution <- predict(quantRF,
                               data = x_p[test,], 
                               type = "quantiles",
                               quantiles = seq(0.01, 0.99, by = 0.01))
  
  
  
  t.quant90 <- cbind( 
    pred.distribution$predictions[, "quantile= 0.05"],
    pred.distribution$predictions[, "quantile= 0.95"])
  # INLA
  
  d <- df[, c("mean_value", "Longitude", "Latitude", covnames0)]
  
  # covnames0 <- NULL
  covnames <- c("b0", covnames0)  # covnames is intercept and covnames0
  
  d$y <- d$mean_value # response
  d$coox <- d$Longitude
  d$cooy <- d$Latitude
  d$b0 <- 1 # intercept
  d$real <- d$y
  dp <- d
  dtraining <- d[training, ]
  dptest <- dp[test, ]
  # Fit model
  lres <- fnFitModelINLA(dtraining, dptest, covnames, TFPOSTERIORSAMPLES = FALSE, formulanew = NULL)
  # Get predictions
  dptest <- fnGetPredictions(lres[[1]], lres[[2]], lres[[3]], dtraining, dptest, covnames, NUMPOSTSAMPLES = 0, cutoff_exceedanceprob = 30)
  coordi = c(mergedall$Longitude, mergedall$Latitude)
  p <- matrix(coordi, ncol = 2)
  tr = p[training,] %>%  
    sf::st_multipoint() 
  
  te = p[test,] %>% 
    st_multipoint()
  
  inla_90= cbind(dptest$pred_ll90,dptest$pred_ul90)
  rf_90 = t.quant90
  plot(inla_90[,1], ylim = c(-30,65), col = "red", typ = "l")
  lines(rf_90[,1])
  lines(inla_90[,2], col = "red")
  lines(rf_90[,2])
  p = ggplot()+geom_sf(data = te, color = "red")+geom_sf(data=tr)
  plot(p)
  #val <- APMtools::error_matrix(validation = dptest$real, prediction = dptest$pred_mean)
  #val
  #val <- APMtools::error_matrix(validation =y_denl_test, prediction = predictions(pred))
  #val
  
}