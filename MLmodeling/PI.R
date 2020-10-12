source("INLA_util.R")
#install.packages("disttree", repos="http://R-Forge.R-project.org")
 
ipak <- function(pkg){
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE , repos='http://cran.muenster.r-project.org')
  sapply(pkg, require, character.only = TRUE)
}
packages <- c( "devtools", "dplyr","distrree","data.table" , "ggplot2" , "RColorBrewer", "xgboost",  "glmnet", "ranger", "randomForest","tidyr" ,"tibble","stargazer", "sf", "CAST", "caret", "quantregForest", "pdp", "h2o", "dismo")
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
 
mergedall$mean_value = sqrt(mergedall$mean_value) # normal 
n =1
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
  varall =df%>%dplyr::select(matches(varstring))
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
  sd_reg<- predict(quantRF,
                   data = x_p[test,],type = "quantiles",what=sd)
  mean_reg <- predict(quantRF,
                      data = x_p[test,], type = "quantiles",what=mean)

  hist(pred.distribution$predictions[5,])
  t.quant90 <- cbind( 
    pred.distribution$predictions[, "quantile= 0.05"],
    pred.distribution$predictions[, "quantile= 0.95"])
  ## distforest
  distf <- distforest(mean_value~.,  data = varall)
  pre =  predict(distf, newdata = x_p[test,], type = "parameter")
  w =  predict(distf, newdata = x_p[test,], type = "weights")
 
  mu_ = pre[["mu"]] 
  sigma_ =pre[["sigma"]]
  dist.q90 = cbind(mu_-1.64*sigma_, mu_+1.64*sigma_)
  
  
  
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
  
 

 
  #plot(inla_90[,1], ylim = c(min(y_denl_test)-1,max(y_denl_test)+1), col = "red", typ = "l")
  df1 = cbind(data.frame(rf_90), data.frame(dist.q90), id = 1:nrow(data.frame(rf_90)),y_denl_test^2)
 
  names(df1) = c("quantile_RF_L90", "quantile_RF_U90", "Distribution_RF_L90","Distribution_RF_U90", "id", "test")
df1 = df1%>%  gather(variable, value, -id, -test )
  df1$value = df1$value^2
ggplot(df1)+aes(x = id, y = value, colour = variable)+geom_line()+
    geom_point(aes(y=test), colour= "black")+
  scale_color_brewer(palette="Spectral")+labs(x = "test points", y = "NO2", colour = "prediction intervals")
,aes(x = 1:lentgh(test), y = test ))
  #l
  
  ines(inla_90[,2], col = "red")
plot(rf_90[,1]  ,ylim = c(min(y_denl_test)-1,max(y_denl_test)+1),  typ = "l")
   lines(rf_90[,2])
  lines(dist.q90[,1], col = "orange")
  lines(dist.q90[,2], col = "orange")
  points(y_denl_test)
  p = ggplot()+geom_sf(data = te, color = "red")+geom_sf(data=tr)
  plot(p)
  #qqnorm(sqrt(df$mean_value))
  #qqline(sqrt(df$mean_value), col = "steelblue", lwd = 2)
  #shapiro.test(sqrt(df$mean_value))
  shapiro.test(pred.distribution$predictions[5,])
  distcrps = crps(y = y_denl_test, family = "norm", mean = mu_, sd = sigma_) 
  plot(distcrps0, distcrps)
  regcrps = crps(y = y_denl_test, family = "norm", mean = as.vector(predictions(mean_reg)), sd =as.vector(sd_reg$predictions)) 
  plot(distcrps, regcrps)
  mean_reg$predictions
  #val <- APMtools::error_matrix(validation = dptest$real, prediction = dptest$pred_mean)
  #val
  #val <- APMtools::error_matrix(validation =y_denl_test, prediction = predictions(pred))
  #val
  
 