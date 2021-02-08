# tried 
# ranger quantreg=T     quantile_RF
# quantregForest        QRF
# ranger predict.all =T, then use quantiles 
# Manually implement quantreg=T with ranger type = "terminalNodes" 
# The differences between them are very small
# ALso reduced random forest. The results seem to be better. 

#install.packages("disttree", repos="http://R-Forge.R-project.org")
 
ipak <- function(pkg){
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE , repos='http://cran.muenster.r-project.org')
  sapply(pkg, require, character.only = TRUE)
}
packages <- c( "devtools", "dplyr","distrree","data.table" , "scoringRules","ggplot2" , "RColorBrewer", 
               "xgboost",  "glmnet", "ranger", "randomForest","tidyr" ,"tibble","stargazer", "sf",
               "quantregForest")
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

# do it for 100 m resolution, get rid of 25m resolution variables. 
if (resolution ==100)
{
  mergedall = mergedall%>%dplyr::select(-c(industry_25,industry_50,road_class_1_25,road_class_1_50,road_class_2_25,road_class_2_50,   road_class_3_25,road_class_3_50))
} 
df = mergedall 
 
#mergedall$mean_value = sqrt(mergedall$mean_value) # normalise for normal distribution 

# to be written into a function
 
set.seed(1)
smp_size <- floor(0.8 * nrow(df)) 
  
training <- sample(seq_len(nrow(df)), size = smp_size)
test = seq_len(nrow(df))[-training] 
  
y_denl = df[,y_var]
y_denl_test = y_denl[test] 
varall =df%>%dplyr::select(matches(varstring))
x_p = df%>%dplyr::select(matches(varstring))%>%dplyr::select(-y_var)

############  
## 1. use ranger quanreg = T
######  
quantRF <- ranger(x = x_p[training,],
                    y = y_denl[training], mtry = NULL, num.trees = 1000,
                    quantreg = T, min.node.size = 10) 
  # compute predictions (mean) for each validation site
pred <- predict(quantRF, data = x_p[test,], what = mean)
  ## predict 0.01, 0.02,..., 0.99 quantiles for validation data
pred.distribution <- predict(quantRF,
                               data = x_p[test,], 
                               type = "quantiles",
                               quantiles = seq(0.01, 0.99, by = 0.01)) # get quantiles
 
sd_reg<- predict(quantRF,
                   data = x_p[test,],type = "quantiles",what=sd)
mean_reg <- predict(quantRF,
                      data = x_p[test,], type = "quantiles",what=mean)
t.quant90 <- cbind( 
    pred.distribution$predictions[, "quantile= 0.05"],
    pred.distribution$predictions[, "quantile= 0.95"])
  
###################### quantregForest: should be the same as the ranger quantreg =T

qrf <- quantregForest(x=x_p[training,],
                        y = y_denl[training], nodesize=10, ntrees=1000)
QRF_U90 <- predict(qrf, x_p[test,], what = 0.95)
QRF_L90 <- predict(qrf, x_p[test,], what = 0.05)
#  plot( conditionalQuantiles,pred.distribution$predictions[, "quantile= 0.95"])
  #hist(pred.distribution$predictions[5,])
 
#################################  use Lasso to aggregate random forest trees, implemented in APMtools 
RF <- ranger(x = x_p[training,],
             y = y_denl[training], mtry = NULL, num.trees = 1000, min.node.size = 10)

allp = predict(RF,x_p[training,],type = "response", predict.all = T)%>%predictions #get all the tree predictions, instead of the mean
rpre= predict(RF,x_p[test,], predict.all=T)%>%predictions # get all the tree predictions


  #allp = predict(RF,x_p[training,],type = "terminalNodes")$predictions
 # length(unique(allp[,39]))
  #plot(  apply(allp, 1, quantile, 0.05),pred.distribution$predictions[, "quantile= 0.05"])
  #apply(allp, 1, quantile, 0.95)-pred.distribution$predictions[, "quantile= 0.95"]
  #plot(  apply(allp, 1, quantile, 0.95),pred.distribution$predictions[, "quantile= 0.95"], col = "red")
  #abline(a=0,b = 1)
   
cvfit = glmnet::cv.glmnet(allp,y_denl[training], 
                              type.measure = "mse", standardize = TRUE, alpha = 1, 
                               lower.limit = 0, nfolds = 10)  

# aggregate using a regularization, here lasso, you can also do elastic net, training alpha or specify alpha between 0 and 1
print(sum(coef(cvfit)[-1]!= 0))
# we can also plot it, using a tool from APMtools
Ls= Lassoselected(cvfit)
Ls_num= as.vector(sapply(Ls, function(x) as.numeric(substr(x, start =2, stop = nchar(x)))))
n = length(test)

########## test: all predictions of ranger vs. quanreg. 
#########somehow i get different quantiles, using the same method. 
#  predict.all =T 
#  get mean for each tree, the quantiles use a random observation y for each node and tree, if a node has many obs. for one tree, it doesnt matter because next tree may select a different y, in this sense it is a bit like bootstrapping.
#  they are closely correlated, but not equal. ###########################################
### try the original quantile rf
### reduced RF used LASSO selected trees.

reduced_rf = rpre[,Ls_num]
y_denl_train = y_denl[training]

num.trees = 1000
y_denl_train = y_denl[training]
    
terminal.nodes <- predict(RF, x_p[training,], type = "terminalNodes")$predictions + 1
random.node.values <- matrix(nrow = max(terminal.nodes), ncol = num.trees)
     
## Select one random obs per node and tree
for (tree in 1:num.trees){
      idx <- sample(1:n, n)
      random.node.values[terminal.nodes[idx, tree], tree] <- y_denl_train[idx]
    }
    terminal.nodes <- predict(RF, x_p[test,], type = "terminalNodes")$predictions + 1
    node.values <- 0 * terminal.nodes
    for (tree in 1:num.trees) {
      node.values[, tree] <-  random.node.values[terminal.nodes[, tree], tree]
    }  
    rf_u90 <-  apply(node.values, 1, quantile, 0.95, na.rm=TRUE)
    
    rrf2_u90 = apply(node.values[,Ls_num], 1, quantile, 0.95,na.rm=TRUE)
    rrf2_l90 = apply(node.values[,Ls_num], 1, quantile, 0.05,na.rm=TRUE)
    rrf2_mean = apply(node.values[,Ls_num], 1, mean,na.rm=TRUE)
    rrf2_sd = apply(node.values[,Ls_num], 1, sd ,na.rm=TRUE)
   
     # aggregating trees using lasso, compare with original random forest, obtained better results
    rrf_u90 = apply(reduced_rf, 1, quantile, 0.95)
    rrf_l90 = apply(reduced_rf, 1, quantile, 0.05)
    
    rrf_mean = apply(reduced_rf, 1,  mean)
    rrf_sd= apply(reduced_rf, 1,   sd)
    
    plot(rrf2_l90, rrf_l90)
    
    error_matrix(y_denl_test,rrf2_mean)
    error_matrix(y_denl_test,rrf_mean)
    error_matrix(y_denl_test,predict(RF,x_p[test,],type = "response")%>%predictions)
    
    #########################################################################
   
    apply(rpre, 1, quantile, 0.95)       
    pred.distribution$predictions[, "quantile= 0.05"]
    cor(apply(rpre, 1, quantile, 0.95), pred.distribution$predictions[, "quantile= 0.95"])
    cor(rf_u90, pred.distribution$predictions[, "quantile= 0.95"])
    cor(rf_u90, apply(rpre, 1, quantile, 0.95))

    plot(apply(rpre, 1, quantile, 0.95), pred.distribution$predictions[, "quantile= 0.95"])
      abline(b= 1, a =0)
   ##############################################################################
    
       
  ##############
  ## distforest
  ################
      
 distf <- distforest(mean_value~.,  data = varall)
  pre =  predict(distf, newdata = x_p[test,], type = "parameter")
  w =  predict(distf, newdata = x_p[test,], type = "weights")
 
  mu_ = pre[["mu"]] 
  sigma_ =pre[["sigma"]]
  dist.q90 = cbind(mu_-1.64*sigma_, mu_+1.64*sigma_)

#compare
  rf_90 = t.quant90
  df1 = cbind(data.frame(rf_90), data.frame(dist.q90), rrf_l90,rrf_u90,QRF_L90,QRF_U90,id = 1:nrow(data.frame(rf_90)),y_denl_test )
 
  names(df1) = c("quantile_RF_L90","quantile_RF_U90", "Distribution_RF_L90","Distribution_RF_U90", "reduced_RF_L90","reduced_RF_U90", "reduced2_RF_L90","reduced2_RF_U90","id", "test")
# plot
    df1 = df1%>%  gather(variable, value, -id, -test, -quantile_RF_L90, -quantile_RF_U90, -Distribution_RF_L90, -Distribution_RF_U90)

  ggplot(df1)+aes(x = id, y = value, colour = variable)+geom_line()+
    geom_point(aes(y=test), colour= "black")+
  scale_color_brewer(palette="Spectral")+labs(x = "test points", y = "NO2", colour = "prediction intervals")
   #
  ggsave("dist_vs_qrf_notsq.png")
 
plot(rf_90[,1]  ,ylim = c(min(y_denl_test)-1,max(y_denl_test)+1),  typ = "l")
   lines(rf_90[,2])
  lines(dist.q90[,1], col = "orange")
  lines(dist.q90[,2], col = "orange")
  points(y_denl_test)

  
  ############

  #qqnorm(sqrt(df$mean_value))
  #qqline(sqrt(df$mean_value), col = "steelblue", lwd = 2)
  #shapiro.test(sqrt(df$mean_value))
shapiro.test(pred.distribution$predictions[5,])
  
  # CRPS evaluration prob. forecast with scoring ranks: https://arxiv.org/pdf/1709.04743.pdf
distcrps = crps(y = y_denl_test, family = "norm", mean = mu_, sd = sigma_) 
rrfcrps = crps(y = y_denl_test, family = "norm", mean = rrf_mean, sd = rrf_sd) 
rrf2crps = crps(y = y_denl_test, family = "norm", mean = rrf2_mean, sd = rrf2_sd) 
  
regcrps = crps(y = y_denl_test, family = "norm", mean = as.vector(predictions(mean_reg)), sd =as.vector(sd_reg$predictions)) 
  
summary(cbind(distcrps,rrfcrps,regcrps, rrf2crps))
plot(distcrps, regcrps)
mean_reg$predictions
  #val <- APMtools::error_matrix(validation = dptest$real, prediction = dptest$pred_mean)
  #val
  #val <- APMtools::error_matrix(validation =y_denl_test, prediction = predictions(pred))
  #val

  ### xgb: rf_lasso is better than rf, but xgb still is the best
  xgb = xgboost(data = as.matrix(x_p[training,]),
                label = y_denl[training],  max_depth =6, gamma=5, eta =0.007, nrounds =1000, lambda = 2, alpha = 0, subsample = 0.7 )
  
  xgbpre = predict(xgb, as.matrix(x_p[test,]))
  error_matrix(y_denl_test, xgbpre)
  
  
 ##################### 
# INLA
 ################
  # source("~/Documents/GitHub/uncertainty/MLmodeling//INLA/INLA_util.R") 
source("https://raw.githubusercontent.com/mengluchu/uncertainty/master/MLmodeling//INLA/INLA_util.R")
  
  covnames0 <- c("nightlight_450", "population_1000", "population_3000",
                 "road_class_1_5000", "road_class_2_100", "road_class_3_300", "trop_mean_filt",
                 "road_class_3_3000", "road_class_1_100", "road_class_3_100"
                 , "road_class_3_5000", "road_class_1_300", "road_class_1_500",
                 "road_class_2_1000", "nightlight_3150", "road_class_2_300", "road_class_3_1000", "temperature_2m_7")
  
  
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
  #plot(inla_90[,1], ylim = c(min(y_denl_test)-1,max(y_denl_test)+1), col = "red", typ = "l")
  
  p = ggplot()+geom_sf(data = te, color = "red")+geom_sf(data=tr)
  plot(p)