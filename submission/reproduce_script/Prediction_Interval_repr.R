# Prediction intervals of QRF, QRFLA, DF, INLA
# figure 4-6 manuscript

#install.packages("disttree", repos="http://R-Forge.R-project.org")

ipak <- function(pkg){
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE , repos='http://cran.muenster.r-project.org')
  sapply(pkg, require, character.only = TRUE)
}
packages <- c( "devtools", "dplyr","disttree","data.table" , "scoringRules","ggplot2" , "RColorBrewer", 
               "xgboost",  "glmnet", "ranger", "randomForest","tidyr" ,"tibble","stargazer", "sf",
               "quantregForest", "disttree")
ipak(packages)

install_github("mengluchu/APMtools") 
library(APMtools)
ls("package:APMtools") 

resolution =100 # resolution of grid
nboot = 12  # number of bootstraps

######
#data preprocessing
#

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
num.trees = 1000  
set.seed(1)
smp_size <- floor(0.8 * nrow(df)) 

######### split training and test ######################
training <- sample(seq_len(nrow(df)), size = smp_size)
test = seq_len(nrow(df))[-training] 
#=======================================================

y_denl = df[,y_var]
y_denl_test = y_denl[test] 
y_denl_train = y_denl[training]

varall =df%>%dplyr::select(matches(varstring))
x_p = df%>%dplyr::select(matches(varstring))%>%dplyr::select(-y_var)
x_train =x_p[training,]
x_test = x_p[test,]

###===========================
## QRF. method 1, (used in the manuscript) use ranger quanreg = T
###===============================

quantRF <- ranger(x = x_train,
                  y = y_denl_train, mtry = NULL, num.trees = 1000,
                  quantreg = T, min.node.size = 10) 
# compute predictions (mean) for each validation site
pred <- predict(quantRF, data = x_test, what = mean)
## predict 0.01, 0.02,..., 0.99 quantiles for validation data
pred.distribution <- predict(quantRF,
                             data = x_test, 
                             type = "quantiles",
                             quantiles = seq(0.01, 0.99, by = 0.01)) # get quantiles

pred.distribution$predictions
sd_reg<- predict(quantRF, data = x_test,type = "quantiles",what=sd)
mean_reg <- predict(quantRF, data = x_test, type = "quantiles",what=mean)
rf_90 <- cbind( 
  pred.distribution$predictions[, "quantile= 0.05"],
  pred.distribution$predictions[, "quantile= 0.95"])
 
#===================================================
##QRF, method 2 (not used in the manuscript): test if the result is the same as the method 1
#-=========================================================
qrf <- quantregForest(x=x_train,
                      y = y_denl_train, nodesize=10, ntrees=1000)
QRF_U90 <- predict(qrf, x_test, what = 0.95)
QRF_L90 <- predict(qrf, x_test, what = 0.05)
 
#====================================================================================
# use QRFLA (used in the manuscript, implemented in APMtools )
#=================================================================================
RF <- ranger(x = x_train,
             y = y_denl_train, mtry = NULL, num.trees = 1000, min.node.size = 10,,importance = "permutation")
RF$variable.importance%>%sort(decreasing = T) 

allp = predict(RF,x_train,type = "response", predict.all = T)%>%predictions #get all the tree predictions, instead of the mean
rpre= predict(RF,x_test, predict.all=T)%>%predictions # get all the tree predictions

cvfit = glmnet::cv.glmnet(allp,y_denl_train, 
                          type.measure = "mse", standardize = TRUE, alpha = 1, 
                          lower.limit = 0, nfolds = 10)  
# aggregate using a regularization, here lasso, you can also do elastic net, training alpha or specify alpha between 0 and 1
print(sum(coef(cvfit)[-1]!= 0))
# we can also plot it, using a tool from APMtools
Ls= lassoselected(cvfit)
Ls_num= as.vector(sapply(Ls, function(x) as.numeric(substr(x, start =2, stop = nchar(x)))))
n = length(test)
# aggregating trees using lasso, compare with original random forest, obtained better results

reduced_rf = rpre[,Ls_num] # 62 trees
rrf_u90 = apply(reduced_rf, 1, quantile, 0.95)
rrf_l90 = apply(reduced_rf, 1, quantile, 0.05)
rrf_mean = apply(reduced_rf, 1,  mean)
rrf_sd= apply(reduced_rf, 1,   sd)
#================================================
# test: all predictions of ranger vs. quanreg. 
#================================================

#  predict.all =T 
#  get mean for each tree, the quantiles use a random observation y for each node and tree, if a node has many obs. for one tree, it doesnt matter because next tree may select a different y, in this sense it is a bit like bootstrapping.
#  they are closely correlated, but not equal. ###########################################
### try the original quantile rf
### reduced RF used LASSO selected trees.

terminal.nodes <- predict(RF, x_train, type = "terminalNodes")$predictions + 1
random.node.values <- matrix(nrow = max(terminal.nodes), ncol = num.trees)

## Select one random obs per node and tree
for (tree in 1:num.trees){
  idx <- sample(1:n, n)
  random.node.values[terminal.nodes[idx, tree], tree] <- y_denl_train[idx]
}
terminal.nodes <- predict(RF, x_test, type = "terminalNodes")$predictions + 1
node.values <- 0 * terminal.nodes
for (tree in 1:num.trees) {
  node.values[, tree] <-  random.node.values[terminal.nodes[, tree], tree]
}  
rf_u90 <-  apply(node.values, 1, quantile, 0.95, na.rm=TRUE)

rrf2_u90 = apply(node.values[,Ls_num], 1, quantile, 0.95,na.rm=TRUE)
rrf2_l90 = apply(node.values[,Ls_num], 1, quantile, 0.05,na.rm=TRUE)
rrf2_mean = apply(node.values[,Ls_num], 1, mean,na.rm=TRUE)
rrf2_sd = apply(node.values[,Ls_num], 1, sd ,na.rm=TRUE)

plot(rrf2_l90, rrf_l90)

error_matrix(y_denl_test,rrf2_mean) # QRF
error_matrix(y_denl_test,rrf_mean) # Lasso + RF 
error_matrix(y_denl_test,predict(RF,x_test,type = "response")%>%predictions) #RF

## compare 

cor(apply(rpre, 1, quantile, 0.95), pred.distribution$predictions[, "quantile= 0.95"])
cor(rf_u90, pred.distribution$predictions[, "quantile= 0.95"])
cor(rf_u90, apply(rpre, 1, quantile, 0.95))

plot(apply(rpre, 1, quantile, 0.95), pred.distribution$predictions[, "quantile= 0.95"])
abline(b= 1, a =0)

#### end test ################### 


##============================================
## distforest (DF), used in manuscript
##=================================================

distf <- distforest(mean_value~.,  data = varall)
pre =  predict(distf, newdata = x_test, type = "parameter")
w =  predict(distf, newdata = x_test, type = "weights")

mu_ = pre[["mu"]] 
sigma_ =pre[["sigma"]]
dist.q90 = cbind(mu_-1.64*sigma_, mu_+1.64*sigma_)

# plot
df1 = cbind(data.frame(rf_90), data.frame(dist.q90), rrf_l90,rrf_u90,rrf2_l90,rrf2_u90,id = 1:nrow(data.frame(rf_90)),y_denl_test)
names(df1) = c("QRF_L90","QRF_U90", "DF_L90","DF_U90", "QRFLA_L90","QRFLA_U90", "reduced2_RF_L90","reduced2_RF_U90","id", "test")

#============
#DF vs QRF, figure 4 manuscript
#===============
dftest = df1%>%  gather(variable, value, -id, -test,-QRFLA_L90,-QRFLA_U90, -reduced2_RF_L90,-reduced2_RF_U90)

ggplot(dftest)+aes(x = id, y = value, colour = variable)+geom_line()+
  geom_point(aes(y=test), colour= "black")+
  scale_color_brewer(palette="Spectral")+labs(x = "test points", y = "NO2", colour = "prediction intervals")


#============
# QRF vs QRFLA, figure 6 manuscript
#===============

dftest = df1%>%  gather(variable, value, -id, -test, -DF_L90, -DF_U90, -reduced2_RF_L90,-reduced2_RF_U90)
ggplot(dftest)+aes(x = id, y = value, colour = variable)+geom_line()+
  geom_point(aes(y=test), colour= "black")+
  scale_color_brewer(palette="Spectral")+labs(x = "test points", y = "NO2", colour = "prediction intervals")


#ggsave("qrf_df.png")
#ggsave("qrf_qrfla.png")

 
# CRPS evaluration prob. forecast with scoring ranks: https://arxiv.org/pdf/1709.04743.pdf
distcrps = crps(y = y_denl_test, family = "norm", mean = mu_, sd = sigma_) 
rrfcrps = crps(y = y_denl_test, family = "norm", mean = rrf_mean, sd = rrf_sd) # lasso+rf
rrf2crps = crps(y = y_denl_test, family = "norm", mean = rrf2_mean, sd = rrf2_sd) 

rrfcrps2 = crps_sample(y = y_denl_test, reduced_rf, method = "edf") # lasso+rf "kde" is better the "edf" and better than assuming gaussian . 
summary(rrfcrps)

rrfcrps3 = crps_sample(y = y_denl_test, rpre, method = "edf") # rf,   "kde" is similar to  "edf" and better than assuming gaussian . 
summary(rrfcrps3)
regcrps = crps(y = y_denl_test, family = "norm", mean = as.vector(predictions(mean_reg)), sd =as.vector(sd_reg$predictions)) # mean and sd of the ranger quantile RF. 

summary(cbind(rrfcrps,regcrps, rrf2crps))
plot(distcrps, regcrps)
mean_reg$predictions
 

#========================  
# INLA
#==========================
source("~/INLA_util.R") 
d = df
covnames = c("nightlight_450", "population_1000", "population_3000",
             "road_class_1_5000", "road_class_2_100", "road_class_3_300",  
             "trop_mean_filt", "road_class_1_100")
d <- df[, c("mean_value", "Longitude", "Latitude", covnames)]

# covnames0 <- NULL
covnames <- c("b0", covnames)  # covnames is intercept and covnames0

d$y <-d$mean_value # response
d$coox <- d$Longitude
d$cooy <- d$Latitude
d$b0 <- 1 # intercept
d$real <- d$y
dp <- d
dtraining <- d[training, ]
dptest <- dp[test, ]

formula = as.formula(paste0('y ~ 0 + ', paste0(covnames, collapse = '+'), " + f(s, model = spde)"))

lres = fnFitModelINLA(d= dtraining, dp = dptest, covnames, formula = formula, TFPOSTERIORSAMPLES = TRUE, family = "gaussian")
res = lres[[1]]
#res = improved.result
stk.full = lres[[2]]
mesh = lres[[3]]
dres = fnGetPredictions(res, stk.full, mesh, d =dtraining, dp=dptest, covnames, NUMPOSTSAMPLES = 0, cutoff_exceedanceprob = 30)

lres = fnFitModelINLA(d= dtraining, dp = dptest, covnames, formula = formula, TFPOSTERIORSAMPLES = TRUE, family = "Gamma")
res = lres[[1]]
#res = improved.result
stk.full = lres[[2]]
mesh = lres[[3]]
dres2 = fnGetPredictions(res, stk.full, mesh, d =dtraining, dp=dptest, covnames, NUMPOSTSAMPLES = 0, cutoff_exceedanceprob = 30)


# Get predictions. NUMPOSTSAMPLES = -1 calculated with estimation data, 0 with prediction data, 1 with inla.posterior.samples()
#  ysim = rnorm(n = nrow(dtest), mean = dres$pred_mean, sd = dres$pred_sd)
#plot(ysim, typ = "l")
inlacrps = crps(y =dptest$real, family = "norm", mean = dres$pred_mean, sd =dres$pred_sd) 
# plot(y_denl_test)
#lines(dres$pred_ul)
# lines(dres$pred_ll, col = "red")

df1 = data.frame(cbind( dres$pred_ll90, dres$pred_ul90,  dres2$pred_ll90, dres2$pred_ul90, id = 1:nrow(data.frame(rf_90)),y_denl_test ))

names(df1) = c( "INLA_L90", "INLA_U90","INLA-G_L90", "INLA-G_U90","id", "test")


#====================
# plot INLA prediction interval, figure 5 manuscript
#=======================


df1 = df1%>%  gather(variable, value, -id, -test)

ggplot(df1)+aes(x = id, y = value, colour = variable)+geom_line()+
  geom_point(aes(y=test), colour= "black")+
  scale_color_brewer(palette="Set1")+labs(x = "test points", y = "NO2", colour = "prediction intervals")
#
#ggsave("INLA_pred.png")

#==============================
# plot training and test, not shown in manuscript
#=============================
coordi = c(mergedall$Longitude, mergedall$Latitude)
p <- matrix(coordi, ncol = 2)
tr = p[training,] %>%  
  sf::st_multipoint() 

te = p[test,] %>% 
  st_multipoint()

inla_90= cbind(dptest$pred_ll90,dptest$pred_ul90)
plot(inla_90[,1], ylim = c(min(y_denl_test)-1,max(y_denl_test)+1), col = "red", typ = "l")

p = ggplot()+geom_sf(data = te, color = "red")+geom_sf(data=tr)
plot(p)