#devtools::install_github("mengluchu/APMtools") 
library(INLA)
library(APMtools)
library(dplyr)
library(h2o)
#h2o.shutdown(prompt = F)
h2o.init(nthreads=-1)

#' Creates triangulated mesh to fit a spatial model using INLA and SPDE
#' 
#' @param coo coordinates to create the mesh
#' @return mesh boject
fnConstructMesh <- function(coo){
  # meshbuilder()
  # offset: size of the inner and outer extensions around the data locations
  (offset1 <- 1/10*max(dist(coo)))
  (offset2 <- 1/2*max(dist(coo)))
  # max.edge: maximum allowed triangle edge lengths in the region and in the extension
  (maxedge1 <- 1/20*max(dist(coo)))
  (maxedge2 <- 1/3*max(dist(coo)))
  # cutoff: minimum allowed distance between points used to avoid building many small triangles around clustered locations
  (cutoff <- 1/10000*max(dist(coo)))
  mesh <- inla.mesh.2d(loc = coo, offset = c(offset1, offset2), cutoff = cutoff, max.edge = c(maxedge1, maxedge2))
  plot(mesh)
  points(coo, col = "red")
  print(mesh$n)
  return(mesh)
}



#' Fits a spatial model using INLA and SPDE
#' 
#' It creates a mesh using coordinates d$coox and d$cooy
#' Formula is passed in argument formula
#' It creates stk.full with data for estimation or data for estimation and prediction (if TFPOSTERIORSAMPLES is TRUE)
#' Calls \code{inla()} and returns list with result and stk.full
#' 
#' @param d Data frame with data for estimation that contains coordinates (coox, cooy), response variable (y) and covariates
#' @param formula Formula for the model
#' If \code{covnames} includes an intercept, \code{d} needs to have column of 1s for the intercept
#' @param dp Data frame with data for prediction that contains coordinates (coox, cooy), and covariates
#' If \code{covnames} includes an intercept, \code{dp} needs to have column of 1s for the intercept
#' If dp is NULL, dp will not used to construct stk.full
#' @param covnames Vector with the names of the intercept and covariates that are in the formula
#' @param TFPOSTERIORSAMPLES Boolean variable to call \code{inla()} with config = TFPOSTERIORSAMPLES.
#' If it config = TRUE we will get a res object with which we could call \code{inla.posterior.samples()}
#' @return list with the results of the fitted model, stk.full and mesh
fnFitModelINLA <- function(d, dp, formula, covnames, TFPOSTERIORSAMPLES){
  # Coordinates locations
  coo <- cbind(d$coox, d$cooy)
  # Mesh
  mesh <- fnConstructMesh(coo)
  # Building the SPDE model on the mesh
  spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
  # Index set
  indexs <- inla.spde.make.index("s", spde$n.spde)
  # Projection matrix
  A <- inla.spde.make.A(mesh = mesh, loc = coo)
  
  # Stack with data for estimation. Effects include intercept and covariates
  stk.e <- inla.stack(tag = "est", data = list(y = d$y), A = list(1, A), effects = list(d[, covnames, drop = FALSE], s = indexs))
  if(is.null(dp)){
  stk.full <- inla.stack(stk.e)
  }else{
  # Prediction coordinate locations and projection matrix
  coop <- cbind(dp$coox, dp$cooy)
  Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
  # stack  
  stk.p <- inla.stack(tag = "pred", data = list(y = NA), A = list(1, Ap), effects = list(dp[, covnames, drop = FALSE], s = indexs))
  stk.full <- inla.stack(stk.e, stk.p)
  }
  
  # Formula that is specified in the arguments
  st1 <- Sys.time()
  res <- inla(formula, data = inla.stack.data(stk.full), control.predictor = list(compute = TRUE, A = inla.stack.A(stk.full)),
              control.compute = list(config = TFPOSTERIORSAMPLES))
  st2 <- Sys.time()
  print(st2-st1)
  return(list(res, stk.full, mesh))
}




#' Computes the linear predictor from one of the samples of an object obtained with \code{inla.posterior.samples()
#' 
#' It retrieves the sample number \code{ite} from the object \code{psamples} that was obtained with \code{inla.posterior.samples()
#' For this sample, it extracts the betas for the coefficients in \code{covnames} and the values of the spatial field
#' Then it calculates beta*covariates + spatial effect
#' 
#' @param psamples Object obtained from \code{inla.posterior.samples() that contains a list with the samples
#' @param ite Number of sample from \code{psamples} that we want to use
#' @param res result object from an \code{inla()} call
#' @param mesh Triangulated mesh that was used to fit the model
#' @param dp Data frame with data for prediction that contains coordinates (coox, cooy), and covariates.
#' If \code{covnames} includes an intercept, \code{dp} needs to have column of 1s for the intercept
#' @param covnames Name of the coefficients in the formula (intercept and other covariates)
#' @return Data frame \code{dp} with added columns \code{pred_mean}, \code{pred_ll}, \code{pred_ul} and \code{excprob}
fnPredictFromPosteriorSample <- function(psamples, ite, res, mesh, dp, covnames){
  # Retrieve elements sample
  (contents <- res$misc$configs$contents)
  # betas for elements of covnames. covnames[1] is the first covariate (b0)
  id_firstcov <- grep(covnames[1], rownames(psamples[[ite]]$latent))
  betas <- psamples[[ite]]$latent[id_firstcov : (id_firstcov + (length(covnames)-1)), ]
  # spatial field
  id_s <- which(contents$tag == "s")
  id_s <- contents$start[id_s]:(contents$start[id_s] + contents$length[id_s] - 1)
  spatialfield <- psamples[[ite]]$latent[id_s]
  # spat <- lapply(ps, function(x) x$latent[id_s])
  # spat <- matrix(unlist(spat), ncol = length(id_s), byrow = T)
  # Multiply model matrix times betas + spatial effect
  coop <- cbind(dp$coox, dp$cooy)
  Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
  predictions <- as.matrix(dp[, covnames]) %*% betas + drop(Ap %*% spatialfield)
  return(predictions)
}

#' Get predictions from a result object obtained by fitting as spatial model using INLA and SPDE
#' 
#' @param res result object from an \code{inla()} call
#' @param stk.full stk.full object constructed during an \code{inla()} call
#' @param mesh Triangulated mesh that was used to fit the model
#' @param covnames Name of the coefficients in the formula (intercept and other covariates)
#' @param d Data frame with data for estimation that contains coordinates (coox, cooy), response variable (y) and covariates
#' If \code{covnames} includes an intercept, \code{d} needs to have column of 1s for the intercept
#' It can be NULL if predictions are added to \code{dp}
#' @param dp Data frame with data for prediction that contains coordinates (coox, cooy), and covariates
#' If \code{covnames} includes an intercept, \code{dp} needs to have column of 1s for the intercept
#' If dp is NULL, dp will not used to construct stk.full
#' It can be NULL if predictions are added to \code{d}
#' @param NUMPOSTSAMPLES number of samples to call \code{inla.posterior.samples()}
#' If NUMPOSTSAMPLES == -1, get the predictions directly from the "est" elements of res and add them to \code{d}
#' If NUMPOSTSAMPLES == 0, get the predictions directly from the "pred" elements of res and add them to \code{dp}
#' If NUMPOSTSAMPLES > 0, get the predictions  using \code{inla.posterior.samples()} and add them to \code{dp}.
#' If NUMPOSTSAMPLES > 0, \code{dp} may or may not have passed previously to \code{inla()}
#' @param cutoff_exceedanceprob cutoff value to compute exceedance probabilities P(theta > cutoff)
#' @return Data frame \code{d} or \code{dp} with added columns \code{pred_mean}, \code{pred_ll}, \code{pred_ul} and \code{excprob}
#' \code{pred_mean} is the posterior mean
#' \code{pred_ll} and \code{pred_ul} are the lower and upper limits of 95%, 90%, 50% credible intervals
#' \code{excprob} is the probability that hte prediction > cutoff value
fnGetPredictions <- function(res, stk.full, mesh, d, dp, covnames, NUMPOSTSAMPLES, cutoff_exceedanceprob){
  if(NUMPOSTSAMPLES == -1){
    index <- inla.stack.index(stk.full, tag = "est")$data
    d$excprob <- sapply(res$marginals.fitted.values[index],
                         FUN = function(marg){1-inla.pmarginal(q = cutoff_exceedanceprob, marginal = marg)})
    d$pred_mean <- res$summary.fitted.values[index, "mean"]
    d$pred_ll <- res$summary.fitted.values[index, "0.025quant"]
    d$pred_ul <- res$summary.fitted.values[index, "0.975quant"]
    d$pred_ll90 <- unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.05, marginal = marg)}))
    d$pred_ul90 <- unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.95, marginal = marg)}))
    d$pred_ll50 <- unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.25, marginal = marg)}))
    d$pred_ul50 <- unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.75, marginal = marg)}))
    dres <- d
  }
  if(NUMPOSTSAMPLES == 0){
    index <- inla.stack.index(stk.full, tag = "pred")$data
    dp$excprob <- sapply(res$marginals.fitted.values[index],
                         FUN = function(marg){1-inla.pmarginal(q = cutoff_exceedanceprob, marginal = marg)})
    dp$pred_mean <- res$summary.fitted.values[index, "mean"]
    dp$pred_ll <- res$summary.fitted.values[index, "0.025quant"]
    dp$pred_ul <- res$summary.fitted.values[index, "0.975quant"]
    dp$pred_ll90 <- unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.05, marginal = marg)}))
    dp$pred_ul90 <- unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.95, marginal = marg)}))
    dp$pred_ll50 <- unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.25, marginal = marg)}))
    dp$pred_ul50 <- unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.75, marginal = marg)}))
    dres <- dp
  }
  
  if(NUMPOSTSAMPLES > 0){
    psamples <- inla.posterior.sample(NUMPOSTSAMPLES, res)
    ps <- sapply(1:NUMPOSTSAMPLES, fnPredictFromPosteriorSample, psamples = psamples, res = res, mesh = mesh, dp = dp, covnames = covnames)
    dp$excprob <- apply(ps, 1, function(x){mean(x > cutoff_exceedanceprob)})
    dp$pred_mean <- rowMeans(ps)
    dp$pred_ll <- apply(ps, 1, function(x){quantile(x, 0.025)})
    dp$pred_ul <- apply(ps, 1, function(x){quantile(x, 0.975)})
    dp$pred_ll90 <- apply(ps, 1, function(x){quantile(x, 0.05)})
    dp$pred_ul90 <- apply(ps, 1, function(x){quantile(x, 0.95)})
    dp$pred_ll50 <- apply(ps, 1, function(x){quantile(x, 0.25)})
    dp$pred_ul50 <- apply(ps, 1, function(x){quantile(x, 0.75)})
    dres <- dp
  }
  return(dres) 
}

#' get h2o cv results from a h2o model and return a df
#' @param h2omodel h2omodel
#' 
get_h2o_pred = function(h2omodel) 
{
  cvpreds_id  = h2omodel@model$cross_validation_holdout_predictions_frame_id$name
  as.data.frame(h2o.getFrame(cvpreds_id))$predict
}


#' Calculation CV predicitons
#' @param d Data frame with data for estimation that contains coordinates (coox, cooy), response variable (y) and covariates
#' @param y_var name of response
#' @param prestring "regex" style string matching for covariates to use
#' @param training the index for training 

fnMLPredictionsCV = function(d, y_var= "y", training,  prestring =  "road|nightlight|population|temp|wind|trop|indu|elev|radi" ){
 
  x_var = d%>%dplyr::select(matches(prestring))%>%names() 
  td = d%>%dplyr::select(c(x_var, y_var))
  td = td[training, ]
  merged_hex = as.h2o(td)
  
  rf = h2o.randomForest(y = y_var , x= x_var, training_frame = merged_hex, ntrees =1000,  nfolds =10, keep_cross_validation_predictions = T, seed =1 )
  xgb = h2o.xgboost(y = y_var , x= x_var, training_frame = merged_hex, nfolds =10, ntrees =300, eta =  0.007, subsample = 0.7, max_depth = 6, gamma = 5, reg_lambda =2, reg_alpha = 0, keep_cross_validation_predictions = T, seed =1 )
  las = h2o.glm(y = y_var , x= x_var, training_frame = merged_hex, alpha = 1,  nfolds =10, keep_cross_validation_predictions = T, seed =1 )

  rf  =get_h2o_pred(rf)
  xgb  =get_h2o_pred(xgb)
  lasso  =get_h2o_pred(las)  
  
  h2o.removeAll()
  
  cbind(rf, xgb, lasso)
}

#' Calcualte predictions from all the obs. 
#' @param d Data frame with data for estimation that contains coordinates (coox, cooy), response variable (y) and covariates
#' @param training the index for training 
#' @param test the index for test
fnMLPredictionsAll = function(d, y_var="y", training, test, prestring =  "road|nightlight|population|temp|wind|trop|indu|elev|radi" ){

  x_var = d%>%dplyr::select(matches(prestring))%>%names() 
  td = d%>%dplyr::select(c(x_var, y_var))
  merged_hex = as.h2o(td[training, ])# train on training set
  
  new_hex = as.h2o(td[test, ]) # predict to test

  rf_all = h2o.randomForest(y = y_var, x= x_var, training_frame = merged_hex, ntrees =1000, seed =1 )
  rf  = as.data.frame(h2o.predict(object = rf_all, newdata = new_hex))$predict 
  xgb_all = h2o.xgboost(y = y_var, x= x_var, training_frame = merged_hex, ntrees =300, eta =  0.007, subsample = 0.7, max_depth = 6, gamma = 5, reg_lambda =2, reg_alpha = 0, seed =1 )
  xgb = as.data.frame(h2o.predict(object = xgb_all, newdata = new_hex))$predict 
  las_all = h2o.glm(y = y_var, x= x_var, training_frame = merged_hex, alpha = 1,  nfolds =1, seed=1 )
  lasso  = as.data.frame(h2o.predict(object = las_all, newdata = new_hex))$predict 
  h2o.removeAll()
  cbind(rf, xgb, lasso)
 
  }



#' Calculates cross-validation measures obtained by fitting a spatial model using INLA and SPDE
#' It uses d and dp because datasets for estimation and prediction have different covariates (predictions by cross-validation and using all data)
#'
#' @param n Number of iteration 
#' @param d Data frame with data for estimation that contains coordinates (coox, cooy), response variable (y) and covariates
#' @param dp Data frame with data for prediction that contains coordinates (coox, cooy), and covariates
#' If \code{covnames} includes an intercept, \code{d} needs to have column of 1s for the intercept
#' @param covnames Vector with the names of the intercept and covariates to be included in the formula
#' @return Vector with the cross-validation results
INLA_crossvali =  function(n, d, formula, covnames){
 
  print(n)
  # Split data
  smp_size <- floor(0.8 * nrow(d)) 
  set.seed(n)
  training <- sample(seq_len(nrow(d)), size = smp_size)
  test <- seq_len(nrow(d))[-training] 
  # Fit model
  dtraining <- d[training, ]
  dptest <- d[test, ]
  
  # INI CODE MENG
  # Add to d[, training] 3 variables with names lasso, rf, xgb that are predictions at locations coox and cooy calculated using cross-validation with the dataset d[, training]
  dtraining <- cbind(dtraining, fnMLPredictionsCV(d=d, training = training))
 
  # Add to dp[test, ] 3 variables with names lasso, rf, xgb that are predictions at locations coox and cooy calculated using all data in d[, training]
  
  dptest <- cbind(dptest,  fnMLPredictionsAll(d=d, training = training, test = test))
   
  # END CODE MENG
  
  
  # Fit model
  lres <- fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE)
  # Get predictions
  dptest <- fnGetPredictions(lres[[1]], lres[[2]], lres[[3]], dtraining, dptest, covnames, NUMPOSTSAMPLES = 0, cutoff_exceedanceprob = 30)

  # Goodness of fit
  val <- APMtools::error_matrix(validation = dptest$real, prediction = dptest$pred_mean)
  val <- c(val, cor = cor(dptest$real, dptest$pred_mean))
  print(val)
  (val <- c(val, covprob95 = mean(dptest$pred_ll <= dptest$real &  dptest$real <= dptest$pred_ul),  # 95% coverage probabilities
                 covprob90 = mean(dptest$pred_ll90 <= dptest$real &  dptest$real <= dptest$pred_ul90),
                 covprob50 = mean(dptest$pred_ll50 <= dptest$real &  dptest$real <= dptest$pred_ul50)))
  
  return(val)
} 

##################################################


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
  rf = h2o.randomForest(y = y_var , x= x_var, training_frame = merged_hex, ntrees =1000,  nfolds =10, keep_cross_validation_predictions = T, seed =1 )
  xgb = h2o.xgboost(y = y_var , x= x_var, training_frame = merged_hex, nfolds =10, ntrees =300, eta =  0.007, subsample = 0.7, max_depth = 6, gamma = 5, reg_lambda =2, reg_alpha = 0, keep_cross_validation_predictions = T, seed =1 )
  las = h2o.glm(y = y_var , x= x_var, training_frame = merged_hex, alpha = 1,  nfolds =10, keep_cross_validation_predictions = T, seed =1 )
  
  ens = h2o.stackedEnsemble(y = y_var, x= x_var, training_frame = merged_hex, base_models = c(rf, xgb, las))
  pred = h2o.predict(object =ens, newdata = dptest)$predict 
  
  val = APMtools::error_matrix(validation = ytest, prediction = pred)
  h2o.removeAll()
  return(val)
}
############## run #
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
 
 
# MODEL 2
# Model stacked generalization. I need predictions with the cross-validated test to fit the model, and predictions with the whole data to predict.
covnames <- c("lasso", "rf", "xgb")

formula <- as.formula("y ~ 0 + f(lasso, model = 'clinear', range = c(0, 1), initial = 1/3) +
                         f(rf, model = 'clinear', range = c(0, 1), initial = 1/3) +
                         f(xgb, model = 'clinear', range = c(0, 1), initial = 1/3) + f(s, model = spde)") 

# Cross-validation
VLA <- lapply(1:2, FUN = INLA_crossvali, d = d, formula = formula, covnames = covnames)
(VLA <- data.frame(LA = rowMeans(data.frame(VLA))))

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
rMAE          0.2224385
rsq           0.6859345
explained_var 0.6869203
cor           0.8322331
covprob95     0.3623711
covprob90     0.3211340
covprob50     0.1453608'

ensembleh2o <- lapply(1:2, FUN =ensemble, d = d )
save(ensembleh2o, file = "ensembleh2o.rda")
save(VLA, file = "VLA.rda")
