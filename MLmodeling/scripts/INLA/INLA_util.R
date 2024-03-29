# functions used in geostatistical modelling in the manuscript
#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#install.packages("INLA")
#devtools::install_github("mengluchu/APMtools") 
library(APMtools)
library(INLA)

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
  #plot(mesh)
 # points(coo, col = "red")
#  print(mesh$n)
  return(mesh)
}

#' Fits a spatial model using INLA and SPDE
#' 
#' It creates a mesh using coordinates d$coox and d$cooy
#' Formula is intercept + covnames + spatial effect unless \code{formulanew} is specified
#' It creates stk.full with data for estimation or data for estimation and prediction (if TFPOSTERIORSAMPLES is TRUE)
#' Calls \code{inla()} and returns list with result and stk.full
#' 
#' @param d Data frame with data for estimation that contains coordinates (coox, cooy), response variable (y) and covariates
#' If \code{covnames} includes an intercept, \code{d} needs to have column of 1s for the intercept
#' @param dp Data frame with data for prediction that contains coordinates (coox, cooy), and covariates
#' If \code{covnames} includes an intercept, \code{dp} needs to have column of 1s for the intercept
#' If dp is NULL, dp will not used to construct stk.full
#' @param covnames Vector with the names of the intercept and covariates to be included in the formula
#' @param TFPOSTERIORSAMPLES Boolean variable to call \code{inla()} with config = TFPOSTERIORSAMPLES.
#' If it config = TRUE we will get a res object with which we could call \code{inla.posterior.samples()}
#' @param formulanew A string with the formula. If formulanew is NULL, the formula will be constructed as intercept + covnames + spatial effect
#' @return list with the results of the fitted model, stk.full and mesh
#' 
#' 
#' 
#' 
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
fnFitModelINLA = function(d, dp, formula, covnames, TFPOSTERIORSAMPLES, family = "gaussian"){
  # Coordinates locations
  coo = cbind(d$coox, d$cooy)
  # Mesh
  mesh = fnConstructMesh(coo)
  # Building the SPDE model on the mesh
  #  spde = inla.spde2.pcmatern(mesh, prior.range = c(500, .5), prior.sigma = c(2, 0.01))
  spde = inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
  # Index set
  indexs = inla.spde.make.index("s", spde$n.spde)
  # Projection matrix
  A = inla.spde.make.A(mesh = mesh, loc = coo)
  
  # Stack with data for estimation. Effects include intercept and covariates
  stk.e = inla.stack(tag = "est", data = list(y = d$y), A = list(1, A), effects = list(d[, covnames, drop = FALSE], s = indexs))
  if(is.null(dp)){
    stk.full = inla.stack(stk.e)
  }else{
    # Prediction coordinate locations and projection matrix
    coop = cbind(dp$coox, dp$cooy)
    Ap = inla.spde.make.A(mesh = mesh, loc = coop)
    # stack  
    stk.p = inla.stack(tag = "pred", data = list(y = NA), A = list(1, Ap), effects = list(dp[, covnames, drop = FALSE], s = indexs))
    stk.full = inla.stack(stk.e, stk.p)
  }
 # cres = list(return.marginals.predictor = TRUE, return.marginals.random = TRUE)
 # cinla = list(strategy = 'adaptive', int.strategy = 'eb')  #
  
  st1 = Sys.time()
  if(family == "gaussian"){
    # Formula that is specified in the arguments
    res = inla(formula, family = "gaussian", data = inla.stack.data(stk.full), 
               control.predictor = list(compute = TRUE, A = inla.stack.A(stk.full), link = 1),
               control.compute = list(config = TFPOSTERIORSAMPLES, return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE), 
               #control.results = list(return.marginals.predictor = TRUE, return.marginals.random = TRUE), 
               control.inla = list(strategy = 'adaptive', int.strategy = 'eb'),
               
               verbose=TRUE)}
  if(family == "Gamma"){
    cres = list(return.marginals.predictor = TRUE, return.marginals.random = TRUE)
    cinla = list(strategy = 'adaptive', int.strategy = 'eb')  #
    st1 = Sys.time()
    res = inla(formula, family = "Gamma", data = inla.stack.data(stk.full),
               control.predictor = list(compute = TRUE, A = inla.stack.A(stk.full), link = 1),
               control.compute = list(config = TFPOSTERIORSAMPLES, return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
               #control.results = cres, control.inla = cinla,
               verbose=TRUE)}
  if(family == "lognormal"){
    cres = list(return.marginals.predictor = TRUE, return.marginals.random = TRUE)
    cinla = list(strategy = 'adaptive', int.strategy = 'eb')  #
    st1 = Sys.time()
    res = inla(formula, family = "lognormal", data = inla.stack.data(stk.full),
               control.predictor = list(compute = TRUE, A = inla.stack.A(stk.full), link = 1),
               control.compute = list(config = TFPOSTERIORSAMPLES, return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
               control.results = cres, control.inla = cinla,
               verbose=TRUE)}
  st2 = Sys.time()
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
#' 

fnGetPredictions = function(res, stk.full, mesh, d, dp, covnames, NUMPOSTSAMPLES, cutoff_exceedanceprob){
  if(NUMPOSTSAMPLES == -1){
    index = inla.stack.index(stk.full, tag = "est")$data
    d$excprob = sapply(res$marginals.fitted.values[index],
                       FUN = function(marg){1-inla.pmarginal(q = cutoff_exceedanceprob, marginal = marg)})
    d$pred_mean = res$summary.fitted.values[index, "mean"]
    d$pred_ll = res$summary.fitted.values[index, "0.025quant"]
    d$pred_ul = res$summary.fitted.values[index, "0.975quant"]
    d$pred_sd = res$summary.fitted.values[index, "sd"]
    d$pred_ll90 = unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.05, marginal = marg)}))
    d$pred_ul90 = unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.95, marginal = marg)}))
    d$pred_ll50 = unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.25, marginal = marg)}))
    d$pred_ul50 = unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.75, marginal = marg)}))
    dres = d
  }
  if(NUMPOSTSAMPLES == 0){
    index = inla.stack.index(stk.full, tag = "pred")$data
    dp$excprob = sapply(res$marginals.fitted.values[index],
                        FUN = function(marg){1-inla.pmarginal(q = cutoff_exceedanceprob, marginal = marg)})
    dp$pred_mean = res$summary.fitted.values[index, "mean"]
    dp$pred_ll = res$summary.fitted.values[index, "0.025quant"]
    dp$pred_ul = res$summary.fitted.values[index, "0.975quant"]
    dp$pred_sd = res$summary.fitted.values[index, "sd"]
    dp$pred_ll90 = unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.05, marginal = marg)}))
    dp$pred_ul90 = unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.95, marginal = marg)}))
    dp$pred_ll50 = unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.25, marginal = marg)}))
    dp$pred_ul50 = unlist(lapply(res$marginals.fitted.values[index], FUN = function(marg){inla.qmarginal(p = 0.75, marginal = marg)}))
    dres = dp
  }
  
  if(NUMPOSTSAMPLES > 0){
    psamples = inla.posterior.sample(NUMPOSTSAMPLES, res)
    ps = sapply(1:NUMPOSTSAMPLES, fnPredictFromPosteriorSample, psamples = psamples, res = res, mesh = mesh, dp = dp, covnames = covnames)
    dp$excprob = apply(ps, 1, function(x){mean(x > cutoff_exceedanceprob)})
    dp$pred_mean = rowMeans(ps)
    dp$pred_ll = apply(ps, 1, function(x){quantile(x, 0.025)})
    dp$pred_ul = apply(ps, 1, function(x){quantile(x, 0.975)})
    dp$pred_sd = res$summary.fitted.values[index, "sd"]
    dp$pred_ll90 = apply(ps, 1, function(x){quantile(x, 0.05)})
    dp$pred_ul90 = apply(ps, 1, function(x){quantile(x, 0.95)})
    dp$pred_ll50 = apply(ps, 1, function(x){quantile(x, 0.25)})
    dp$pred_ul50 = apply(ps, 1, function(x){quantile(x, 0.75)})
    dres = dp
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
  xgb = h2o.xgboost(y = y_var , x= x_var, training_frame = merged_hex, nfolds =10, ntrees =1000, eta =  0.007, subsample = 0.7, max_depth = 6, gamma = 5, reg_lambda =2, reg_alpha = 0, keep_cross_validation_predictions = T, seed =1 )
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
  xgb_all = h2o.xgboost(y = y_var, x= x_var, training_frame = merged_hex, ntrees =1000, eta =  0.007, subsample = 0.7, max_depth = 6, gamma = 5, reg_lambda =2, reg_alpha = 0, seed =1 )
  xgb = as.data.frame(h2o.predict(object = xgb_all, newdata = new_hex))$predict 
  las_all = h2o.glm(y = y_var, x= x_var, training_frame = merged_hex, alpha = 1,  nfolds =1, seed=1 )
  lasso  = as.data.frame(h2o.predict(object = las_all, newdata = new_hex))$predict 
  h2o.removeAll()
  cbind(rf, xgb, lasso)
  
}


#' Calculates cross-validation measures obtained by fitting a spatial model using INLA and SPDE
#'
#' @param n Number of iteration 
#' @param d Data frame with data for estimation that contains coordinates (coox, cooy), response variable (y) and covariates
#' If \code{covnames} includes an intercept, \code{d} needs to have column of 1s for the intercept
#' @param covnames Vector with the names of the intercept and covariates to be included in the formula
#' @param spatialsample if doing spatial sampling for validation, i.e. sample less in dense areas but more in sparse areas so that the points sampled are even in space. This may however not be needed in air pollution mapping, where the ground stations are dense in places where it should be. 

#' @return Vector with the cross-validation results
#' 
#' 
INLA_crossvali =  function(n, d, dp, 
                           formula=as.formula(paste0('y ~ 0 + ', paste0(covnames, collapse = '+'), " + f(s, model = spde)")),
                           covnames, typecrossvali = "non-spatial", family = c("gaussian","Gamma","lognormal")){
  print(n)
  # Split data
  smp_size = floor(0.2 * nrow(d)) 
  set.seed(n)
  if(typecrossvali == "non-spatial"){
    test <- sample(seq_len(nrow(d)), size = smp_size)
    }
  if(typecrossvali == "spatial"){
    # The validation data needs to spatially represent the whole region where the prevalence is predicted
    # We use locations of a spatially representative sample of the prediction surface
    # To obtain a valid data set, X% of the observations are sampled without replacement where
    # each observation has a probability of selection proportional to the area of the Voronoi polygon
    # surrounding its location, that is, the area closest to the location relative to the surrounding points
    p <- matrix(coordi, ncol = 2)
    v <- dismo::voronoi(p) # extent?
    prob_selection <- area(v)/sum(area(v))
    test <- sample(seq_len(nrow(d)), size = smp_size, prob = prob_selection, replace = FALSE)
  }
  training = seq_len(nrow(d))[-test] 
  # Fit model
  dtraining = d[training, ]
  dptest = dp[test, ]
  
  lres = fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family =family)
    # Get predictions
  dptest = fnGetPredictions(lres[[1]], lres[[2]], lres[[3]], dtraining, dptest, covnames, NUMPOSTSAMPLES = 0, cutoff_exceedanceprob = 30)
  # Goodness of fit
  val = APMtools::error_matrix(validation = dptest$real, prediction = dptest$pred_mean)
  val = c(val, cor = cor(dptest$real, dptest$pred_mean))
  inlacrps = crps(y =dptest$real, family = "norm", mean = dptest$pred_mean, sd =dptest$pred_sd) 
  
  (val = c(val, covprob95 = mean(dptest$pred_ll <= dptest$real &  dptest$real <= dptest$pred_ul),  # 95% coverage probabilities
           covprob90 = mean(dptest$pred_ll90 <= dptest$real &  dptest$real <= dptest$pred_ul90),
           covprob50 = mean(dptest$pred_ll50 <= dptest$real &  dptest$real <= dptest$pred_ul50),
           meancrps = mean(inlacrps),
           mediancrps = median(inlacrps)))
  return(val)
} 

#' Calculates cross-validation measures obtained by fitting a spatial model using INLA and SPDE
#' It uses d and dp because datasets for estimation and prediction have different covariates (predictions by cross-validation and using all data)
#'
#' @param n Number of iteration 
#' @param d Data frame with data for estimation that contains coordinates (coox, cooy), response variable (y) and covariates
#' @param dp Data frame with data for prediction that contains coordinates (coox, cooy), and covariates
#' @param spatialsample if doing spatial sampling for validation, i.e. sample less in dense areas but more in sparse areas so that the points sampled are even in space. This may however not be needed in air pollution mapping, where the ground stations are dense in places where it should be. 
#' If \code{covnames} includes an intercept, \code{d} needs to have column of 1s for the intercept
#' @param covnames Vector with the names of the intercept and covariates to be included in the formula
#' @return Vector with the cross-validation results
INLA_stack_crossvali =  function(n, d, formula, covnames, spatialsample = F,family = "gaussian"){
  
  print(n)
  # Split data
  smp_size <- floor(0.2 * nrow(d)) 
  set.seed(n)
  
  
  if(!(spatialsample)){
    test <- sample(seq_len(nrow(d)), size = smp_size)
  }
  if(spatialsample){
    # The validation data needs to spatially represent the whole region where the prevalence is predicted
    # We use locations of a spatially representative sample of the prediction surface
    # To obtain a valid data set, X% of the observations are sampled without replacement where
    # each observation has a probability of selection proportional to the area of the Voronoi polygon
    # surrounding its location, that is, the area closest to the location relative to the surrounding points
    p <- matrix(c(d$coox, d$cooy), ncol = 2)
    v <- dismo::voronoi(p) # extent?
    prob_selection <- area(v)/sum(area(v))
    test<- sample(seq_len(nrow(d)), size = smp_size, prob = prob_selection, replace = FALSE)
  }
  
  
  training <- seq_len(nrow(d))[-test] 
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
  lres <-  fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family = family)
  # Get predictions
  dptest <- fnGetPredictions(lres[[1]], lres[[2]], lres[[3]], dtraining, dptest, covnames, NUMPOSTSAMPLES = 0, cutoff_exceedanceprob = 3)
 
  # Goodness of fit
  val <- APMtools::error_matrix(validation = dptest$real, prediction = dptest$pred_mean)
  val <- c(val, cor = cor(dptest$real, dptest$pred_mean))
  #print(val)
  inlacrps = crps(y =dptest$real, family = "norm", mean = dptest$pred_mean, sd =dptest$pred_sd) 
  
  val = c(val, covprob95 = mean(dptest$pred_ll <= dptest$real &  dptest$real <= dptest$pred_ul),  # 95% coverage probabilities
           covprob90 = mean(dptest$pred_ll90 <= dptest$real &  dptest$real <= dptest$pred_ul90),
           covprob50 = mean(dptest$pred_ll50 <= dptest$real &  dptest$real <= dptest$pred_ul50),
           meancrps = mean(inlacrps),
           mediancrps = median(inlacrps))
   
  
  return(val)
} 

############stack without inla ######################################
#' Stack with h2o, not using inla, return error metrics. 
#'  @param  d: dataframe with predictors and response
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
  xgb = h2o.xgboost(y = y_var , x= x_var, training_frame = merged_hex, nfolds =10, ntrees =1000, eta =  0.007, subsample = 0.7, max_depth = 6, gamma = 5, reg_lambda =2, reg_alpha = 0, keep_cross_validation_predictions = T, seed =1 )
  las = h2o.glm(y = y_var , x= x_var, training_frame = merged_hex, alpha = 1,  nfolds =10, keep_cross_validation_predictions = T, seed =1 )
  
  ens = h2o.stackedEnsemble(y = y_var, x= x_var, training_frame = merged_hex, base_models = c(rf, xgb, las))
  dtest <- fnMLPredictionsAll(d=d, training = training, test = test)
  dtesth2o = as.h2o(dtest)
  #print(dtest)
  pred = as.data.frame(h2o.predict(object =ens, newdata = dtesth2o))$predict 

  val = APMtools::error_matrix(validation = ytest, prediction = pred)
  
  h2o.removeAll()
  return(val)
}