rm(list = ls())

#devtools::install_github("mengluchu/APMtools") 
library(INLA)
library(APMtools)
library(dismo) # for area() inside INLA_crossvali()
library(caret)
library("quantregForest")
#' Creates triangulated mesh to fit a spatial model using INLA and SPDE
#' 
#' @param coo coordinates to create the mesh
#' @return mesh boject
fnConstructMesh = function(coo){
  # meshbuilder()
  # offset: size of the inner and outer extensions around the data locations
  (offset1 = 1/8*max(dist(coo)))
  (offset2 = 1/8*max(dist(coo)))
  # max.edge: maximum allowed triangle edge lengths in the region and in the extension
  (maxedge1 = 1/30*max(dist(coo)))
  (maxedge2 = 1/5*max(dist(coo)))
  # cutoff: minimum allowed distance between points used to avoid building many small triangles around clustered locations
  (cutoff = 1/10000*max(dist(coo)))
  mesh = inla.mesh.2d(loc = coo, offset = c(offset1, offset2), cutoff = cutoff, max.edge = c(maxedge1, maxedge2))
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
fnFitModelINLA = function(d, dp, formula, covnames, TFPOSTERIORSAMPLES, family = ""){
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
  cres = list(return.marginals.predictor = TRUE, return.marginals.random = TRUE)
  cinla = list(strategy = 'adaptive', int.strategy = 'eb')  #
  st1 = Sys.time()
  if(family == "gaussian"){
  # Formula that is specified in the arguments
  res = inla(formula, family = "gaussian", data = inla.stack.data(stk.full), 
              control.predictor = list(compute = TRUE, A = inla.stack.A(stk.full), link = 1),
              control.compute = list(config = TFPOSTERIORSAMPLES, return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE), 
              control.results = cres, control.inla = cinla,
              verbose=TRUE)}
  if(family == "Gamma"){
    cres = list(return.marginals.predictor = TRUE, return.marginals.random = TRUE)
    cinla = list(strategy = 'adaptive', int.strategy = 'eb')  #
    st1 = Sys.time()
    res = inla(formula, family = "Gamma", data = inla.stack.data(stk.full),
               control.predictor = list(compute = TRUE, A = inla.stack.A(stk.full), link = 1),
               control.compute = list(config = TFPOSTERIORSAMPLES, return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
               control.results = cres, control.inla = cinla,
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
fnPredictFromPosteriorSample = function(psamples, ite, res, mesh, dp, covnames){
  # Retrieve elements sample
  (contents = res$misc$configs$contents)
  # betas for elements of covnames. covnames[1] is the first covariate (b0)
  id_firstcov = grep(covnames[1], rownames(psamples[[ite]]$latent))
  betas = psamples[[ite]]$latent[id_firstcov : (id_firstcov + (length(covnames)-1)), ]
  # spatial field
  id_s = which(contents$tag == "s")
  id_s = contents$start[id_s]:(contents$start[id_s] + contents$length[id_s] - 1)
  spatialfield = psamples[[ite]]$latent[id_s]
  # spat = lapply(ps, function(x) x$latent[id_s])
  # spat = matrix(unlist(spat), ncol = length(id_s), byrow = T)
  # Multiply model matrix times betas + spatial effect
  coop = cbind(dp$coox, dp$cooy)
  Ap = inla.spde.make.A(mesh = mesh, loc = coop)
  predictions = as.matrix(dp[, covnames]) %*% betas + drop(Ap %*% spatialfield)
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
fnGetPredictions = function(res, stk.full, mesh, d, dp, covnames, NUMPOSTSAMPLES, cutoff_exceedanceprob){
  if(NUMPOSTSAMPLES == -1){
    index = inla.stack.index(stk.full, tag = "est")$data
    d$excprob = sapply(res$marginals.fitted.values[index],
                        FUN = function(marg){1-inla.pmarginal(q = cutoff_exceedanceprob, marginal = marg)})
    d$pred_mean = res$summary.fitted.values[index, "mean"]
    d$pred_ll = res$summary.fitted.values[index, "0.025quant"]
    d$pred_ul = res$summary.fitted.values[index, "0.975quant"]
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
    dp$pred_ll90 = apply(ps, 1, function(x){quantile(x, 0.05)})
    dp$pred_ul90 = apply(ps, 1, function(x){quantile(x, 0.95)})
    dp$pred_ll50 = apply(ps, 1, function(x){quantile(x, 0.25)})
    dp$pred_ul50 = apply(ps, 1, function(x){quantile(x, 0.75)})
    dres = dp
  }
  return(dres) 
}



#' Calculates cross-validation measures obtained by fitting a spatial model using INLA and SPDE
#' It uses d and dp because datasets for estimation and prediction have different covariates (predictions by cross-validation and using all data)
#'
#' @param n Number of iteration 
#' @param d Data frame with data for estimation that contains coordinates (coox, cooy), response variable (y) and covariates
#' @param dp Data frame with data for prediction that contains coordinates (coox, cooy), and covariates
#' If \code{covnames} includes an intercept, \code{d} needs to have column of 1s for the intercept
#' @param covnames Vector with the names of the intercept and covariates to be included in the formula
#' @param typecrossvali string that denotes if cross-validation is spatial ("crossvalispatial") or not ("crossvalinotspatial")
#' @return Vector with the cross-validation results
INLA_crossvali =  function(n, d, dp, formula, covnames, typecrossvali = "non-spatial", family = ""){
  print(n)
  # Split data
  smp_size = floor(0.8 * nrow(d)) 
  set.seed(n)
  if(typecrossvali == "non-spatial"){
    training = sample(seq_len(nrow(d)), size = smp_size)
    test = seq_len(nrow(d))[-training] 
  }
  if(typecrossvali == "spatial"){
    # The validation data needs to spatially represent the whole region where the prevalence is predicted
    # We use locations of a spatially representative sample of the prediction surface
    # To obtain a valid data set, X% of the observations are sampled without replacement where
    # each observation has a probability of selection proportional to the area of the Voronoi polygon
    # surrounding its location, that is, the area closest to the location relative to the surrounding points
    p <- matrix(c(d$coox, d$cooy), ncol = 2)
    v <- dismo::voronoi(p) # extent?
    prob_selection <- area(v)/sum(area(v))
    test <- sample(seq_len(nrow(d)), size = nrow(d)-smp_size, prob = prob_selection, replace = TRUE)
    train = sample(seq_len(nrow(d)), size = smp_size)
    training = seq_len(nrow(d))[-train] 
  }
  # Fit model
  dtraining = d[training, ]
  dptest = dp[test, ]
     
   
  # Fit model
  if(family == "gaussian"){
  lres = fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family = "gaussian")}
  if(family == "Gamma"){
    lres = fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family = "Gamma")}
  if(family == "lognormal"){
    lres = fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family = "lognormal")}
  # Get predictions
  dptest = fnGetPredictions(lres[[1]], lres[[2]], lres[[3]], dtraining, dptest, covnames, NUMPOSTSAMPLES = 0, cutoff_exceedanceprob = 30)
  # Goodness of fit
  val = APMtools::error_matrix(validation = dptest$real, prediction = dptest$pred_mean)
  val = c(val, cor = cor(dptest$real, dptest$pred_mean))
  (val = c(val, covprob95 = mean(dptest$pred_ll <= dptest$real &  dptest$real <= dptest$pred_ul),  # 95% coverage probabilities
            covprob90 = mean(dptest$pred_ll90 <= dptest$real &  dptest$real <= dptest$pred_ul90),
            covprob50 = mean(dptest$pred_ll50 <= dptest$real &  dptest$real <= dptest$pred_ul50)))
  
    
  return(val)
} 




##################################################


#setwd("C:/Users/Usuario/Desktop/Projects/2021/KAUST/INLA")
d = read.csv("https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/DENL17_uc.csv")

#d = read.csv("dat2.csv", header = T)
head(d, 3)




#=======================================
# Data for estimation. Create variables y with the response, coox and cooy with the coordinates, and b0 with the intercept (vector of 1s)
d$y = sqrt(d$mean_value) # response transform sqrt (GAUSSIAN distribution)
#d$y = d$mean_value #     # For GAMMA distribution
d$coox = d$Longitude
d$cooy = d$Latitude
d$b0 = 1 # intercept
d$real = d$y
# Variables for stacked generalization
# d$lasso = d$lasso10f_pre
# d$rf = d$rf10f_pre
# d$xgb = d$xgb10f_pre
d$urbantype = as.factor(d$urbantype)

#======================================
# Data for prediction
dp = d
# Variables for stacked generalization
# dp$lasso = dp$lasso_all_pre
# dp$rf = dp$rf_all_pre
# dp$xgb = dp$xgb_all_pre




# MODEL 1
# Model with covariates selected with lasso
covnames = c("b0", "nightlight_450", "population_1000", "population_3000",
              "road_class_1_5000", "road_class_2_100", "road_class_3_300", "trop_mean_filt",
              "road_class_3_3000", "road_class_1_100", "road_class_3_100",
              "road_class_3_5000", "road_class_1_300", "road_class_1_500",
              "road_class_2_1000", "nightlight_3150", "road_class_2_300", "road_class_3_1000", 
             "temperature_2m_7","urbantype")

# # MODEL 2 10 from lasso and add urbantype
#covnames = c("b0", "nightlight_450", "population_1000", "population_3000",
#             "road_class_1_5000", "road_class_2_100", "road_class_3_300", "trop_mean_filt",
#             "road_class_3_3000", "road_class_1_100", "road_class_3_100",
#             "urbantype")

covnames = c("b0", "nightlight_450", "population_1000", "population_3000", 
             "road_class_1_5000", "road_class_2_100", "road_class_3_300",  
             "trop_mean_filt", "road_class_1_100", "urbantype")
#"road_class_3_5000", "road_class_1_300", "road_class_1_500",
#"road_class_2_1000", "nightlight_3150", "road_class_2_300", "road_class_3_1000",
#"temperature_2m_7", 





formula = as.formula(paste0('y ~ 0 + ', paste0(covnames, collapse = '+'), " + f(s, model = spde)"))



# Call inla()
lres = fnFitModelINLA(d, dp = dp, covnames, formula = formula, TFPOSTERIORSAMPLES = TRUE, family = "gaussian")

# Metrics (for comparative purposes)
lres[[1]]$waic$waic
sum(lres[[1]]$cpo$failure, na.rm = TRUE)

# # If cpo sum is > 1
# improved.result = inla.cpo(lres[[1]])
# 
# # Cheking
# sum(improved.result$cpo$failure, na.rm = TRUE)



slcpo <- function(m, na.rm = TRUE) {
  - sum(log(m$cpo$cpo), na.rm = na.rm)
}
slcpo(lres[[1]])


res = lres[[1]]
#res = improved.result
stk.full = lres[[2]]
mesh = lres[[3]]

# Get predictions. NUMPOSTSAMPLES = -1 calculated with estimation data, 0 with prediction data, 1 with inla.posterior.samples()
dres = fnGetPredictions(res, stk.full, mesh, d, dp, covnames, NUMPOSTSAMPLES = 0, cutoff_exceedanceprob = 30)

# Goodness of fit
APMtools::error_matrix(validation = dres$real, prediction = dres$pred_mean)
cor(dres$real, dres$pred_mean)
mean(dres$pred_ll <= dres$real &  dres$real <= dres$pred_ul)
mean(dres$pred_ll90 <= dres$real &  dres$real <= dres$pred_ul90)
mean(dres$pred_ll50 <= dres$real &  dres$real <= dres$pred_ul50)



# Cross-validation
VLA = lapply(1:20, FUN = INLA_crossvali, d = d, dp = dp, formula = formula, covnames = covnames, 
             typecrossvali = "non-spatial", family = "gaussian")

(VLA = data.frame(LA = rowMeans(data.frame(VLA))))
names(dp)

# Posterior predictive distribution
N = length(d$mean_value)
x = 1:N

sim = inla.posterior.sample(N, res)

s = sim[[1]]$latent
tail(s[,1], n=3)

plot(s[,1], typ = "l")
str(sim)



library(ranger)
library(dplyr)
n =1
y_var = "mean_value"
prestring =  "road|nightlight|population|temp|wind|trop|indu|elev|radi"
varstring = paste(prestring,y_var,sep="|")

for (n in 1:2){
  covnames0 <-covnames  
  
  set.seed(n)
  smp_size <- floor(0.8 * nrow(d)) 
  
  training <- sample(seq_len(nrow(d)), size = smp_size)
  test = seq_len(nrow(d))[-training] 
  
  y_denl = d[,y_var]
  y_denl_test = y_denl[test] 
  x_p = d%>%dplyr::select(matches(varstring))%>%dplyr::select(-y_var)
  
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
  
  d <- d[, c("mean_value", "Longitude", "Latitude", covnames0)]
  
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
  lres <- fnFitModelINLA(d, dp = dp, covnames, formula = formula, TFPOSTERIORSAMPLES = TRUE, family = "gaussian")

  # Get predictions
  dptest <- fnGetPredictions(lres[[1]], lres[[2]], lres[[3]], dtraining, dptest, covnames, NUMPOSTSAMPLES = 0, cutoff_exceedanceprob = 30)
  coordi = c(d$Longitude, d$Latitude)
  p <- matrix(coordi, ncol = 2)
     
  inla_90= cbind(dptest$pred_ll90,dptest$pred_ul90)
  rf_90 = t.quant90
  plot(inla_90[,1], ylim = c(-30,65), col = "red", typ = "l")
  lines(rf_90[,1])
  lines(inla_90[,2], col = "red")
  lines(rf_90[,2])}