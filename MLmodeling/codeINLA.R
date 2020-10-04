#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)
#devtools::install_github("mengluchu/APMtools") 
library(APMtools)

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
fnFitModelINLA <- function(d, dp, covnames, TFPOSTERIORSAMPLES, formulanew = NULL){
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
  
  # Formula
  # covnames includes the intercept and the covariates
  if(is.null(formulanew)){
  formula <- as.formula(paste0('y ~ 0 + ', paste0(covnames, collapse = '+'), " + f(s, model = spde)"))
  }
  else{
    formula <- formulanew
  }
  
  # Call inla(). Add config = TRUE if we want to call inla.posterior.sample()
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
#' \code{pred_ll} and \code{pred_ul} are the lower and upper limits of 95% credible intervals
#' \code{excprob} is the probability that hte prediction > cutoff value
fnGetPredictions <- function(res, stk.full, mesh, d, dp, covnames, NUMPOSTSAMPLES, cutoff_exceedanceprob){
  if(NUMPOSTSAMPLES == -1){
    index <- inla.stack.index(stk.full, tag = "est")$data
    d$pred_mean <- res$summary.fitted.values[index, "mean"]
    d$pred_ll <- res$summary.fitted.values[index, "0.025quant"]
    d$pred_ul <- res$summary.fitted.values[index, "0.975quant"]
    d$excprob <- sapply(res$marginals.fitted.values[index],
                         FUN = function(marg){1-inla.pmarginal(q = cutoff_exceedanceprob, marginal = marg)})
    dres <- d
  }
  if(NUMPOSTSAMPLES == 0){
    index <- inla.stack.index(stk.full, tag = "pred")$data
    dp$pred_mean <- res$summary.fitted.values[index, "mean"]
    dp$pred_ll <- res$summary.fitted.values[index, "0.025quant"]
    dp$pred_ul <- res$summary.fitted.values[index, "0.975quant"]
    dp$excprob <- sapply(res$marginals.fitted.values[index],
                         FUN = function(marg){1-inla.pmarginal(q = cutoff_exceedanceprob, marginal = marg)})
    dres <- dp
  }
  
  if(NUMPOSTSAMPLES > 0){
    psamples <- inla.posterior.sample(NUMPOSTSAMPLES, res)
    ps <- sapply(1:NUMPOSTSAMPLES, fnPredictFromPosteriorSample, psamples = psamples, res = res, mesh = mesh, dp = dp, covnames = covnames)
    dp$pred_mean <- rowMeans(ps)
    dp$pred_ll <- apply(ps, 1, function(x){quantile(x, 0.025)})
    dp$pred_ul <- apply(ps, 1, function(x){quantile(x, 0.975)})
    dp$excprob <- apply(ps, 1, function(x){mean(x > cutoff_exceedanceprob)})
    dres <- dp
  }
  return(dres) 
}



#' Calculates cross-validation measures obtained by fitting a spatial model using INLA and SPDE
#'
#' @param n Number of iteration 
#' @param d Data frame with data for estimation that contains coordinates (coox, cooy), response variable (y) and covariates
#' If \code{covnames} includes an intercept, \code{d} needs to have column of 1s for the intercept
#' @param covnames Vector with the names of the intercept and covariates to be included in the formula
#' @return Vector with the cross-validation results
INLA_crossvali =  function(n, d, covnames){
  print(n)
  # Split data
  smp_size <- floor(0.8 * nrow(d)) 
  set.seed(n)
  training <- sample(seq_len(nrow(d)), size = smp_size)
  test <- seq_len(nrow(d))[-training] 
  # Fit model
  # Data for prediction
  dp <- d
  dtraining <- d[training, ]
  dptest <- dp[test, ]
  # Fit model
  lres <- fnFitModelINLA(dtraining, dptest, covnames, TFPOSTERIORSAMPLES = FALSE, formulanew = NULL)
  # Get predictions
  dptest <- fnGetPredictions(lres[[1]], lres[[2]], lres[[3]], dtraining, dptest, covnames, NUMPOSTSAMPLES = 0, cutoff_exceedanceprob = 30)
  # Goodness of fit
  val <- APMtools::error_matrix(validation = dptest$real, prediction = dptest$pred_mean)
  val <- c(val, cor = cor(dptest$real, dptest$pred_mean))
  (val <- c(val, covprob = mean(dptest$pred_ll <= dptest$real &  dptest$real <= dptest$pred_ul))) # 95% coverage probabilities
  return(val)
} 

##################################################


 
d <- read.csv("~/Documents/GitHub/uncertainty/data_vis_exp/DENL17_uc.csv")

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

# Data for prediction
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

# Cross-validation
VLA <- lapply(1:20, FUN = INLA_crossvali, d = d, covnames = covnames)
(VLA <- data.frame(LA = rowMeans(data.frame(VLA))))

 
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