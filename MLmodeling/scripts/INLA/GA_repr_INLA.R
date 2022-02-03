rm(list = ls())
#remove.packages("INLA")
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#devtools::install_github("mengluchu/APMtools") 
library(INLA)
 
inla.setOption(inla.mode="experimental")

ipak <- function(pkg){
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE , repos='http://cran.muenster.r-project.org')
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("rgdal", "mapview","cowplot", "viridis","devtools", "dplyr","data.table" , "ggplot2" , "RColorBrewer", "tidyr" ,"tibble",  "sf", "dismo","scoringRules")
ipak(packages)

install_github("mengluchu/APMtools") 
library(APMtools)
#ls("package:APMtools") 
# INLA FUNCTIONS
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
               inla.mode = "experimental",
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
  if(family == "beta"){
    cres = list(return.marginals.predictor = TRUE, return.marginals.random = TRUE)
    cinla = list(strategy = 'adaptive', int.strategy = 'eb')  #
    st1 = Sys.time()
    res = inla(formula, family = "beta", data = inla.stack.data(stk.full),
               control.predictor = list(compute = TRUE, A = inla.stack.A(stk.full), link = 1),
               control.compute = list(config = TFPOSTERIORSAMPLES, return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
               control.results = cres, control.inla = cinla,
               verbose=TRUE)}
  st2 = Sys.time()
  print(st2-st1)
  return(list(res, stk.full, mesh, coo, coop))
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
INLA_crossvali =  function(n, d, dp, formula, covnames, typecrossvali = "non-spatial", family = "gaussian"){
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
  
  # INI CODE MENG
  # Add to d[, training] 3 variables with names lasso, rf, xgb that are predictions at locations coox and cooy calculated using cross-validation with the dataset d[, training]
  #dtraining = fnMLPredictionsCV(d[, training])
  # Add to dp[test, ] 3 variables with names lasso, rf, xgb that are predictions at locations coox and cooy calculated using all data in d[, training]
  #dptest = fnMLPredictionsAll(dp[, test])
  # END CODE MENG
  
  
  # Fit model
  if(family == "gaussian"){
    lres = fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family = "gaussian")
  }
  if(family == "Gamma"){
    lres = fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family = "Gamma")
  }
  if(family == "lognormal"){
    lres = fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family = "lognormal")}
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


INLA_cvsp =  function(n, d, dp, formula, covnames, typecrossvali = "non-spatial", family = "gaussian"){
  print(n)
  # Split data
  X<-split(d,  d$grp)
  dptest= X[[n]] 
  
  dtraining = anti_join(d,X[[1]])
  
  # Fit model
  if(family == "gaussian"){
    lres = fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family = "gaussian")
  }
  if(family == "Gamma"){
    lres = fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family = "Gamma")
  }
  if(family == "lognormal"){
    lres = fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family = "lognormal")}
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


# for block CV
plotspcv= function(result){
  
  VLAdf = as.data.frame(result) 
  names(VLAdf)= paste0("CV", 1:20)
  
  #data.frame(xgbcv,h2o.cross_validation_fold_assignment(xgb))
  #h2o.cross_validation_fold_assignment(xgb) 
  
  
  merged2 = na.omit(stc%>%dplyr::select(c(colnames(merged), grp)))
  m3 = distinct(merged2)
  plot(m3)
  m3$r2 = unlist(ifelse(VLAdf[7,]> 0.1, VLAdf["rsq",], 0.1))
  
  
  m3$medcrps = unlist(VLAdf['mediancrps',]) 
  
  
  mapviewOptions(
    basemaps = c("OpenStreetMap.Mapnik","Esri.OceanBasemap")
    , raster.palette = colorRampPalette(rev(brewer.pal(9, "PiYG")))
    , vector.palette = colorRampPalette(brewer.pal(9, "PuBuGn"))
    , na.color = "gray"
    , layers.control.pos = "topright"
    , viewer.suppress = TRUE # open browser
  )
  
  plot(m3['r2'])
  mapview(m3['r2'],layer.name= "R2",col.regions =rev(brewer.pal(11, "PiYG")), at = seq(0.1, 0.9, 0.1)) + mapview(locations_sf["urbantype_chara"],col.regions = c(brewer.pal(3, "Paired")[2],brewer.pal(3, "Accent")[1:2],brewer.pal(3, "Dark2")[1]), layer.name= "NO2", cex = 4)
  #mapshot(a, file = paste0(getwd(), "/R2map.png")) 
  
  mapview(m3['medcrps'],layer.name= "median_CRPS",col.regions =rev(brewer.pal(11, "PiYG")), at = seq(0, 5.5, 0.5)) + mapview(locations_sf["urbantype_chara"],col.regions = c(brewer.pal(3, "Paired")[2],brewer.pal(3, "Accent")[1:2],brewer.pal(3, "Dark2")[1]), layer.name= "NO2", cex = 4)
}

# data
#shapefile <- readOGR(dsn = "~/Documents/GitHub/uncertainty/shapefiles/", layer = "shape2")

 
####################
# INLA
####################
covnames = c("b0", "nightlight_450", "population_1000", "population_3000",
             "road_class_1_5000", "road_class_2_100", "road_class_3_300",  
             "trop_mean_filt", "road_class_1_100")

 
# data used
mergedall = read.csv("https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/DENL17_uc.csv")

load("~/Documents/GitHub/uncertainty/data_vis_exp/missingstation.rda") #msname
mergedall =mergedall%>%filter(!(AirQualityStation %in% msname)) #474


if (resolution ==100)
{
  mergedall = mergedall%>%dplyr::select(-c(industry_25,industry_50,road_class_1_25,road_class_1_50,road_class_2_25,road_class_2_50,   road_class_3_25,road_class_3_50))
}  

merged= mergedall%>%dplyr::select(matches(varstring))%>% na.omit() # there is actually no na in this file, but for now RF and LA doesnt deal with missing data, leave out for quick examination 

# ####################################
# merged is the final dataset to use
######################################


#======================================
# Data for prediction
#====================

d2= merged%>%dplyr::select(covnames)%>%scale()%>%data.frame
d2$b0 = 1 # intercept
d2$y = d$mean_value #     # For GAMMA distribution
d2$coox = d$Longitude
d2$cooy = d$Latitude

d2$real = d$y

#spatial and non-spatial model
formula = as.formula(paste0('y ~ 0 + ', paste0(covnames, collapse = '+'), " + f(s, model = spde)"))



lres <- fnFitModelINLA(d2, dp = d2, covnames, formula = formula, TFPOSTERIORSAMPLES = TRUE, family = "gaussian")
lres
  
#compare with nonspatial INLA model
formula2 = as.formula(paste0('y ~ 0 + ', paste0(covnames, collapse = '+')))
lres2 <- fnFitModelINLA(d2, dp = d2, covnames, formula = formula2, TFPOSTERIORSAMPLES = TRUE, family = "gaussian")
lres2

#==========================================
# spatial and non-spatial INLA model comparison, figure 3 and 4 supplementary
#============================================
lres[[1]]$summary.fixed
lres2[[1]]$summary.fixed

#### test function fnGetPredictions
#dtraining = d[1:400, ]
#dptest = dp[401:474, ]
#res = lres[[1]]
#fnGetPredictions(lres[[1]], lres[[2]], lres[[3]], dtraining, dptest, covnames, NUMPOSTSAMPLES = 0, cutoff_exceedanceprob = 30)
###


#===============================================
# model hyperparameter, figure 5 supplementary
#==================================================
par(mfrow=c(2,2))
# Intercept
plot(lres[[1]]$marginals.fix[[1]], type='l', xlab=expression(Intercept), ylab="Density", cex.lab=1.6, cex.axis=1.4)

# Precision for Gaussian observations
plot(lres[[1]]$marginals.hy[[1]], type='l',xlab=expression(tau[y]), ylab="Density", cex.lab=1.6, cex.axis=1.4)

# Marginal for tau of the spatial random field
plot(lres[[1]]$marginals.hy[[2]], type='l',xlab=expression(tau[x]), ylab="Density", cex.lab=1.6, cex.axis=1.4)

# Marginal for kappa of the spatial random field
plot(lres[[1]]$marginals.hy[[3]], type='l',xlab=expression(kappa), ylab="Density", cex.lab=1, cex.axis=1)


#=========================================
#     Get predicted data on grid
#=========================================

index.pred <- inla.stack.index(lres[[2]], "pred")$data

pred_mean <- lres[[1]]$summary.fitted.values[index.pred, "mean"]
pred_ll <- lres[[1]]$summary.fitted.values[index.pred, "0.025quant"]
pred_ul <- lres[[1]]$summary.fitted.values[index.pred, "0.975quant"]

#==============================
#     Plot of predictions
#==============================

data_pred = data.frame(pred_mean, pred_ll, pred_ul, lres[[5]][, 1], lres[[5]][, 2])
colnames(data_pred) <- c("pred_mean", "pred_ll", "pred_ul", "Longitude", "Latitude")



#==================================
# predictions, figure 9 manuscript
#=================================
resdf = data.frame(observation =d$y, prediction_mean = data_pred$pred_mean, prediction_upper = data_pred$pred_ul,prediction_lower = data_pred$pred_ll, lat =d$Latitude, lon = d$Longitude)

gres = gather(resdf, "key", "value",-lat, -lon)

ggplot(gres, aes(lon,lat)) + 
  ggtitle("Observation vs. INLA prediction") +
  theme(plot.title = element_text(hjust = 0, size = 10),strip.text.x = element_text(size = 15, colour = "black", angle = 0))+
  geom_point(aes(colour=value)) +facet_wrap(~key)+ 
  geom_polygon(data = shapefile, colour = "black",aes(x = long, y = lat, group = group), fill = NA) + 
  scale_color_gradientn(name = expression(paste(NO[2],~mu, g/m^{3})),
                        colours = c(viridis(100, begin = 0.3, end = 0.9),rev( magma(100, begin = 0.3))), limits = c(0,48))

ggsave("pred.png",height = 10, width = 10)


##============================================
## differences between prediction and observations, figure 2, supplementary
##=============================================

res_dif = data.frame(obs_predmean = d$y- data_pred$pred_mean,   obs_predul = d$y- data_pred$pred_ul, obs_predll = d$y- data_pred$pred_ll,  lat =d$Latitude, lon = d$Longitude)
gres = gather(res_dif, "key", "value",-lat, -lon)

ggplot(gres, aes(lon,lat)) + 
  ggtitle("Observation - INLA prediction") +
  theme(plot.title = element_text(hjust = 0, size = 10),strip.text.x = element_text(size = 15, colour = "black", angle = 0))+
  geom_point(aes(colour=value)) +facet_wrap(~key)+ 
  geom_polygon(data = shapefile, colour = "black",aes(x = long, y = lat, group = group), fill = NA) + 
  scale_color_gradientn(name = expression(paste(NO[2],~mu, g/m^{3})),
                        colours = rev(brewer.pal(10,"Spectral")))

ggsave("dif_obs_INLA.png",height = 6, width = 10)

#======================================================
#    Posterior mean and sd of the spatial random field
#======================================================

#==========================
#    figure 8 manuscript
#==========================
rang <- apply(lres[[3]]$loc[, c(1, 2)], 2, range)
proj <- inla.mesh.projector(lres[[3]], xlim = rang[, 1], ylim = rang[, 2], dims = c(300, 300))

mean_field <- inla.mesh.project(proj, lres[[1]]$summary.random$s$mean)
sd_field   <- inla.mesh.project(proj, lres[[1]]$summary.random$s$sd)


dat <- expand.grid(x = proj$x, y = proj$y)
dat$mean_field <- as.vector(mean_field)
dat$sd_field <- as.vector(sd_field)


gmean <- ggplot(dat, aes(x = x, y = y, fill = mean_field)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw() +geom_polygon(data = shapefile, aes(x = long, y = lat, group = group), colour = "black", fill = NA)

ggsave("~/Documents/GitHub/uncertainty/meanrf.png",height = 12, width = 8)

gsd <-   ggplot(dat, aes(x = x, y = y, fill = sd_field)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw() +geom_polygon(data = shapefile, aes(x = long, y = lat, group = group), colour = "black", fill = NA)

grid.arrange(gmean, gsd)  
#ggsave("~/Documents/GitHub/uncertainty/field_mean_sd.png")


#=====================
# non-spatial Cross-validation, INLA and INLA-G
#========================
VLAg = lapply(1:20, FUN = INLA_crossvali, d = d, dp = d, formula = formula, covnames = covnames, 
              typecrossvali = "non-spatial", family = "gaussian")

VLAma = lapply(1:20, FUN = INLA_crossvali, d = d, dp = d, formula = formula, covnames = covnames, 
               typecrossvali = "non-spatial", family = "Gamma")

(VLA = data.frame(LA = rowMeans(data.frame(VLAma))))
(VLAg = data.frame(LA = rowMeans(data.frame(VLAmag))))

#=====================================
#SPblock
#==================

#get the spatial points and make the grid
locations_sf = st_as_sf(mergedall, coords = c("Longitude","Latitude"), crs=4642)

locations_sf$Latitude = mergedall$Latitude
locations_sf$Longitude = mergedall$Longitude
grid1 = st_make_grid(locations_sf, 2) #"2" is the cell-size. 64 grids if 1 degree #20 grids
grid1%>%plot  

stc = st_join(st_sf(grid1), locations_sf,join=st_intersects )
stc$grp = sapply(st_equals(stc), max)
merged = mergedall%>%dplyr::select(matches(varstring))%>% na.omit() # there is actually no na in this file, but for now RF and LA doesnt deal with missing data, leave out for quick examination 

d = data.frame(na.omit(stc))%>%dplyr::select(c(colnames(merged), 'grp', "Longitude", "Latitude"))%>%
  mutate(grp = as.factor(grp))

d$y = d$mean_value #     # For GAMMA distribution
d$coox = d$Longitude
d$cooy = d$Latitude
d$b0 = 1 # intercept
d$real = d$y
 
n = d$grp%>%unique()%>%length() # number of grid cells

VLA = lapply(1:n, FUN = INLA_cvsp, d = d, dp = d, formula = formula, covnames = covnames,  typecrossvali = "non-spatial", family = "gaussian")

apply(data.frame(VLA),1, mean, na.rm=T)


plotspcv(VLA)


####
#SP2 constomed CV
#=====================

 
#- "tr_hp": close to roads, high population
#- "tr_mlp": close to roads, middle-low population 
#- "f": far away from roads 

 
INLA_crossvali2 =  function(n, test, training, d, dp, formula, covnames, typecrossvali = "non-spatial", family = "gaussian"){
  
  dtraining = d[training, ]
  dptest = dp[test, ]
  
  if(family == "gaussian"){
    lres = fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family = "gaussian")
  }
  if(family == "Gamma"){
    lres = fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family = "Gamma")
  }
  if(family == "lognormal"){
    lres = fnFitModelINLA(dtraining, dptest, formula, covnames, TFPOSTERIORSAMPLES = FALSE, family = "lognormal")}
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

 
# Variables for stacked generalization
# d$lasso = d$lasso10f_pre
# d$rf = d$rf10f_pre
# d$xgb = d$xgb10f_pre
d$Countrycode  = as.factor(mergedall$Countrycode)
d$MeasurementType  = as.factor(mergedall$MeasurementType)
d$AirQualityStationType = as.factor(mergedall$AirQualityStationType)
d$AirQualityStationArea = as.factor(mergedall$AirQualityStationArea)
d$urbantype = as.factor(mergedall$urbantype)
 
 
sp2_cv =  function(n, df_type= c("tr_hp", "tr_mlp", "f") , df_model, y_var) {
  
  set.seed(n)
  a= df_model%>%filter((road_class_2_100 > 0 | road_class_1_100 > 0|road_class_3_100>quantile(road_class_3_100, .75)) & population_1000 > quantile(population_1000, 0.75))  
  #traffic_lmpop 
  b= df_model%>%filter((road_class_2_100 > 0 | road_class_1_100 > 0 |road_class_3_100 > quantile(road_class_3_100, .75)) & population_1000 < quantile(population_1000, 0.5))   
  
  #fartr_highpop
  c= df_model%>%filter((road_class_2_100 == 0 & road_class_1_100 == 0 & road_class_3_100 < quantile(road_class_3_100, .5)))  
  print(nrow(a))
  print(nrow(b))
  print(nrow(c))
  #85 65 177
  
  methodID = switch(df_type,  "tr_hp"=1,"tr_mlp" =2,"f"=3  ) 
  totest = switch(methodID,a,b,c)
  
  
  others = setdiff(df_model, totest)
  orderedall=rbind(totest, others) #order data
  nrow(orderedall)
  
  #each time use 7% of the data satisfying with the conditions for testing. 
  #test_size = floor(0.07*nrow(df_model)) # 30  is about 20% of traffic, use a consistent size about 7% of data
  test_size = floor(0.07*nrow(df_model)) # 8  is about 6% of traffic, use a consistent size about 2% of data
  
  test = sample(nrow(totest), size = test_size) # sample 20% from e.g. traffic and then use others as training 
  training = setdiff(seq_len(nrow(df_model)), test)

  INLA_crossvali2(d = orderedall, dp =orderedall, formula = formula, covnames = covnames,   typecrossvali = "non-spatial", family = "gaussian",training=training, test=test)
 
}  

tr_hp = lapply(1:nboot, df_type ="tr_hp",  df_model =merged, y_var = y_var, sp2_cv)%>%data.frame() 
tr_lmp= lapply(1:nboot, df_type ="tr_mlp", df_model =merged, y_var = y_var, sp2_cv)%>%data.frame() 
far= lapply(1:nboot, df_type ="f", df_model =merged, y_var = y_var, sp2_cv)%>%data.frame() 

#tr_hp = lapply(1:nboot, df_type ="tr_hp",  df_model =d, y_var = y_var, sp2_cv)%>%data.frame() 
#tr_lmp= lapply(1:nboot, df_type ="tr_mlp", df_model =d, y_var = y_var, sp2_cv)%>%data.frame() 
#far= lapply(1:nboot, df_type ="f", df_model =d, y_var = y_var, sp2_cv)%>%data.frame() 

F1 = function(m, pre, f=quote(summary), nvaria) {apply(pre[, seq(m, ncol(pre), by =nvaria)], 1, f)}

nv = 1# number of algorithms.
cv_traffic= data.frame(sapply(1:nv, F1, tr_hp, mean,nv)) 
names(cv_traffic) = paste0(c("INLA"), "_tr_hp")

cv_bg = data.frame(sapply(1:nv, F1, tr_lmp, mean,nv)) 
names(cv_bg) =  paste0(c("INLA"),"_tr_lmp")

cv_far = data.frame(sapply(1:nv, F1, far, mean,nv)) 
names(cv_far) =  paste0(c( "INLA"),"_far")
cbind(cv_traffic, cv_bg, cv_far)
 

