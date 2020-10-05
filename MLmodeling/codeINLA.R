source("INLA_util")

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
  smp_size <- floor(0.2 * nrow(d)) 
  set.seed(n)
  if(typecrossvali == "crossvalinotspatial"){
    test <- sample(seq_len(nrow(d)), size = smp_size)
  }
  if(typecrossvali == "crossvalispatial"){
    # The validation data needs to spatially represent the whole region where the prevalence is predicted
    # We use locations of a spatially representative sample of the prediction surface
    # To obtain a valid data set, X% of the observations are sampled without replacement where
    # each observation has a probability of selection proportional to the area of the Voronoi polygon
    # surrounding its location, that is, the area closest to the location relative to the surrounding points
    p <- matrix(c(d$coox, d$cooy), ncol = 2)
    v <- dismo::voronoi(p) # extent?
    prob_selection <- area(v)/sum(area(v))
    test <- sample(seq_len(nrow(d)), size = smp_size, prob = prob_selection, replace = FALSE)
  }
  
  training <- seq_len(nrow(d))[-test] 
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
  
  (val <- c(val, covprob95 = mean(dptest$pred_ll <= dptest$real &  dptest$real <= dptest$pred_ul),  # 95% coverage probabilities
            covprob90 = mean(dptest$pred_ll90 <= dptest$real &  dptest$real <= dptest$pred_ul90),
            covprob50 = mean(dptest$pred_ll50 <= dptest$real &  dptest$real <= dptest$pred_ul50)))
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