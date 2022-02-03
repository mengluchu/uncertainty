ipak <- function(pkg){
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE , repos='http://cran.muenster.r-project.org')
  sapply(pkg, require, character.only = TRUE)
}
packages <- c( "devtools", "dplyr","data.table" , "ggplot2" , "RColorBrewer", "xgboost",  "glmnet", "ranger", "randomForest","tidyr" ,"tibble","stargazer", "sf", "CAST", "caret", "quantregForest", "pdp", "h2o", "dismo","scoringRules")
ipak(packages)
install_github("mengluchu/APMtools") 
library(APMtools)
ls("package:APMtools") 
#source("https://raw.githubusercontent.com/mengluchu/uncertainty/master/MLmodeling//INLA/INLA_util.R")
source("~/Documents/GitHub/uncertainty/MLmodeling/scripts/INLA/INLA_util.R")
mergedall = read.csv("https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/DENL17_uc.csv")
library("mapview")
library("RColorBrewer")

resolution =100 # resolution of grid
nboot = 20  # number of bootstraps
y_var = "mean_value"
prestring =  "road|nightlight|population|temp|wind|trop|indu|elev|radi"
varstring = paste(prestring,y_var,sep="|")

if (resolution ==100)
{
  mergedall = mergedall%>%dplyr::select(-c(industry_25,industry_50,road_class_1_25,road_class_1_50,road_class_2_25,road_class_2_50,   road_class_3_25,road_class_3_50))
}  

merged = mergedall%>%dplyr::select(matches(varstring))%>% na.omit() # there is actually no na in this file, but for now RF and LA doesnt deal with missing data, leave out for quick examination 
table(mergedall$urbantype_chara)


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

rf_la_spcv =  function(n,d, y_var,typecrossvali = "crossvalinotspatial", coordi = c(mergedall$Longitude, mergedall$Latitude)) {
   
  d$id = 1 :nrow(d) # assign a new id because the rowname does not match the index. because rf_lasso_LUR takes indices for traing and testing.
  X<-split(d,  d$grp)
  dptest= X[[n]] # test id 
  test = dptest$id 
  training = setdiff(d$id, test) 
  #sum(is.na(d[training,]))
  rf_Lasso_LUR(d, numtrees =  2000, mtry = NULL,  y_varname= y_var, training=training, test=test, grepstring =prestring)
} 


xgb_spcv =  function(n,d, y_var,typecrossvali = "crossvalinotspatial", coordi = c(mergedall$Longitude, mergedall$Latitude)) {
  
  d$id = 1 :nrow(d) # assign a new id because the rowname does not match the index. because rf_lasso_LUR takes indices for traing and testing.
  X<-split(d,  d$grp)
  dptest= X[[n]] # test id 
  test = dptest$id 
  training = setdiff(d$id, test) 
  xgboost_LUR(d, y_varname= y_var, training=training, test=test, grepstring =varstring,  max_depth =6, gamma=5, eta =0.007, nrounds =1000, xgb_lambda = 2, xgb_alpha = 0, subsample = 0.7 )
  
} 

qrf_spcv =  function(n,d, y_var,typecrossvali = "crossvalinotspatial", coordi = c(mergedall$Longitude, mergedall$Latitude)) {
  
  d$id = 1 :nrow(d) # assign a new id because the rowname does not match the index. because rf_lasso_LUR takes indices for traing and testing.
  X<-split(d,  d$grp)
  dptest= X[[n]] # test id 
  test = dptest$id 
  training = setdiff(d$id, test) 
  q_rf_LUR(d, numtrees =  2000, mtry = NULL,  y_varname= y_var, training=training, test=test, grepstring =varstring)} 

} 

 

# data
locations_sf = st_as_sf(mergedall, coords = c("Longitude","Latitude"), crs=4642)

locations_sf$Latitude = mergedall$Latitude
locations_sf$Longitude = mergedall$Longitude
grid1 = st_make_grid(locations_sf, 2) # 64 grids #20 grids
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
# Variables for stacked generalization
# d$lasso = d$lasso10f_pre
# d$rf = d$rf10f_pre
# d$xgb = d$xgb10f_pre
names(d)


covnames = c("b0", "nightlight_450", "population_1000", "population_3000", 
             "road_class_1_5000", "road_class_2_100", "road_class_3_300",  
             "trop_mean_filt", "road_class_1_100")

#, "Countrycode",  "urbantype"

formula = as.formula(paste0('y ~ 0 + ', paste0(covnames, collapse = '+'), " + f(s, model = spde)"))

#====
# Cross-validation
n = d$grp%>%unique()%>%length()

VLA = lapply(1:20, FUN = INLA_cvsp, d = d, dp = d, formula = formula, covnames = covnames,  typecrossvali = "non-spatial", family = "gaussian")
rfla= lapply(1:20, d =d, y_var = y_var,rf_la_spcv)
xgb= lapply(1:20, d =d, y_var = y_var,xgb_spcv)
qrf= lapply(1:20, d =d, y_var = y_var,qrf_spcv)

#rfla = data.frame(RFLA = rowMeans(data.frame(rfla)))
 

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
apply(data.frame(xgb),1, mean, na.rm=T)
apply(data.frame(rfla),1, mean, na.rm=T)
apply(data.frame(VLA),1, mean, na.rm=T)


plotspcv(VLA)
plotspcv(rfla)


 