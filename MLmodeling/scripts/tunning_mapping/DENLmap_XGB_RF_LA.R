# crossvalidation on server: R -e "rmarkdown::render('/data/lu01/Glo_crossvali.Rmd', output_file = '/data/lu01/glo_cossvali.html')"
xgbname = "/data/lu01/denl_xgb.tif"
rfname = "/data/lu01/denl_rf.tif"
laname = "/data/lu01/denl_la.tif"
resolution = 100
xgbname = "/data/lu01/mean_NLDE_xgb.tif" 
set.seed(1)
ipak <- function(pkg){
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos='http://cran.muenster.r-project.org')
  sapply(pkg, require, character.only = TRUE)
}
packages <- c( "raster","rgdal" ,"dplyr", "devtools",  "xgboost","data.table" , "randomForest", "glmnet" ) 
ipak(packages) 
install_github("mengluchu/APMtools")
library(APMtools)

# prepare data for mapping
glo = read.csv("/data/GLOBAL_LUR/glo4variables.csv")
 
y_var = "mean"
prestring = "road|nightlight|population|temp|wind|trop|indu|elev|radi"
varstring = paste(prestring,y_var,sep="|")
inde_var = glo %>%dplyr::select(matches(varstring))
if (resolution == 100){
  inde_var = inde_var%>%select(-c(industry_25,industry_50,road_class_1_25,road_class_1_50,road_class_2_25,road_class_2_50, road_class_3_25, road_class_3_50))
  }
# prepare rasters for mapping

lf = list.files("/data/gghdc/gap/output/2020_09_de_nl/100m/laea", pattern = paste0("(^",prestring,")",".*\\.map$"), full.names =T)
#lf = list.files("/data/gghdc/gap/output/2020_06_phoenutr/1/laea", pattern = paste0("(",prestring,")",".*\\.map$"), full.names =T)
#lf = list.files("/Volumes/t/2020_06_world/109624/laea", pattern = ".map$", full.names = T)
#lf = lf[-1]
sr = stack(lf)
#use this to merge roads if needed: sr[[names(sr)[grepl("road_class",names(sr))]]]
names(sr)
y_var = paste0(y_var, "_value") 

predicLA_RF_XGBtiles(df = inde_var, rasstack = sr, yname = y_var, varstring = prestring, xgbname=xgbname, rfname = rfname, laname = laname, ntree = 1000, mtry = 34,  
                     nrounds = 3000, eta = 0.007, gamma =5,max_depth = 6, xgb_alpha = 0, xgb_lambda = 2, subsample=0.7 )

###########################################################################
# visualize locally
###############################################################
rf1= raster("~/Downloads/non_backups/rfmph.tif")
#xgb1= raster("~/Downloads/non_backups/xgbph.tif")
#la1= raster("~/Downloads/non_backups/Lamnph.tif")
xgb2= raster("~/Downloads/xgbph_slowlr.tif")
xgb3= raster("~/Downloads/lr003lam002.tif")
xgb_lr003_alp1= raster("~/Downloads/densexgb.tif")
xgb_lr003_alp0= raster("~/Downloads/densexgb_alpha0.tif")
xgb_lr003_alp10= raster("~/Downloads/densexgb_alpha10.tif")

names(xgb_lr003_alp0)= "alpha0"
names(xgb_lr003_alp1)= "alpha1"
library(tmap)
library(sf)
xgb1 = raster("~/Downloads/phoen_xgb_py.map")
myTheme = rasterTheme(region = rev(brewer.pal(7, "Spectral")))
levelplot(rf1, at =seq(5, 50, by =3) , par.settings = myTheme)
levelplot(xgb1, par.settings = myTheme) # 
#quantile(global_annual$value_mean, c(0.01, 0.999, 0.99999))
#quite clear LAsso is not robust against exterme values. Question is do we need to remove very high values, like larger than 200, before aggregation? the maximum is 98 for the meanIt's going to be 62 for 99% and 84 for 99.9%
levelplot(xgb2,   par.settings = myTheme)
library(APMtools)
data(global_annual)
locations_sf = st_as_sf(global_annual, coords = c("long","lat"))

rf1 = projectRaster(rf1, crs=   CRS("+init=epsg:4326"))
xgb2 = projectRaster(xgb2, crs=   CRS("+init=epsg:4326"))
xgb3 = projectRaster(xgb3, crs=   CRS("+init=epsg:4326"))
xgb_lr003_alp1=projectRaster(xgb_lr003_alp1, crs=   CRS("+init=epsg:4326"))
xgb_lr003_alp0=projectRaster(xgb_lr003_alp0, crs=   CRS("+init=epsg:4326"))
xgb_lr003_alp10=projectRaster(xgb_lr003_alp10, crs=   CRS("+init=epsg:4326"))


osm_valuemean = tm_shape(rf1)+
  tm_raster(names(rf1),palette = "-Spectral", n = 9,alpha = 0.9)+
  tm_shape(xgb2)+
  tm_raster(names(xgb2),palette = "-Spectral", n = 9,alpha = 0.9)+
  tm_shape(xgb3)+
  tm_raster(names(xgb3),palette = "-Spectral", n = 9,alpha = 0.9)+
  tm_shape(xgb_lr003_alp1)+
  tm_raster(names(xgb_lr003_alp1),palette = "-Spectral", n = 9,alpha = 0.9)+
  tm_shape(locations_sf)+
  tm_dots( "value_mean", col = "value_mean", size = 0.05,title = "NO2 value",
           popup.vars = c("value_mean" ))+
  tm_view(basemaps = c('OpenStreetMap'))
tmap_save(osm_valuemean, "~/Downloads/rfxgb_phe.html")


osm_valuemean = tm_shape(xgb_lr003_alp0)+
  tm_raster(names(xgb_lr003_alp0),palette = "-Spectral", n = 9,alpha = 0.9, breaks=c(0, 5,10,15,20,25,30,35,40,45,50, Inf)) +tm_layout(legend.show = F)+
  tm_shape(xgb_lr003_alp1)+
  tm_raster(names(xgb_lr003_alp1),palette = "-Spectral", n = 9,alpha = 0.9,breaks=c(0, 5,10,15,20,25,30,35,40,45,50, Inf))+tm_layout(legend.show = F)+
  tm_shape(xgb_lr003_alp10)+tm_raster(names(xgb_lr003_alp10),palette = "-Spectral", n = 9,alpha = 0.9,breaks=c(0, 5,10,15,20,25,30,35,40,45,50, Inf))+tm_layout(legend.show = T,main.title = "NO2")+
  tm_shape(locations_sf)+ 
  tm_dots( "value_mean", col = "value_mean", size = 0.05,title = "NO2 value",
           popup.vars = c("value_mean" )) 
tmap_save(osm_valuemean, "~/Downloads/alpha0110.html")

hist(global_annual$value_mean)
quantile(global_annual$value_mean, c(0.01, 0.999))

