# crossvalidation on server: R -e "rmarkdown::render('/data/lu01/Glo_crossvali.Rmd', output_file = '/data/lu01/glo_cossvali.html')"

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
