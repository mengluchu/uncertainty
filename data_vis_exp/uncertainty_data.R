library(dplyr)
library("mapview")
# GDAL version >= 3.1.0 | setting mapviewOptions(fgb = TRUE)
library("leafem")
library("leafgl")
library("leaflet")
library("sf")
# Linking to GEOS 3.8.0, GDAL 3.1.0, PROJ 7.0.0
library("stars")
library("raster")
library("RColorBrewer")

mapviewOptions(
  basemaps = c("OpenStreetMap.Mapnik","Esri.OceanBasemap")
  , raster.palette = colorRampPalette(rev(brewer.pal(9, "Greys")))
  , vector.palette = colorRampPalette(brewer.pal(9, "YlGnBu"))
  , na.color = "magenta"
  , layers.control.pos = "topright"
  , viewer.suppress = TRUE # open browser
)

merged = read.csv("~/Documents/GitHub/Global mapping/glo4variables.csv")
EUmeta = read.csv("~/Downloads/non_backups/global_mapping_data/PanEuropean_metadata.csv")
merged = read.csv("~/Documents/GitHub/Global mapping/glo_hr.csv")

selectcountry=function (merged, EUmeta)
{
EU = merge(merged, EUmeta, by=c("Longitude", "Latitude"), all.x =T)
EU = EU%>%filter(!is.na(Countrycode))
EU%>%filter(Countrycode=="DE"|Countrycode=="NL")%>%select(-X) # with X (index) we cant do spread
}
#EU 
#locations_sf = st_as_sf(EU, coords = c("Longitude","Latitude"), crs=4642)
 # https://leaflet-extras.github.io/leaflet-providers/preview/
#mapview(locations_sf)

#NL,DE


DENL_2017 = selectcountry(merged, EUmeta )
DENL_2017spread = spread(DENL_2017, hours, wkd_hr_value) 

 
write.csv(DENL_2017spread,"/Users/menglu/Documents/GitHub/mobiair/DENL17_hr_spread.csv")



#write.csv(DENL_2017,"/Users/menglu/Documents/GitHub/DENL17_uc.csv")
 
# two records with missing values in predictors are removed 
DENL_2017=DENL_2017%>%filter(DENL_2017!=-10000)
locations_sf = st_as_sf(DENL_2017, coords = c("Longitude","Latitude"), crs=4642)
mapview(locations_sf)

#write.csv(DENL_2017,"/Users/menglu/Documents/GitHub/DENL17_hr.csv")
library(raster)
library(dplyr)
library(ranger)
library(APMtools)
data(global_annual) 
mergedall = read.csv("https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/DENL17_uc.csv")
names(mergedall)
gr = raster("~/Downloads/s5p2019_handpoly.tif")
gr1 = raster("~/Downloads/s5p2019_handpoly_aug_dec.tif")
g = list.files("/Volumes/Meng_Mac/S5p/", full.names = T) 

 

for ( i in g)
{
  ras = raster(i)
  mergedall = mergeraster2file(mergedall,ras, c("Longitude", "Latitude"), paste0("S5p_", substr(i, 29, 30)))
}
mergedall = mergeraster2file(mergedall,gr, c("Longitude", "Latitude"), "S5p_b_Aug")
mergedall = mergeraster2file(mergedall,gr1, c("Longitude", "Latitude"), "S5p_a_Aug")


y_var = "mean_value"
prestring =  "S5p|road|nightlight|population|trop|indu|temp|wind|elev|radi|urbantype$"
varstring = paste(prestring,y_var,sep="|")
#
merged = mergedall%>%dplyr::select(matches(varstring))%>%na.omit()
names(merged)
set.seed(1)
r = ranger(mean_value~., merged, importance = "permutation")
r
sort(r$variable.importance, decreasing = T)
getwd()

 
write.csv(merged, "~/Documents/GitHub/uncertainty/data_vis_exp/DENL17_S5p_uc.csv")
 cor(merged)
 plot(unlist(apply(merged,2,median)), typ = "l")




