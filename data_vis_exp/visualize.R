library(dplyr)
library("mapview")
 
library("sf")
# Linking to GEOS 3.8.0, GDAL 3.1.0, PROJ 7.0.0
library("stars")
library("raster")
library("RColorBrewer")
install_github("mengluchu/APMtools") 


library(APMtools)
mapviewOptions(
  basemaps = c("OpenStreetMap.Mapnik","Esri.OceanBasemap")
  , raster.palette = colorRampPalette(rev(brewer.pal(9, "YlGnBu")))
  , vector.palette = colorRampPalette(brewer.pal(9, "PuBuGn"))
  , na.color = "gray"
  , layers.control.pos = "topright"
  , viewer.suppress = TRUE # open browser
)

DENL_2017= read.csv("/Users/menglu/Documents/GitHub/uncertainty/DENL17_uc.csv")

locations_sf = st_as_sf(DENL_2017, coords = c("Longitude","Latitude"), crs=4642)
mapview(locations_sf)
r = raster("/Users/menglu/Documents/GitHub/uncertainty/mean_NLDE.tif")
crop_center(r, 150,  vis= T)
