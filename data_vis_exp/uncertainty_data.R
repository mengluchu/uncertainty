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

EU = merge(merged, EUmeta, by=c("Longitude", "Latitude"), all.x =T)
EU = EU%>%filter(!is.na(Countrycode))

#EU 
locations_sf = st_as_sf(EU, coords = c("Longitude","Latitude"), crs=4642)
 # https://leaflet-extras.github.io/leaflet-providers/preview/
mapview(locations_sf)

#NL,DE
DENL_2017 = EU%>%filter(Countrycode=="DE"|Countrycode=="NL")

# two records with missing values in predictors are removed 
DENL_2017%>%filter(DENL_2017!=-10000)
locations_sf = st_as_sf(DENL_2017, coords = c("Longitude","Latitude"), crs=4642)
mapview(locations_sf)

write.csv(DENL_2017,"/Users/menglu/Documents/GitHub/DENL17_uc.csv")
