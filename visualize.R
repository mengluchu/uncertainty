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

DENL_2017= read.csv("/Users/menglu/Documents/GitHub/uncertainty/DENL17_uc.csv")

locations_sf = st_as_sf(DENL_2017, coords = c("Longitude","Latitude"), crs=4642)
mapview(locations_sf)