# First download the data in your working directory
download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip", "countries.zip")
# Then unzip
unzip("countries.zip")
# Load maptools
library(maptools)
# Read in the shapefile
world <- readShapeSpatial("ne_10m_admin_0_countries.shp")
# Plot France
plot(world[world$ADMIN=="Germany",1]) 
plot(world[world$ADMIN=="Netherlands",1], add = T) 

# Column 1 is the id column, column "ADMIN" contains the country names
 