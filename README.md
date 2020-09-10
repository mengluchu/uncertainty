# uncertainty
### Uncertainty assessment comparing machine learning and statistical models. 
#### Data
Columns in the DENL17_uc.csv:

NO2 concentrations at various temporal aggregations (averaged according to local times):
•	Wkd_day_value: weekday (Monday to Friday) daytime (08:00 am – 22:00 pm) 
•	Wkd_night_value: weekday (Monday to Friday) nighttime (00: 00 – 07:00 am, 23:00 pm)
•	Wnd_day_value: weekend daytime
•	Wnd_night_value: weekend nighttime
•	Mean_value: annual mean ***This variable is the response variable currently in use for the model comparison***

Predictors are provided at each location, including

•	Road_class1_buffer: total highways length within different buffer rings (25m, 50m, 100m, 300m, 500 m, 800 m, 1000m, 3000 m) from OpenStreetMaps (OSM). You will see the radius size in meters as part of variable names, this applies to all the buffered variables.
•	Road_class2_buffer: primary roads
•	Road_class3_buffer: local roads
•	Industry_buffer: industrial area in different buffer rings. 
•	Nightlight_buffer: Earth nightlight from VIIRS satellite instrument, in different buffer rings.
•	Wind speed: wind speed at 2m height from ERA5-Land climate reanalysis, monthly average
•	Temperature: temperature at 2m height from ERA5-Land climate reanalysis, monthly average.
•	Population: population in different buffer rings from human settlement layers. GHS-POP R2019A
•	DEM: from MERIT
•	Radiation: global horizontal sun radiation

