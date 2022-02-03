
#==========================================
# test distributin, figure 1 manuscript
#============================================
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(goft) # for gamma test

mergedall = read.csv("https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/DENL17_uc.csv")
file_url <- "https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/missingstation.rda?raw=true"
load(url(file_url)) # remove stations contain more less than 25% of data

d =mergedall%>%filter(!(AirQualityStation %in% msname)) #474
d$y = d$mean_value # 
shapiro.test(d$y)
gamma_test(d$y) #goft
hist = ggplot(d, aes(x=y))+xlab( expression(NO[2]~(~mu~g/m^3))) + 
  geom_histogram(binwidth = 1, aes(y=..density..), colour="black", fill="white")+geom_density(alpha=.1, fill="#FF6666") 

original <- ggplot(d, aes(Longitude, Latitude)) +
  geom_point(aes(colour=AirQualityStationType)) +
  #  scale_colour_gradient(name = expression(Level~of~NO[2]), low = "yellow", high = "red") + 
  theme(plot.title = element_text(hjust = 0))+ 
  scale_color_manual(values =brewer.pal(3,"Set2"), name = expression(NO[2]~(~mu~g/m^3)))+
  geom_polygon(data = shapefile, aes(x = long, y = lat, group = group), colour = "black", fill = NA)
original
#ggsave("~/Documents/GitHub/uncertainty/histogram_NO2.png")
qq=ggqqplot(d$y)
#ggsave("~/Documents/GitHub/uncertainty/ggplot_NO2.png")
ggplot(d, aes(y))+geom_histogram()

g=grid.arrange(original, hist,qq, ncol = 3)
ggsave("~/Documents/GitHub/uncertainty/histqq_NO2.png",g, height = 6, width = 18)
#p-value < 0.05 implying that the distribution of the data are significantly different from normal distribution. 
