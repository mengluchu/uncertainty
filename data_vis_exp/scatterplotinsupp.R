 
library(dplyr)
library(tidyr)
library(ggplot2)
 
 
scatterplot = function(variabledf, y_name, fittingmethod = "lm", ncol = 6) {
  # quote(y_var)
  variabledf %>% gather(VAR, predictors, -y_name) %>% 
    ggplot(aes_string(x = "predictors", y = y_name)) + geom_point() + 
    facet_wrap(~VAR, scales = "free", ncol=ncol) + xlab ("predictor variables") +
    ylab ( "mean NO2" ) + theme(
     strip.text = element_text( size=14),
      axis.title.x = element_text( size=20, face="bold"))+theme(
      axis.title.y = element_text( size=20, face="bold"))+
    stat_smooth(method = fittingmethod)  
}

 
 
resolution = 100 
# variables
y_var = "mean_value"
prestring =  "road|nightlight|population|temp|wind|trop|indu|elev|radi"
varstring = paste(prestring,y_var,sep="|")
 
mergedall = read.csv("https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/DENL17_uc.csv")

if (resolution ==100)
{
  mergedall = mergedall%>%dplyr::select(-c(industry_25,industry_50,road_class_1_25,road_class_1_50,road_class_2_25,road_class_2_50,   road_class_3_25,road_class_3_50))
}  
# data to use
merged = mergedall%>%dplyr::select(matches(varstring))%>% na.omit() # there is actually no na in this file, but for now RF and LA doesnt deal with missing data, leave out for quick examination 


merged%>%scatterplot(y_name = "mean_value", fittingmethod = "gam")
ggsave("~/Documents/GitHub/uncertainty/data_vis_exp/scat.png",width = 16, height = 16)
 