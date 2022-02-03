# following Regional Working Paper 2014 A harmonised
#definition of
#cities and rural
#areas: the new degree of urbanisation
#Lewis Dijkstra and Hugo Poelman
#https://ec.europa.eu/regional_policy/sources/docgener/work/2014_01_new_urban.pdf

# urban correspond to suburban and city centre correspond to urban. 
mergedall = read.csv("https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/DENL17_uc.csv")
mergedall = mergedall%>%dplyr::select(-c(X, X.1, X.2, X.3))
sep = c( -0.001, 300*3.14, 1500*3.14,    max(mergedall$population_1000))
mergedall$urbantype = as.factor(cut(mergedall$population_1000, sep, label = c(1:(length(sep)-1))))
mergedall$urbantype_chara = as.factor(cut(mergedall$population_1000, sep, label = c("rural", "suburban", "urban")))
j = 0
for(i in 1:nrow(mergedall)){
  if(mergedall$AirQualityStationArea[i] == mergedall$urbantype_chara[i])
    j = j+1
}
j
table(mergedall$urbantype_chara)
write.csv(mergedall, "~/Documents/GitHub/uncertainty/data_vis_exp/DENL17_uc.csv")
mergedall$AirQualityStationArea