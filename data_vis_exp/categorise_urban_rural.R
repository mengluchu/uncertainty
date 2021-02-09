mergedall = read.csv("https://raw.githubusercontent.com/mengluchu/uncertainty/master/data_vis_exp/DENL17_uc.csv")
mergedall = mergedall%>%dplyr::select(-c(X, X.1, X.2, X.3))
sep = c( -0.001, 300, 1500,    max(mergedall$population_1000))
mergedall$urbantype = as.factor(cut(mergedall$population_1000, sep, label = c(1:(length(sep)-1))))
mergedall$urbantype_chara = as.factor(cut(mergedall$population_1000, sep, label = c("rural", "suburban", "urban")))
j = 0
for(i in 1:nrow(mergedall)){
  if(mergedall$AirQualityStationArea[i] == mergedall$urbantype_chara[i])
    j = j+1
}
j
write.csv(mergedall, "~/Documents/GitHub/uncertainty/data_vis_exp/DENL17_uc.csv")
head(mergedall)