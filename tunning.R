library(caret)
library(xgboost)
DENL_2017= read.csv("/Users/menglu/Documents/GitHub/uncertainty/DENL17_uc.csv")
y_var = "mean_value"
prestring =  "road|nightlight|population|temp|wind|trop|indu|elev|radi"
varstring = paste(prestring,y_var,sep="|")
DENL_2017 = DENL_2017%>%dplyr::select(matches(varstring))

xgboostgrid = expand.grid(nrounds = seq(200, 3000, by = 200), max_depth = 3:6, eta = seq(0.001, 0.1, by = 0.05),
                          colsample_bytree = 1,
                          gamma = 3, 
                          
                          min_child_weight = 1,
                          subsample = c(0.7,0.9)) 
#gamma: Minimum loss reduction required to make a further partition on a leaf node of the tree. The larger gamma is, the more conservative the algorithm will be.

trainControl <- trainControl(method="cv", number=5, savePredictions = "final", allowParallel = T) #5 - folds
# train the model
model <- train(mean_value~., data=DENL_2017, method="xgbTree", trControl=trainControl, tuneGrid =xgboostgrid)

model

ggplot(model)
ggsave("xgboosttunning.png")
 