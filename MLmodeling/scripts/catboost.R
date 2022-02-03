install.packages('devtools')
devtools::install_github('catboost/catboost', subdir = 'catboost/R-package')
library(catboost)

catboost_LUR = function(variabledf, numtrees = 300, eta = 0.007, y_varname= "mean_value", training, test, grepstring, ...) {
   
  pre_mat = subset_grep(variabledf[training, ], grepstring)
  pre_mat = pre_mat%>%dplyr::select(-y_varname)
  label = variabledf[training, y_varname] 
  pool = catboost.load_pool(data = pre_mat, label = label)
    
  y_test = variabledf[test, y_varname]
  x_test = variabledf[test, ]
  x_test = catboost.load_pool(data = x_test )
  
  model = catboost.train(pool,  NULL,
                 params = list(loss_function = 'RMSE', eval_metric="RMSE",
                               iterations = numtrees, learning_rate = eta, metric_period=300))
  pred = catboost.predict(model, x_test)
   
  return(error_matrix(y_test, pred))
}
library(dplyr)
library(APMtools)
df = merged

catboost_crossvali = function(df, n){ 
smp_size <- floor(0.8 * nrow(df)) 
set.seed(n)
training<- sample(seq_len(nrow(df)), size = smp_size)
test = seq_len(nrow(df))[-training]
"road|nightlight|population|temp|wind|trop|indu|elev|radi"
catboost_LUR(df, numtrees =  3000 ,y_varname= y_var, training=training, test=test, grepstring =varstring)
} 
nboot = 20
V2 = lapply(1:nboot, df = merged,catboost_crossvali)
V2 = data.frame(XGB = rowMeans(data.frame(V2)))
V2 