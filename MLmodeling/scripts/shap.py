#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 13:47:48 2021

@author: menglu
"""
import shap
#import h2o
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import matplotlib.pyplot as plt
import xgboost as xgb
from sklearn.model_selection import train_test_split
import random
from importlib.metadata import version
from h2o.automl import H2OAutoML
import h2o
from h2o.estimators.random_forest import H2ORandomForestEstimator
from h2o.estimators import H2OXGBoostEstimator
from sklearn.metrics import r2_score, explained_variance_score, mean_absolute_error, mean_absolute_percentage_error
from sklearn.ensemble import RandomForestRegressor


def performance (actural, pred):
    r2 = r2_score(actural, pred)
    expvar = explained_variance_score(actural, pred)
    mae = mean_absolute_error(actural, pred)
    #RMSE
    rmse = np.sqrt(mean_squared_error(actural, pred))
    # for between data-set comparison
    mape = mean_absolute_percentage_error(actural, pred)
    m ={"R2":r2, "expvar":expvar,"MAE": mae, "RMSE": rmse, "MAPE": mape}
    return (m)

res = 100
random.seed(1)
  
''' 
# AMS dataset
ap = pd.read_csv("~/Documents/GitHub/AP_AMS/J1.csv")
X_regex='pop|nightlight|trop|ele|wind|temp|ind|GH|road|sdb|ssb|fea|ltc|mtb|building'
predvar = "Vierweekse_11"

#X=ap.drop(columns=list(ap.filter(regex=X_regex))) # drop way  
''' 

# Uncertainty dataset
ap = pd.read_csv("~/Documents/GitHub/uncertainty/data_vis_exp/DENL17_traf.csv")
X_regex="pop|nightlight|trop|ele|wind|temp|ind|GH|road|CorrectedTraffic"

X_regex2="pop|nightlight|trop|ele|wind|temp|ind|GH|road"
predvar = "wkd_day_value"
Traffic = pd.read_excel("~/Documents/work/Students/Foeke/CorrectedTraffic.xlsx")
ap2 = ap.join(Traffic[:-1])


 
#ap = ap.sample(frac=1)
if res == 100:
    ap = ap.drop(ap.filter(regex='_25$|_50$').columns, axis = 1)
ap2 = ap2.sample(frac=1) 
X = ap2.filter(regex=X_regex2)
X.columns
 
#X["tra_100"]=X.Daily_Traffic_5km*X.road_class_2_100
y= ap2[predvar] # nothing with 10， 11， 12 13 ,1, 5, 7，
 

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.01, random_state=7)

# one way of running the xgboost
d_train = xgb.DMatrix(X_train, label=y_train)
d_test = xgb.DMatrix(X_test, label=y_test)
#Train the model
params = {
    "eta": 0.002, 
    "subsample": 0.7,
    'max_depth': 6, 
    'alpha':0.1
}

model = xgb.train(params, d_train, 2000, evals = [(d_test, "test")], verbose_eval=100, early_stopping_rounds=20)
fig, axes = plt.subplots(nrows = 1,ncols = 1,constrained_layout=True) 
xgb.plot_importance(model, grid=False, importance_type='gain', title='Feature importance', max_num_features=20)


#Differneces between SHAP and XGB variable importance.
shap_values = shap.TreeExplainer(model).shap_values(X) # mean importance

fig, axes = plt.subplots(nrows = 1, ncols = 1, constrained_layout=True) 
shap.summary_plot(shap_values, X.columns, plot_type="bar")

#fig.savefig('/Volumes/Meng_Mac/obia/temp/tree.png')  

# Negative and positive contributions, and local interpretability: the effect of each individual sample.  
fig, axes1 = plt.subplots(nrows = 1,ncols = 1,constrained_layout=True) 
shap.summary_plot(shap_values, X,   cmap="BrBG",max_display=10 ) 
 
fig.savefig('/Volumes/Meng_Mac/xgbshap.png')  



rf = RandomForestRegressor(n_estimators = 1000, random_state = 42)
# Train the model on training data
rf.fit(X_train, y_train)
rfshap_values = shap.TreeExplainer(rf).shap_values(X) # mean importance
shap.summary_plot(rfshap_values, X,   cmap="BrBG",max_display=10 ) 
 
fig.savefig('/Volumes/Meng_Mac/xgbshap.png')  



#fig.savefig("shap.png")
#road_class_2_100 affects a few predictions by a large amount, while population_1000 affects all predictions by a smaller amount.
shap.dependence_plot("BufferfoSpaJoin_5.CorrectedTraffic2km", shap_values, X, interaction_index = "road_class_2_100" )
shap.dependence_plot("BufferfoSpaJoin_5.CorrectedTraffic1km", shap_values, X)
plt.savefig('/Volumes/Meng_Mac/dep1km.png')
 
#shap.dependence_plot("tra_100", shap_values, X) 
plt.close("all")

''' AMS
shap.dependence_plot("ltcBuA", shap_values, X) # For 9， negative effect
shap.dependence_plot("features3", shap_values, X) # For 8, building density and volumnn density are consitent, interact with population
shap.dependence_plot("building_density", shap_values, X) # For 4, weird interaction with road_2_5000
shap.dependence_plot("sdbAre", shap_values, X) # for 6, interact with population
shap.dependence_plot("ssbElo", shap_values, X) #For 3

shap.dependence_plot("mtbAli", shap_values, X)
'''

#h2o

h2o.connect()
#if you cannot connect to h2o
#h2o.init(ip="127.0.0.1")
train, test = train_test_split(ap,test_size=0.1,random_state=1234)
#df_hex = h2o.H2OFrame(ap)
train_hex = h2o.H2OFrame(train)
test_hex = h2o.H2OFrame(test)

predictors =list(ap.filter(regex='pop|nightlight|trop|ele|wind|temp|ind|GH|road|sdb|ssb|fea|ltc'))
target="Vierweekse_2"

#prediction and add to the actual
def actual_predict(model,test_hex,target):
    y_pred = model.predict(test_hex).as_data_frame()
    y_actual = test_hex[target].as_data_frame()
    df_actual_predict = pd.concat([y_actual,y_pred],axis=1)
    df_actual_predict.columns = ['actual','pred']
    return(df_actual_predict)


RF_modl = H2ORandomForestEstimator(
        model_id = 'RF',
        ntrees = 300,
        nfolds=10,
        min_rows=100,
        seed=1234)

RF_modl.train(predictors,target,training_frame=train_hex)



XGB_modl = H2OXGBoostEstimator(
        model_id = 'XGB_modl',
        fold_assignment = "Modulo", #fixed fold, for comparing with other models, you can also try "Auto", "Modulo", and "Stratified"
        max_depth = 5, 
         eta = 0.002, 
        subsample=0.8,
        keep_cross_validation_predictions=True,
        ntrees = 1000,
        nfolds=10,
        min_rows=100,
        seed=1234)
XGB_modl.train(predictors,target,training_frame=train_hex)
#a = XGB_modl.cross_validation_predictions()
#len(a)

XGB_actual_predict = actual_predict(XGB_modl,test_hex,target)


# performance on the validation data for each fold. 
performance(XGB_actual_predict['actual'], XGB_actual_predict['pred'])

# we can compare with model performance on training 
XGB_modl.model_performance()
 

XGB_roc_auc_value

exm = XGB_modl.explain(test_hex)
print(exm.keys())
#or
#exm = XGB_modl.explain_row(test_hex, row_index = 0)
exm["varimp"] 

#['residual_analysis', 'varimp', 'shap_summary', 'pdp', 'ice']

RF_actual_predict = actual_predict(RF_modl,test_hex,target)
RF_actual_predict.head()

aml = H2OAutoML(max_models=20, seed=1)
aml.train(x=predictors, y=target, training_frame=train_hex)

# View the AutoML Leaderboard
lb = aml.leaderboard
lb.head(rows=lb.nrows)  




# run XGB using  the scikit-learn API
#params = {'max_depth': 6, 'eta': 0.002, 'num_class': 1, 'n_estimators':500}
#model = xgb.XGBModel(**params)
#model.fit(X_train, y_train)


#model.fit(X, y)

# fig, axes = plt.subplots(nrows = 2,ncols = 3,constrained_layout=True) 

# axes = axes.ravel()
# var = ["Daily_Traffic_500m", "Daily_Traffic_2km", "road_class_1_1000", "population_1000", "road_class_2_100","road_class_2_100"]
# for i in range(len(var)):
#     shap.plots.partial_dependence(
#         var[i], model.predict, X100,ice=False,
#         model_expected_value=True, feature_expected_value=True,ax=ax1,show=False, 
#     )
#     print(i)

np.corrcoef(ap['Daily_Traffic_500m'],ap['wkd_day_value'])
 
shao.plots.scatter(ap['Daily_Traffic_100m'], ap['wkd_day_value'], "o")
plt.plot(ap['Daily_Traffic_2km'], ap['wkd_day_value'], "o")
 
plt.show()

'Daily_Traffic_100m',
       'Daily_Traffic_500m', 'Daily_Traffic_1km', 'Daily_Traffic_2km',
 

   
       
# explain the GAM model with SHAP
explainer_ebm = shap.Explainer(clf.predict, )
shap_values_ebm = explainer_ebm(X)

# make a standard partial dependence plot with a single SHAP value overlaid
fig,ax = shap.partial_dependence_plot(
    "RM", model_ebm.predict, X, model_expected_value=True,
    feature_expected_value=True, show=False, ice=False,
    shap_values=shap_values_ebm[sample_ind:sample_ind+1,:]
)

# explain the GAM model with SHAP
explainer_xgb = shap.Explainer(clf, X)
shap_values_xgb = explainer_xgb(X)

# make a standard partial dependence plot with a single SHAP value overlaid
fig,ax = shap.partial_dependence_plot(
    "RM", model_xgb.predict, X, model_expected_value=True,
    feature_expected_value=True, show=False, ice=False,
    shap_values=shap_values_ebm[sample_ind:sample_ind+1,:]
)