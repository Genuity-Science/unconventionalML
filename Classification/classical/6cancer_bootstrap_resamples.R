# Run all classical methods for bootstrap resampling on 6cancer dataset. Bootstrap resampling 
# refers to taking all the data and resampling 100 cuts of 80% training and
# 20% testing.  
# Author: Richard Li. Adapted from code by Nicholas Cilfone and Omar Gamel

rm(list = ls())
gc()
library(feather)
library(caret)
library(stringr)
library(HandTill2001)
library(MLmetrics)
library(PRROC)
library(R.matlab)

Sys.setenv(TZ="US/Eastern")

# Load Data 
classes_levels = c("lihc","brca","lgg","coad",'kidn','lung')
methods = c('glmnet','glmnet','svmLinear2','rf','nb') #('xgbLinear','svmRadial','svmLinear')
method_disp = c('lasso','ridge','svm','rf','nb')
lassoGrid=expand.grid(.alpha=1, .lambda=seq(0, 100, by = 0.1))
ridgeGrid=expand.grid(.alpha=0, .lambda=seq(0, 100, by = 0.1))
rfGrid=expand.grid(.mtry=c(1:50))#, .ntree=c(1000, 1500, 2000))
tuneGrids=list(lassoGrid,ridgeGrid,NULL,rfGrid,NULL)
cvCtrl = trainControl(method = "cv", number = 10,classProbs = TRUE)

n_pc = 13
n_splits = 100

cat("Dataset: 6cancer \n")
mat = readMat("Data/6cancer_bootstrap_resamples.mat")
traindatas = unlist(mat$traindatas,recursive=FALSE)
testdatas = unlist(mat$testdatas,recursive=FALSE)

### Parameters
# Num PCs
pos_class = NULL

# Stats Mat
info=data.frame(dataset=character(),method=character(),tr_acc=double(),tst_acc=double(),
                tr_bacc=double(),tst_bacc=double(),tr_auroc=double(), tst_auroc=double(),
                tr_prec=double(), tst_prec=double(), tr_recall=double(), tst_recall=double(),
                tr_F1=double(),tst_F1=double(),stringsAsFactors=FALSE)

n=1;
cm_train_list = list()
cm_test_list = list()
response_train_list = list()
response_test_list = list()
# Run through # of Splits
for (j in 1:n_splits) {
  
  pc_train = traindatas[[j]][,2:(n_pc+1)]
  pc_test = testdatas[[j]][,2:(n_pc+1)]
  class_train = as.factor(traindatas[[j]][,1])
  class_test = as.factor(testdatas[[j]][,1])
  levels(class_train) = classes_levels
  levels(class_test) = classes_levels
  
  cat("\n---------------\nCut",as.character(j),"...\n---------------\n")
  colnames(z_pc_train) = paste('PC_', seq(1,ncol(z_pc_train)), sep = '')
  colnames(z_pc_test) = paste('PC_', seq(1,ncol(z_pc_test)), sep = '')
  
  for (m in 1:length(methods)) {
    method=methods[m]
    print(method_disp[[m]])
    
    ######### Model Runs
    fit = train(z_pc_train,class_train,method=method,trControl = cvCtrl,tuneGrid = tuneGrids[[m]])
    
    # Class returns the class of the prediction 
    pred_train = predict(fit,z_pc_train)
    pred_test = predict(fit,z_pc_test)
    
    # Class Stats
    cm_train = confusionMatrix(pred_train, class_train, positive = pos_class)
    cm_test = confusionMatrix(pred_test, class_test, positive = pos_class)
    
    # calculate response probabilities
    response_train = predict(fit,z_pc_train,type="prob")
    response_test = predict(fit,z_pc_test,type="prob")
    response_train_list[[n]] = response_train
    response_test_list[[n]] = response_test
    
    # ROC train and test
    roc.train = auc(multcap(response = class_train, predicted = data.matrix(response_train)))
    roc.test = auc(multcap(response = class_test, predicted = data.matrix(response_test)))
    info[n,'dataset']='6cancer' 
    info[n,'method']=method_disp[[m]]
    info[n,'tr_acc']=cm_train$overall["Accuracy"]
    info[n,'tst_acc']=cm_test$overall["Accuracy"]
    cm_train_list[[n]] = cm_train
    cm_test_list[[n]] = cm_test
    
    
    # need code that gets same data, vector or matrix, just a comma difference in notation!
    info[n,'tr_bacc']=mean(cm_train$byClass[,'Balanced Accuracy'])
    info[n,'tst_bacc']=mean(cm_test$byClass[,'Balanced Accuracy'])
    info[n,'tr_prec']=mean(cm_train$byClass[,'Precision'])
    info[n,'tst_prec']=mean(cm_test$byClass[,'Precision'])
    info[n,'tr_recall']=mean(cm_train$byClass[,'Recall'])
    info[n,'tst_recall']=mean(cm_test$byClass[,'Recall'])
    info[n,'tr_F1']=mean(cm_train$byClass[,'F1'])
    info[n,'tst_F1']=mean(cm_test$byClass[,'F1'])
    
    info[n,'tr_auroc']=roc.train
    info[n,'tst_auroc']=roc.test
    n = n+1

  }
  # save periodically, just in case
  if ( j %% 25 == 0) {
      save.image("tmp.RData")
  }
}
save_vars = list(cm_train_list,cm_test_list,info,response_train_list,response_test_list)
names(save_vars) = c("cm_train_list","cm_test_list","info","response_train_list","response_test_list")
saveRDS(save_vars,paste("6cancer_bootstrap_resamples_save.RDS",sep=""))
