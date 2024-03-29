rm(list = ls())
gc()
library(feather)
library(caret)
library(stringr)
library(HandTill2001)
library(MLmetrics)
#library(PRROC)
library(R.matlab)

Sys.setenv(TZ="US/Eastern")

# Load Data -- LUM A vs. LUM B
files=c("ERpn","lumAB")
base_dir = '/boston_ailab/users/rli/quantum/Data/'
positive_classes = c("Positive","Luminal_A")
classes_levels = list(c("Negative","Positive"),c("Luminal_B","Luminal_A"))
all_fracs = list(c(0.06,seq(0.1,0.95,0.05)), c(0.18,seq(0.2,0.95,0.05)))

args = commandArgs(trailingOnly=T)
ii = as.numeric(args[1])
jj = as.numeric(args[2])
f = files[[ii]]
fracs = all_fracs[[ii]]
pos_class = positive_classes[[ii]]
class_levels = classes_levels[[ii]]

base_dir = paste('/boston_ailab/users/rli/quantum/Data/',f,'_splits/',sep="")

methods = c('glmnet','glmnet','svmLinear2','rf','nb') #('xgbLinear','svmRadial','svmLinear')
method_disp = c('lasso','ridge','svm','rf','nb')
lassoGrid=expand.grid(.alpha=1, .lambda=seq(0, 100, by = 0.1))
ridgeGrid=expand.grid(.alpha=0, .lambda=seq(0, 100, by = 0.1))
rfGrid=expand.grid(.mtry=c(1:50))#, .ntree=c(1000, 1500, 2000))
tuneGrids=list(lassoGrid,ridgeGrid,NULL,rfGrid,NULL)
cvCtrl = trainControl(method = "cv", number = 10,classProbs = TRUE)

# Stats Mat
n_pc = 44
n_splits = 50


#fracs = c(0.06,seq(0.1,0.95,0.05))
#fracs = c(0.18,seq(0.2,0.95,0.05))
info=data.frame(frac=double(),method=character(),tr_acc=double(),tst_acc=double(),
              tr_bacc=double(),tst_bacc=double(),tr_auroc=double(), tst_auroc=double(),
              tr_prec=double(), tst_prec=double(), tr_recall=double(), tst_recall=double(),
              tr_F1=double(),tst_F1=double(),stringsAsFactors=FALSE)
n=1;
cm_train_list = list()
cm_test_list = list()
cm_val_list = list()
cm_exptest_list = list()
response_train_list = list()
response_test_list = list()
response_val_list = list()
response_exptest_list = list()

  cat(f, " Frac: ",fracs[[jj]],"\n")

  mat = readMat(paste(base_dir,"frac_",fracs[[jj]],"_data_resamples.mat",sep=""))
  traindatas = unlist(mat$traindatas,recursive=FALSE)
  exptestdatas = unlist(mat$exptestdatas,recursive=FALSE)
  valdatas = unlist(mat$valdatas,recursive=FALSE)
  testdatas = unlist(mat$testdatas,recursive=FALSE)
  
  # Run through # of Splits
  for (j in 1:n_splits) {
    x_train = traindatas[[j]][,2:(n_pc+1)]
    x_test = testdatas[[j]][,2:(n_pc+1)]
    x_val = valdatas[[j]][,2:(n_pc+1)]
    x_exptest = exptestdatas[[j]][,2:(n_pc+1)]

    class_train = as.factor(traindatas[[j]][,1])
    class_exptest = as.factor(exptestdatas[[j]][,1])
    class_val = as.factor(valdatas[[j]][,1])
    class_test = as.factor(testdatas[[j]][,1])
    levels(class_train) = class_levels
    levels(class_test) = class_levels
    levels(class_val) = class_levels
    levels(class_exptest) = class_levels
    
    cat("\n---------------\nCut",as.character(j),"...\n---------------\n")
    
    x_train_mean = apply(x_train, 2, mean)
    x_train_sd = apply(x_train, 2, sd)
    # Train
    z_x_train = sweep(x_train, 2, x_train_mean, "-")
    z_x_train = sweep(z_x_train, 2, x_train_sd, "/")
    # Test
    z_x_test = sweep(x_test, 2, x_train_mean, "-")
    z_x_test = sweep(z_x_test, 2, x_train_sd, "/")
    # Valid
    z_x_val = sweep(x_val, 2, x_train_mean, "-")
    z_x_val = sweep(z_x_val, 2, x_train_sd, "/")
    # Valid
    z_x_exptest = sweep(x_exptest, 2, x_train_mean, "-")
    z_x_exptest = sweep(z_x_exptest, 2, x_train_sd, "/")
    
    colnames(z_x_train) = paste('PC_', seq(1,ncol(z_x_train)), sep = '')
    colnames(z_x_test) = paste('PC_', seq(1,ncol(z_x_test)), sep = '')
    colnames(z_x_val) = paste('PC_', seq(1,ncol(z_x_val)), sep = '')
    colnames(z_x_exptest) = paste('PC_', seq(1,ncol(z_x_exptest)), sep = '')
    for (m in 1:length(methods)) {
      method=methods[m]
      print(method_disp[[m]])
      
      ######### Model Runs
      fit = train(z_x_train,class_train,method=method,trControl = cvCtrl,tuneGrid = tuneGrids[[m]])
      
      # Class returns the class of the prediction 
      pred_train = predict(fit,z_x_train)
      pred_test = predict(fit,z_x_test)
      pred_val = predict(fit,z_x_val)
      pred_exptest = predict(fit,z_x_exptest)
      
      # Class Stats
      cm_train = confusionMatrix(pred_train, class_train, positive = pos_class)
      cm_test = confusionMatrix(pred_test, class_test, positive = pos_class)
      cm_val = confusionMatrix(pred_val, class_val, positive = pos_class)
      cm_exptest = confusionMatrix(pred_exptest, class_exptest, positive = pos_class)
      
      # calculate response probabilities
      response_train = data.matrix(predict(fit,z_x_train,type="prob"))
      response_test = data.matrix(predict(fit,z_x_test,type="prob"))
      response_val = data.matrix(predict(fit,z_x_val,type="prob"))
      response_exptest = data.matrix(predict(fit,z_x_exptest,type="prob"))
      response_train[is.nan(response_train)] = 0.5
      response_test[is.nan(response_test)] = 0.5
      response_val[is.nan(response_val)] = 0.5
      response_exptest[is.nan(response_exptest)] = 0.5
#      response_train

      response_train_list[[n]] = response_train
      response_test_list[[n]] = response_test
      response_val_list[[n]] = response_val
      response_exptest_list[[n]] = response_exptest
      
      # ROC train and test
      roc.train = auc(multcap(response = class_train, predicted = response_train))
      roc.test = auc(multcap(response = class_test, predicted = response_test))
      roc.val = auc(multcap(response = class_val, predicted = response_val))
      roc.exptest = auc(multcap(response = class_exptest, predicted = response_exptest))
     
      info[n,'frac'] = fracs[[jj]] 
      info[n,'method']=method_disp[[m]]
      info[n,'tr_acc']=cm_train$overall["Accuracy"]
      info[n,'tst_acc']=cm_test$overall["Accuracy"]
      info[n,'val_acc']=cm_val$overall["Accuracy"]
      info[n,'exptst_acc']=cm_exptest$overall["Accuracy"]
      cm_train_list[[n]] = cm_train
      cm_test_list[[n]] = cm_test
      cm_val_list[[n]] = cm_val
      cm_exptest_list[[n]] = cm_exptest
      
      # need code that gets same data, vector or matrix, just a comma difference in notation!
      info[n,'tr_bacc']=mean(cm_train$byClass['Balanced Accuracy'])
      info[n,'tst_bacc']=mean(cm_test$byClass['Balanced Accuracy'])
      info[n,'val_bacc']=mean(cm_val$byClass['Balanced Accuracy'])
      info[n,'exptst_bacc']=mean(cm_exptest$byClass['Balanced Accuracy'])
      info[n,'tr_prec']=mean(cm_train$byClass['Precision'])
      info[n,'tst_prec']=mean(cm_test$byClass['Precision'])
      info[n,'val_prec']=mean(cm_val$byClass['Precision'])
      info[n,'exptst_prec']=mean(cm_exptest$byClass['Precision'])
      info[n,'tr_recall']=mean(cm_train$byClass['Recall'])
      info[n,'tst_recall']=mean(cm_test$byClass['Recall'])
      info[n,'val_recall']=mean(cm_val$byClass['Recall'])
      info[n,'exptst_recall']=mean(cm_exptest$byClass['Recall'])
      info[n,'tr_F1']=mean(cm_train$byClass['F1'])
      info[n,'tst_F1']=mean(cm_test$byClass['F1'])
      info[n,'val_F1']=mean(cm_val$byClass['F1'])
      info[n,'exptst_F1']=mean(cm_exptest$byClass['F1'])
      info[n,'tr_auroc']=roc.train
      info[n,'tst_auroc']=roc.test
      info[n,'val_auroc']=roc.val
      info[n,'exptst_auroc']=roc.exptest
      n = n+1 
    }
  }
  save_vars = list(cm_train_list,cm_test_list,cm_val_list,cm_exptest_list,info,response_train_list,response_test_list,response_val_list,response_exptest_list)
  names(save_vars) = c("cm_train_list","cm_test_list","cm_val_list","cm_exptest_list","info","response_train_list","response_test_list","response_val_list","response_exptest_list")
  saveRDS(save_vars,paste(f,"_frac_", as.character(fracs[[jj]]), "_save.RDS",sep=""))
