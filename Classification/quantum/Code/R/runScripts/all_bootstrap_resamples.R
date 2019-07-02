# Run all classical methods for bootstrap resampling. Bootstrap resampling 
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

# Load Data -- LUM A vs. LUM B
files=c("brcaMatchedTN","ERpn","kirckirp","luadlusc") #"lumAB"
base_dir = '~/Dropbox-Work/Wuxi/Data/'
positive_classes = c("tumor","Positive","kirc","luad") # "Luminal_A"
classes_levels = list(c("normal","tumor"),c("Negative","Positive"),c("kirp","kirc"),c("lusc","luad")) # c("Luminal_A","Luminal_B") 

# Parameters for CV
methods = c('glmnet','glmnet','svmLinear2','rf','nb') #('xgbLinear','svmRadial','svmLinear')
method_disp = c('lasso','ridge','svm','rf','nb')
lassoGrid=expand.grid(.alpha=1, .lambda=seq(0, 100, by = 0.1))
ridgeGrid=expand.grid(.alpha=0, .lambda=seq(0, 100, by = 0.1))
rfGrid=expand.grid(.mtry=c(1:50))#, .ntree=c(1000, 1500, 2000))
tuneGrids=list(lassoGrid,ridgeGrid,NULL,rfGrid,NULL)
cvCtrl = trainControl(method = "cv", number = 10,classProbs = TRUE)

n_pc = 44
n_splits = 100

for (ii in 1:length(files)) {
  cat("Dataset: ",files[[ii]],"\n")
  mat = readMat(paste("Data/",files[[ii]],"_bootstrap_resamples.mat",sep=""))
  traindatas = unlist(mat$traindatas,recursive=FALSE)
  testdatas = unlist(mat$testdatas,recursive=FALSE)

  ### Parameters
  # Num PCs
  pos_class = positive_classes[[ii]]
  
  # Intialize stats dataframe
  info=data.frame(dataset=character(),method=character(),tr_acc=double(),tst_acc=double(),
                tr_bacc=double(),tst_bacc=double(),tr_auroc=double(), tst_auroc=double(),
                tr_prec=double(), tst_prec=double(), tr_recall=double(), tst_recall=double(),
                tr_F1=double(),tst_F1=double(),stringsAsFactors=FALSE)

  # initialize counter
  n=1;

  # initialize list of confusion matrices
  cm_train_list = list()
  cm_test_list = list()

  # initialize list of response probabilities
  response_train_list = list()
  response_test_list = list()
  # Run through # of Splits
  for (j in 1:n_splits) {
    
    pc_train = traindatas[[j]][,2:(n_pc+1)]
    pc_test = testdatas[[j]][,2:(n_pc+1)]
    class_train = as.factor(traindatas[[j]][,1])
    class_test = as.factor(testdatas[[j]][,1])
    levels(class_train) = classes_levels[[ii]]
    levels(class_test) = classes_levels[[ii]]
    
    cat("\n---------------\nCut",as.character(j),"...\n---------------\n")
    pc_train_mean = apply(pc_train, 2, mean)
    pc_train_sd = apply(pc_train, 2, sd)
    # Train
    z_pc_train = sweep(pc_train, 2, pc_train_mean, "-")
    z_pc_train = sweep(z_pc_train, 2, pc_train_sd, "/")
    # Test
    z_pc_test = sweep(pc_test, 2, pc_train_mean, "-")
    z_pc_test = sweep(z_pc_test, 2, pc_train_sd, "/")
    # Valid
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
      info[n,'dataset']=files[[ii]]      
      info[n,'method']=method_disp[[m]]
      info[n,'tr_acc']=cm_train$overall["Accuracy"]
      info[n,'tst_acc']=cm_test$overall["Accuracy"]
      cm_train_list[[n]] = cm_train
      cm_test_list[[n]] = cm_test
      
      info[n,'tr_bacc']=mean(cm_train$byClass['Balanced Accuracy'])
      info[n,'tst_bacc']=mean(cm_test$byClass['Balanced Accuracy'])
      info[n,'tr_prec']=mean(cm_train$byClass['Precision'])
      info[n,'tst_prec']=mean(cm_test$byClass['Precision'])
      info[n,'tr_recall']=mean(cm_train$byClass['Recall'])
      info[n,'tst_recall']=mean(cm_test$byClass['Recall'])
      info[n,'tr_F1']=mean(cm_train$byClass['F1'])
      info[n,'tst_F1']=mean(cm_test$byClass['F1'])
      
      info[n,'tr_auroc']=roc.train
      info[n,'tst_auroc']=roc.test
      n = n+1

    }
    if ( j %% 25 == 0) {
        save.image("tmp.RData")
    }
  }
  save_vars = list(cm_train_list,cm_test_list,info,response_train_list,response_test_list)
  names(save_vars) = c("cm_train_list","cm_test_list","info","response_train_list","response_test_list")
  saveRDS(save_vars,paste(files[[ii]],"_bootstrap_resamples_save.RDS",sep=""))
  # should consider outputting mean and standard deviation so that matches format used previously 
}
