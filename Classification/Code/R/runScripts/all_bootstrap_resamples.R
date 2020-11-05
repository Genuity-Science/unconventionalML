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
args = commandArgs(trailingOnly=TRUE)
ii = as.numeric(args[1])
st_idx = as.numeric(args[2])
stop_idx = as.numeric(args[3])
#ii = 5
#st_idx = 1
#stop_idx = 20
print(args)
# Load Data -- LUM A vs. LUM B
# load('/Users/ncilfone/Documents/Data_Files/TCGA/lumAB/data5_lumAB_all.RData')
files=c("brcaMatchedTN","ERpn","kirckirp","luadlusc" ,"lumAB")
base_dir = '/boston_ailab/users/rli/quantum/Data/'
save_dir = '/boston_ailab/users/rli/quantum/Results/'
positive_classes = c("tumor","Positive","kirc","luad","Luminal_A")
classes_levels = list(c("normal","tumor"),c("Negative","Positive"),c("kirp","kirc"),c("lusc","luad"),c("Luminal_B","Luminal_A"))

methods = c('glmnet','glmnet','svmLinear2','rf','nb') #('xgbLinear','svmRadial','svmLinear')
method_disp = c('lasso','ridge','svm','rf','nb')
lassoGrid=expand.grid(.alpha=1, .lambda=seq(0, 100, by = 0.1))
ridgeGrid=expand.grid(.alpha=0, .lambda=seq(0, 100, by = 0.1))
rfGrid=expand.grid(.mtry=c(5:30))#, .ntree=c(1000, 1500, 2000))
tuneGrids=list(lassoGrid,ridgeGrid,NULL,rfGrid,NULL)
cvCtrl = trainControl(method = "cv", number = 10,classProbs = TRUE, allowParallel = TRUE)
n_splits = 100
n_pc = 44
#for (ii in 1:1) { #length(files)) {
  cat("Dataset: ",files[[ii]],"\n")
  
  ### Parameters
  # Num PCs
  pos_class = positive_classes[[ii]]
  
  n=1;
  # Stats Mat
  info=data.frame(dataset=character(),method=character(),tr_acc=double(),tst_acc=double(),
                  tr_bacc=double(),tst_bacc=double(),tr_auroc=double(), tst_auroc=double(),
                  tr_prec=double(), tst_prec=double(), tr_recall=double(), tst_recall=double(),
                  tr_F1=double(),tst_F1=double(),stringsAsFactors=FALSE)
  
  cm_train_list = list()
  cm_test_list = list()
  response_train_list = list()
  response_test_list = list()
  
  mat = readMat(paste(base_dir,files[[ii]],'_bootstrap_resamples.mat',sep=''))
  traindatas = unlist(mat$traindatas,recursive=F)
  testdatas = unlist(mat$testdatas,recursive=F)

  # Run through # of Splits
  for (j in seq(st_idx,n_splits,stop_idx)) {
    cat("Split: ",j,"\n") 
#    mat = readMat(paste(files[[ii]], "_bootstrap_resamples/resample_", as.character(j), "_data.mat",sep=""))

#    pc_train = mat$traindata[,-1]
#    pc_test = mat$testdata[,-1]
#    class_train = as.factor(mat$traindata[,1])
#    class_test = as.factor(mat$testdata[,1])
    pc_train = traindatas[[j]][,2:(n_pc+1)]
    pc_test = testdatas[[j]][,2:(n_pc+1)]
    class_train = as.factor(traindatas[[j]][,1])
    class_test = as.factor(testdatas[[j]][,1])
    levels(class_train) = classes_levels[[ii]]
    levels(class_test) = classes_levels[[ii]]

    z_pc_train = pc_train
    z_pc_test = pc_test
#    cat("\n---------------\nNormalize PC Data...\n---------------\n")
#    pc_train_mean = apply(pc_train, 2, mean)
#    pc_train_sd = apply(pc_train, 2, sd)
#    # Train
#    z_pc_train = sweep(pc_train, 2, pc_train_mean, "-")
#    z_pc_train = sweep(z_pc_train, 2, pc_train_sd, "/")
#    # Test
#    z_pc_test = sweep(pc_test, 2, pc_train_mean, "-")
#    z_pc_test = sweep(z_pc_test, 2, pc_train_sd, "/")
#    # Valid
    colnames(z_pc_train) = paste('PC_', seq(1,ncol(z_pc_train)), sep = '')
    colnames(z_pc_test) = paste('PC_', seq(1,ncol(z_pc_test)), sep = '')
#    
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
#    class_train = as.factor(mat$traindata[,1])
#    class_test = as.factor(mat$testdata[,1])
      cm_test = confusionMatrix(pred_test, class_test, positive = pos_class)
      
      # calculate response probabilities
      response_train = predict(fit,z_pc_train,type="prob")
      response_test = predict(fit,z_pc_test,type="prob")
      response_train_list[[n]] = response_train
      response_test_list[[n]] = response_test
      
      # AUC train and test
      auc.train = ModelMetrics::auc(as.numeric(class_train), predict(fit, newdata = z_pc_train, type = 'prob')[,1])
      auc.test = ModelMetrics::auc(as.numeric(class_test), predict(fit, newdata = z_pc_test, type = 'prob')[,1])
      
      # # ROC train and test
      # roc.train = auc(multcap(response = class_train, predicted = data.matrix(response_train)))
      # roc.test = auc(multcap(response = class_test, predicted = data.matrix(response_test)))
      
      info[n,'method']=method_disp[[m]]
      info[n,'tr_acc']=cm_train$overall["Accuracy"]
      info[n,'tst_acc']=cm_test$overall["Accuracy"]
      cm_train_list[[n]] = cm_train
      cm_test_list[[n]] = cm_test
      
      
      # need code that gets same data, vector or matrix, just a comma difference in notation!
      info[n,'tr_bacc']=mean(cm_train$byClass['Balanced Accuracy'])
      info[n,'tst_bacc']=mean(cm_test$byClass['Balanced Accuracy'])
      info[n,'tr_prec']=mean(cm_train$byClass['Precision'])
      info[n,'tst_prec']=mean(cm_test$byClass['Precision'])
      info[n,'tr_recall']=mean(cm_train$byClass['Recall'])
      info[n,'tst_recall']=mean(cm_test$byClass['Recall'])
      info[n,'tr_F1']=mean(cm_train$byClass['F1'])
      info[n,'tst_F1']=mean(cm_test$byClass['F1'])
      
      info[n,'tr_auroc']=auc.train
      info[n,'tst_auroc']=auc.test
      n = n+1 
    }
  }
  save_vars = list(cm_train_list,cm_test_list,info,response_train_list,response_test_list)
  names(save_vars) = c("cm_train_list","cm_test_list","info","response_train_list","response_test_list")
  saveRDS(save_vars,paste(save_dir,files[[ii]],"_max_pcs_multirun_save",st_idx, ".RDS",sep=""))
