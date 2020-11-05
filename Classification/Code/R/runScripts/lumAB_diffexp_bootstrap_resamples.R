rm(list = ls())
gc()
library(feather)
library(caret)
library(stringr)
library(HandTill2001)
library(MLmetrics)

options(stringsAsFactors=F)
Sys.setenv(TZ="US/Eastern")
args = commandArgs(trailingOnly=TRUE)
st_idx = as.numeric(args[1])
stop_idx = as.numeric(args[2])

base_dir = '/boston_ailab/users/rli/quantum/Data/updated/'
save_dir = '/boston_ailab/users/rli/quantum/Results/'
pos_class = 'Luminal_A'
classes_levels = c("Luminal_B","Luminal_A")
f = 'lumAB'

methods = c('glmnet','glmnet','svmLinear2','rf','nb') #('xgbLinear','svmRadial','svmLinear')
method_disp = c('lasso','ridge','svm','rf','nb')
lassoGrid=expand.grid(.alpha=1, .lambda=seq(0, 100, by = 0.1))
ridgeGrid=expand.grid(.alpha=0, .lambda=seq(0, 100, by = 0.1))
rfGrid=expand.grid(.mtry=c(5:30))#, .ntree=c(1000, 1500, 2000))
tuneGrids=list(lassoGrid,ridgeGrid,NULL,rfGrid,NULL)
cvCtrl = trainControl(method = "cv", number = 10,classProbs = TRUE, allowParallel = TRUE)
n_splits = 100
n_pc = 44

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

cat("\n---------------\nLoading train data...\n---------------\n")
train_data = read_feather(paste(base_dir,"data5_",f,"_train.feather",sep=""))

cat("\n---------------\nLoading test data...\n---------------\n")
test_data = read_feather(paste(base_dir,"data5_",f,"_test.feather",sep=""))

block_path = '/boston_ailab/users/rli/quantum/Data/'

ids = read.table(paste(block_path,'lumAB_top44_diffexp_genes.txt',sep=''))
all_data = rbind(train_data,test_data)
# fix factor names

x_data = all_data[ids$V1]
y_data = all_data$cancer
y_data = gsub(" ", "_", y_data, fixed = TRUE)
y_data = as.factor(y_data)
splits = as.matrix(read.table('/boston_ailab/users/rli/quantum/Data/lumAB_bootstrap_resamples.txt'))

# Run through # of Splits
for (j in seq(st_idx,n_splits,stop_idx)) {
  cat("Split: ",j,"\n") 

  class_train = y_data[splits[j,]]
  class_test = y_data[-splits[j,]]
  # Split Data
  data_train = x_data[splits[j,],]
  data_test = x_data[-splits[j,],]

  train_mean = apply(data_train, 2, mean)
  train_sd = apply(data_train, 2, sd)
  # Train
  x_train = sweep(data_train, 2, train_mean, "-")
  x_train = sweep(x_train, 2, train_sd, "/")
  # Test
  x_test = sweep(data_test, 2, train_mean, "-")
  x_test = sweep(x_test, 2, train_sd, "/")

  for (m in 1:length(methods)) {
    method=methods[m]
    print(method_disp[[m]])
    
    ######### Model Runs
    fit = train(x_train,class_train,method=method,trControl = cvCtrl,tuneGrid = tuneGrids[[m]])
    
    # Class returns the class of the prediction 
    pred_train = predict(fit,x_train)
    pred_test = predict(fit,x_test)
    
    # Class Stats
    cm_train = confusionMatrix(pred_train, class_train, positive = pos_class)
    cm_test = confusionMatrix(pred_test, class_test, positive = pos_class)
    
    # calculate response probabilities
    response_train = predict(fit,x_train,type="prob")
    response_test = predict(fit,x_test,type="prob")
    response_train_list[[n]] = response_train
    response_test_list[[n]] = response_test
    
    # AUC train and test
    auc.train = ModelMetrics::auc(as.numeric(class_train), predict(fit, newdata = x_train, type = 'prob')[,1])
    auc.test = ModelMetrics::auc(as.numeric(class_test), predict(fit, newdata = x_test, type = 'prob')[,1])
    
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
  saveRDS(save_vars,paste(save_dir,"lumAB_diffexp_gene_max_pcs_multirun_save",st_idx, ".RDS",sep=""))
