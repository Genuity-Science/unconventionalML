# Output RDS file with annealing-type metrics for bootstrap resamples. Run after 
# doing analysis code in Matlab. Used to be consistent in how the performance
# metrics are calculated for the classical algorithms. Saves as .RDS file
# Author: Richard Li

rm(list=ls())
gc()

library(caret)
library(HandTill2001)
library(R.matlab)

n_splits = 100
base_dir = '~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/6cancer_bootstrap_resamples/'

# change the string to change method. Assumes that all saved in the same format; e.g.,
# ${method_str}_nsols20_pred_for_R.mat, where method_str can be "dw","sa","rand". 
# Any user-defined string is fine, as long as as this fileformat
method_str = 'dw' # 'dw','sa', 'field', or 'rand'.
sfx =  '_cinit8_nsols20_ntotsols1000'
# load data
mat = readMat(paste(base_dir,method_str,sfx,"_pred_for_R.mat",sep=''))

# define positive classes and levels
pos_class=NULL
classes = c("lihc","brca","lgg","coad",'kidn','lung')
info=data.frame(dataset=character(),method=character(),tr_acc=double(),tst_acc=double(),
                tr_bacc=double(),tst_bacc=double(),tr_auroc=double(), tst_auroc=double(),
                tr_prec=double(), tst_prec=double(), tr_recall=double(), tst_recall=double(),
                tr_F1=double(),tst_F1=double(),stringsAsFactors=FALSE)

response_trains = unlist(mat$y.pred.trains,recursive=FALSE)
response_tests = unlist(mat$y.pred.tests,recursive=FALSE)
for (j in 1:n_splits) {
 
  # response_train and response_test are the predicted probabilities for each 
  # class. Of size Nx6, where 6 is the number of classes
  response_train = response_trains[[j]]
  response_test = response_tests[[j]]

  # class test and class train are the actual classes
  class_test = factor(mat$y.tests[,j])
  class_train = factor(mat$y.trains[,j])
  
  pred_train = factor(max.col(response_train))
  pred_test = factor(max.col(response_test))
  if (nlevels(pred_train) != 6) {
    print(paste("Pred train iteration ", j, " does not have enough levels"))
  }
  if (nlevels(pred_test) != 6) {
    print(paste("Pred test iteration ", j, " does not have enough levels"))
  }
  colnames(response_train) = classes
  colnames(response_test) = classes
  levels(pred_train) = classes
  levels(pred_test) = classes
  levels(class_train) = classes
  levels(class_test) = classes
  
  # confusion matrix to calculate various performance metrics
  cm_train = confusionMatrix(pred_train, class_train, positive = pos_class)
  cm_test = confusionMatrix(pred_test, class_test, positive = pos_class)
  
  # ROC train and test
  roc.train = auc(multcap(response = class_train, predicted = data.matrix(response_train)))
  roc.test = auc(multcap(response = class_test, predicted = data.matrix(response_test)))
  
  info[j,'dataset'] = '6cancer'
  info[j,'method'] = method_str
  info[j,'tr_acc']=cm_train$overall["Accuracy"]
  info[j,'tst_acc']=cm_test$overall["Accuracy"]
  
  info[j,'tr_bacc']=mean(cm_train$byClass[,'Balanced Accuracy'])
  info[j,'tst_bacc']=mean(cm_test$byClass[,'Balanced Accuracy'])
  info[j,'tr_prec']=mean(cm_train$byClass[,'Precision'])
  info[j,'tst_prec']=mean(cm_test$byClass[,'Precision'])
  info[j,'tr_recall']=mean(cm_train$byClass[,'Recall'])
  info[j,'tst_recall']=mean(cm_test$byClass[,'Recall'])
  info[j,'tr_F1']=mean(cm_train$byClass[,'F1'])
  info[j,'tst_F1']=mean(cm_test$byClass[,'F1'])
  
  info[j,'tr_auroc']=roc.train
  info[j,'tst_auroc']=roc.test
}
saveRDS(info,paste(base_dir,"6cancer_bootstrap_resamples_", method_str, sfx,".RDS",sep=''))
