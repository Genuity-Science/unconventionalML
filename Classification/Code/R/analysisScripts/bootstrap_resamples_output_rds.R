# Output RDS file with metrics for bootstrap resamples. Run after
# doing analysis code in Matlab. Used to be consistent in how the performance
# metrics are calculated for the classical algorithms. Saves as .RDS file
# Author: Richard Li
rm(list=ls())
gc()

library(caret)
library(HandTill2001)
library(R.matlab)

datasets = c("brcaMatchedTN","ERpn","kirckirp","luadlusc","lumAB","lumAB_gene")

# change the string to change method. Assumes that all saved in the same format; i.e.,
# ${method_str}_nsols20_pred_for_R.mat, where method_str can be "dw","sa","rand". 
# Any user-defined string is fine, as long as as this fileformat
meth = "field" #'dw' 'sa' 'rand' 'field'
sfx = "" # '_nsols_20_ntotsols_1000'

base_dir = '~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/'
all_info = data.frame()
# define positive classes and levels
positive_classes = c("tumor","Positive","kirc","luad","Luminal_A","Luminal_A")
classes_levels = list(c("normal","tumor"),c("Negative","Positive"),c("kirp","kirc"),c("lusc","luad"),c("Luminal_B","Luminal_A"),c("Luminal_B","Luminal_A"))
n_splits = 100

for (n in 1:length(datasets)) {
  dir=paste(base_dir,datasets[[n]],"_bootstrap_resamples/",sep='')
  mat = readMat(paste(dir,meth,sfx,'_pred_for_R.mat',sep=""))
  pos_class=positive_classes[[n]]
  classes = classes_levels[[n]]
  
  info=data.frame(dataset=character(),method=character(),tr_acc=double(),tst_acc=double(),
                  tr_bacc=double(),tst_bacc=double(),tr_auroc=double(), tst_auroc=double(),
                  tr_prec=double(), tst_prec=double(), tr_recall=double(), tst_recall=double(),
                  tr_F1=double(),tst_F1=double(),stringsAsFactors=FALSE)
  
  for (j in 1:n_splits) {
    
    response_train = mat$y.pred.trains[,j]
    response_test = mat$y.pred.tests[,j]
    class_train = factor(mat$y.trains[,j])
    class_test = factor(mat$y.tests[,j])
    
    pred_train = factor(response_train >= 0.5)
    pred_test = factor(response_test >= 0.5)
    levels(pred_train) = classes
    levels(pred_test) = classes
    levels(class_train) = classes
    levels(class_test) = classes
    
    cm_train = confusionMatrix(pred_train, class_train, positive = pos_class)
    cm_test = confusionMatrix(pred_test, class_test, positive = pos_class)
    
    if(pos_class == classes[[1]]) {
      auc_train_pred = 1-response_train
      auc_test_pred = 1-response_test
    }  else {
      auc_train_pred = response_train
      auc_test_pred = response_test
    }
    # ROC train and test
    roc.train = auc(bincap(response = class_train, predicted = auc_train_pred,true=pos_class))
    roc.test = auc(bincap(response = class_test, predicted = auc_test_pred,true=pos_class))
    
    info[j,'dataset'] = datasets[[n]]
    info[j,'method']=meth
    info[j,'tr_acc']=cm_train$overall["Accuracy"]
    info[j,'tst_acc']=cm_test$overall["Accuracy"]
    info[j,'tr_bacc']=mean(cm_train$byClass['Balanced Accuracy'])
    info[j,'tst_bacc']=mean(cm_test$byClass['Balanced Accuracy'])
    info[j,'tr_prec']=mean(cm_train$byClass['Precision'])
    info[j,'tst_prec']=mean(cm_test$byClass['Precision'])
    info[j,'tr_recall']=mean(cm_train$byClass['Recall'])
    info[j,'tst_recall']=mean(cm_test$byClass['Recall'])
    info[j,'tr_F1']=mean(cm_train$byClass['F1'])
    info[j,'tst_F1']=mean(cm_test$byClass['F1'])
    info[j,'tr_auroc']=roc.train
    info[j,'tst_auroc']=roc.test
  }
  all_info = rbind(all_info,info)
}
saveRDS(all_info,paste(base_dir,"bootstrap_resamples_",meth,sfx,".RDS",sep=""))
