# Output statistics for incremental decrease for annealing-type algorithms.

rm(list = ls())
gc()

library(R.matlab)
library(HandTill2001)
library(e1071)
#library(MLmetrics)
library(caret)

# specify the name of the dataset. Different datasets have different fractions
# depending on the size of the dataset 
d = "ERpn"

if (d == "luadlusc") {
  all_fracs = seq(0.1,0.95,0.05)  
  pos = "luad"
  classes = c("lusc","luad")
} else if (d=="lumAB") {
  low_fracs = seq(0.18, 0.2, 0.02)
  high_fracs = seq(0.25, 0.95, 0.05)
  all_fracs = c(low_fracs, high_fracs)
  pos = "Luminal_A"
  classes = c("Luminal_A","Luminal_B")
} else{
  all_fracs = c(0.06,seq(0.1,0.95,0.05))
  pos = "Positive"
  classes = c("Negative","Positive")
}

base_dir = paste('~/Dropbox-Work/Wuxi/Results/', d, '_splits/preds_for_R/',sep='')

n_splits = 50
info_list = list()
meth = 'dw' # 'dw','sa','rand',or 'field'
sfx = '_nsols20_ntotsols1000' # e.g.,'_nsols20_ntotsols1000'

for (i in 1:length(all_fracs)) {
  # load file with predictions of labels
  mat = readMat(paste(base_dir,meth,"_frac_",all_fracs[[i]],sfx,"_preds_for_R.mat",sep=""))
  
  # initialize counter
  n=1

  # initialize empty data frame with relevant fields
  info=data.frame(frac=double(),tr_acc=double(),tst_acc=double(),val_acc=double(),exptst_acc=double(),
                tr_bacc=double(),tst_bacc=double(),val_bacc=double(),exptst_bacc=double(),
                tr_auroc=double(), tst_auroc=double(),val_auroc=double(),exptst_auroc=double(),
                tr_prec=double(), tst_prec=double(), val_prec=double(), exptst_prec=double(),
                tr_recall=double(), tst_recall=double(), val_recall=double(), exptst_recall=double(),
                tr_F1=double(),tst_F1=double(),val_F1=double(),exptst_F1=double(),stringsAsFactors=FALSE)

  # loop through the number of cuts (or splits)
  for (j in 1:n_splits) {
    # get labels and treat as factors
    y_train = factor(mat$y.trains[,j])
    y_test = factor(mat$y.tests[,j])
    y_val = factor(mat$y.valids[,j])
    y_exptest = factor(mat$y.exptests[,j])
    
    # re-label factors so that matches classes
    levels(y_train) = classes
    levels(y_test) = classes
    levels(y_val) = classes
    levels(y_exptest) = classes
    
    # predict classes based on raw probabilities
    response_train = mat$y.pred.trains[,j]
    pred_train = response_train >= 0.5
    pred_train = factor(pred_train)
    levels(pred_train) = classes
    
    response_test = mat$y.pred.tests[,j]
    pred_test = response_test >= 0.5
    pred_test = factor(pred_test)
    levels(pred_test) = classes
    
    response_val = mat$y.pred.vals[,j]
    pred_val = response_val >= 0.5
    pred_val = factor(pred_val)
    levels(pred_val) = classes
    
    response_exptest = mat$y.pred.exptests[,j]
    pred_exptest = response_exptest >= 0.5
    pred_exptest = factor(pred_exptest)
    levels(pred_exptest) = classes
    
    # Confusion matrix
    cm_train = confusionMatrix(pred_train, y_train, positive = pos)
    cm_test = confusionMatrix(pred_test, y_test, positive = pos)
    cm_val = confusionMatrix(pred_val, y_val, positive = pos)
    cm_exptest = confusionMatrix(pred_exptest, y_exptest, positive= pos)

    # ROC train and test    
    roc.train = HandTill2001::auc(bincap(response = y_train, predicted = 1-response_train,true=pos))
    roc.test = HandTill2001::auc(bincap(response = y_test, predicted = 1-response_test,true=pos))
    roc.val = HandTill2001::auc(bincap(response = y_val, predicted = 1-response_val,true=pos))
    roc.exptest = HandTill2001::auc(bincap(response = y_exptest, predicted = 1-response_exptest,true=pos))
    
    # Export some summary data    
    info[n,'frac']=all_fracs[i]
    info[n,'tr_acc']=cm_train$overall["Accuracy"]
    info[n,'tst_acc']=cm_test$overall["Accuracy"]
    info[n,'val_acc']=cm_val$overall["Accuracy"]
    info[n,'exptst_acc']=cm_exptest$overall["Accuracy"]
    
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
    
    n=n+1
  }
  info_list[[i]] = info
}
saveRDS(info_list,paste(base_dir,d,"_split_pca_performance_",meth,sfx,".RDS",sep=""))
