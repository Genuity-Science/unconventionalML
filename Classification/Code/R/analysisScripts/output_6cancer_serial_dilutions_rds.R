# Output RDS file with annealing-type metrics for incremental decrease of data 
# for multinomial data. Run after doing analysis code in Matlab. Used to be consistent 
# in how the performance metrics are calculated for the classical algorithms. 
# Saves as .RDS file
# Author: Richard Li

rm(list=ls())
gc()

library(caret)
library(HandTill2001)
library(R.matlab)

n_splits = 50
base_dir = '~/Dropbox-Work/Wuxi/Results/6cancer_splits/preds_for_R/'

# change the string to change method. Assumes that all saved in the same format; e.g.,
# ${method_str}_nsols20_pred_for_R.mat, where method_str can be "dw","sa","rand". 
# Any user-defined string is fine, as long as as this fileformat
d = '6cancer'
meth = 'sa' # 'dw','sa', 'field', or 'rand'.
sfx = '_nsols20_ntotsols1000_b1_0d3'

# specify the fractions of data to use 
all_fracs = c(seq(0.02,0.2,0.02),seq(0.25,0.95,0.05))

info_list = list()
pos_class=NULL
classes = c("lihc","brca","lgg","coad",'kidn','lung')

for (i in 1:length(all_fracs)) {
    info=data.frame(frac=double(),tr_acc=double(),tst_acc=double(),
                    tr_bacc=double(),tst_bacc=double(),tr_auroc=double(), tst_auroc=double(),
                    tr_prec=double(), tst_prec=double(), tr_recall=double(), tst_recall=double(),
                    tr_F1=double(),tst_F1=double(),stringsAsFactors=FALSE)

    mat = readMat(paste(base_dir,meth,"_frac_",all_fracs[[i]],sfx,"_preds_for_R.mat",sep=""))
    response_trains = unlist(mat$y.pred.trains,recursive=FALSE)
    response_tests = unlist(mat$y.pred.tests,recursive=FALSE)
    response_exptests = unlist(mat$y.pred.exptests,recursive=FALSE)
    response_vals = unlist(mat$y.pred.vals,recursive=FALSE)
    
    for (j in 1:n_splits) {
 
      # response_train and response_test are the predicted probabilities for each 
      # class. Of size Nx6, where 6 is the number of classes
      response_train = response_trains[[j]]
      response_test = response_tests[[j]]
      response_val = response_vals[[j]]
      response_exptest = response_exptests[[j]]

      # class test and class train are the actual classes
      class_test = factor(mat$y.tests[,j])
      class_train = factor(mat$y.trains[,j])
      class_val = factor(mat$y.valids[,j])
      class_exptest = factor(mat$y.exptests[,j])

      pred_train = factor(max.col(response_train))
      pred_test = factor(max.col(response_test))
      pred_val = factor(max.col(response_val))
      pred_exptest = factor(max.col(response_exptest))

      if (nlevels(pred_train) != 6) {
        print(paste("Frac ", as.character(all_fracs[[i]]), " Pred train iteration ", j, " does not have enough levels"))
      }
      if (nlevels(pred_test) != 6) {
        print(paste("Frac ", as.character(all_fracs[[i]]), " Pred test iteration ", j, " does not have enough levels"))
      }
      if (nlevels(pred_val) != 6) {
        print(paste("Frac ", as.character(all_fracs[[i]]), " Pred val iteration ", j, " does not have enough levels"))
      }
      if (nlevels(pred_exptest) != 6) {
        print(paste("Frac ", as.character(all_fracs[[i]]), " Pred exptest iteration ", j, " does not have enough levels"))
      }
       
      colnames(response_train) = classes
      colnames(response_test) = classes
      colnames(response_val) = classes
      colnames(response_exptest) = classes

      levels(pred_train) = classes
      levels(pred_test) = classes
      levels(pred_exptest) = classes
      levels(pred_val) = classes
      
      levels(class_train) = classes
      levels(class_test) = classes
      levels(class_val) = classes
      levels(class_exptest) = classes
  
      # confusion matrix to calculate various performance metrics
      cm_train = confusionMatrix(pred_train, class_train, positive = pos_class)
      cm_test = confusionMatrix(pred_test, class_test, positive = pos_class)
      cm_val = confusionMatrix(pred_val, class_val, positive = pos_class)
      cm_exptest = confusionMatrix(pred_exptest, class_exptest, positive = pos_class)
  
      # ROC train and test
      roc.train = auc(multcap(response = class_train, predicted = data.matrix(response_train)))
      roc.test = auc(multcap(response = class_test, predicted = data.matrix(response_test)))
      roc.val = auc(multcap(response = class_val, predicted = data.matrix(response_val)))
      roc.exptest = auc(multcap(response = class_exptest, predicted = data.matrix(response_exptest)))
  
      info[j,'frac'] = all_fracs[[i]]
      info[j,'tr_acc']=cm_train$overall["Accuracy"]
      info[j,'tst_acc']=cm_test$overall["Accuracy"]
      info[j,'val_acc']=cm_val$overall["Accuracy"]
      info[j,'exptst_acc']=cm_exptest$overall["Accuracy"]

      info[j,'tr_bacc']=mean(cm_train$byClass[,'Balanced Accuracy'])
      info[j,'tst_bacc']=mean(cm_test$byClass[,'Balanced Accuracy'])
      info[j,'val_bacc']=mean(cm_val$byClass[,'Balanced Accuracy'])
      info[j,'exptst_bacc']=mean(cm_exptest$byClass[,'Balanced Accuracy'])

      info[j,'tr_prec']=mean(cm_train$byClass[,'Precision'])
      info[j,'tst_prec']=mean(cm_test$byClass[,'Precision'])
      info[j,'val_prec']=mean(cm_val$byClass[,'Precision'])
      info[j,'exptst_prec']=mean(cm_exptest$byClass[,'Precision'])

      info[j,'tr_recall']=mean(cm_train$byClass[,'Recall'])
      info[j,'tst_recall']=mean(cm_test$byClass[,'Recall'])
      info[j,'val_recall']=mean(cm_val$byClass[,'Recall'])
      info[j,'exptst_recall']=mean(cm_exptest$byClass[,'Recall'])

      info[j,'tr_F1']=mean(cm_train$byClass[,'F1'])
      info[j,'tst_F1']=mean(cm_test$byClass[,'F1'])
      info[j,'val_F1']=mean(cm_val$byClass[,'F1'])
      info[j,'exptst_F1']=mean(cm_exptest$byClass[,'F1'])

      info[j,'tr_auroc']=roc.train
      info[j,'tst_auroc']=roc.test
      info[j,'val_auroc']=roc.val
      info[j,'exptst_auroc']=roc.exptest

    }
    info_list[[i]] = info
}
saveRDS(info_list,paste(base_dir,d,"_split_pca_performance_",meth,sfx,".RDS",sep=""))
