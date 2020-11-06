library(deepnet)
library(R.matlab)
library(pROC)
library(caret)
library(feather)


softmax <- function(x) {
    x2 = x - max(x)
    y = exp(x2)/sum(exp(x2))
}

sig <- function(x) {1/(1+exp(-x))}
predict_rbm <- function(r1,x_test,class_levels) {
    energy_mat = matrix(nrow=nrow(x_test),ncol=length(class_levels))
    for (m in 1 : length(class_levels)) {
      x_test[class_levels] = 0
      x_test[class_levels[m]] = 1
      x_test_tmp = as.matrix(x_test)
      
      h = rbm.up(r1,as.matrix(x_test_tmp))
      h[h<0.5] = 0 
      h[h>=0.5] = 1
      e_cross = diag(h %*% r1$W %*% t(x_test_tmp))
      vb_e = x_test_tmp %*% r1$B
      hb_e = h %*% r1$C
      energy_mat[,m] = hb_e + vb_e + e_cross 
    }
    prob_df = as.data.frame(t(apply(energy_mat/sqrt(nrow(x_test))/2,1,softmax)))
#    prob_df = data.frame(sig(energy_mat[,1] - energy_mat[,2]),sig(energy_mat[,2]-energy_mat[,1])) # only works for two class_levels
    names(prob_df) = class_levels
    return(prob_df)
}


cv_rbm <- function(x,y,hyper_params,pos_class) {
# do cross-fold validation with k=10
  k = 3
  flds <- createFolds(y, k = k, list = TRUE, returnTrain = TRUE)
  res = matrix(0,nrow=nrow(hyper_params),ncol=k)
  for (n in 1 : length(flds)) {
    x_train = x[flds[[n]],]
    y_train = y[flds[[n]]]
    x_test = x[-flds[[n]],]
    y_test = y[-flds[[n]]]
    for (m in 1 : nrow(hyper_params)) {
      r1<-rbm.train(as.matrix(x_train),hyper_params[m,'hidden_size'],numepochs=100,batchsize=hyper_params[m,'batch_size'],cd=hyper_params[m,'k'])
      response_test = predict_rbm(r1,x_test,levels(y_test))
      pred_test = as.factor(apply(response_test,1,which.max))
      levels(pred_test) = colnames(response_test)
      res[m,n] = mean(pred_test==y_test)     
#      cm_test = confusionMatrix(pred_test, y_test, positive = pos_class)
#      res[m,n] = =mean(cm_test$byClass['Balanced Accuracy'])
    }
  }
  # find best average performance
  mean_acc = apply(res,1,mean)
  best_idx = which.max(mean_acc)

  # rerun with best hyperparameters
  
  r1<-rbm.train(as.matrix(x),hyper_params[best_idx,'hidden_size'],numepochs=100,batchsize=hyper_params[best_idx,'batch_size'],cd=hyper_params[best_idx,'k'])
}
# Load Data -- LUM A vs. LUM B
files=c("ERpn","lumAB","ERpn_gene","lumAB_gene")
positive_class_levels = c("Positive","Luminal_A","Positive","Luminal_A")
class_levels_levels = list(c("Negative","Positive"),c("Luminal_B","Luminal_A"),c("Negative","Positive"),c("Luminal_B","Luminal_A"))
all_fracs = list(c(0.06,seq(0.1,0.95,0.05)), c(0.18,seq(0.2,0.95,0.05)),c(0.02,seq(0.05,0.95,0.05)), c(seq(0.05,0.95,0.05)))

args = commandArgs(trailingOnly=T)
#ii = 2 
#jj = 1
ii = as.numeric(args[1])
jj = as.numeric(args[2])
f = files[[ii]]
fracs = all_fracs[[ii]]
pos_class = positive_class_levels[[ii]]
class_levels = class_levels_levels[[ii]]
base_dir = paste('/boston_ailab/users/rli/quantum/Data/',f,'_splits/',sep="")

save_dir = '/boston_ailab/users/rli/quantum/Results/'

info=data.frame(frac=double(),method=character(),tr_acc=double(),tst_acc=double(),
              tr_bacc=double(),tst_bacc=double(),tr_auroc=double(), tst_auroc=double(),
              tr_prec=double(), tst_prec=double(), tr_recall=double(), tst_recall=double(),
              tr_F1=double(),tst_F1=double(),stringsAsFactors=FALSE)
n=1;
n_pc = 44
n_splits = 50
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

hyper_params =  expand.grid(k=c(1,3,5),hidden_size=seq(20,40,5),batch_size=2^(seq(3,5)))
# Run through # of Splits
for (j in 1:n_splits) {
  x_train = traindatas[[j]][,2:(n_pc+1)]
  x_test = testdatas[[j]][,2:(n_pc+1)]
  x_val = valdatas[[j]][,2:(n_pc+1)]
  x_exptest = exptestdatas[[j]][,2:(n_pc+1)]
  colnames(x_train) = paste('PC_', seq(1,ncol(x_train)), sep = '')
  colnames(x_test) = paste('PC_', seq(1,ncol(x_test)), sep = '')
  colnames(x_val) = paste('PC_', seq(1,ncol(x_val)), sep = '')
  colnames(x_exptest) = paste('PC_', seq(1,ncol(x_exptest)), sep = '')
  x_train = as.data.frame(x_train)
  x_test = as.data.frame(x_test)
  x_val = as.data.frame(x_val)
  x_exptest = as.data.frame(x_exptest)

  if (grepl("gene",files[[ii]])) {
      mins = apply(x_train,2,min)
      maxs = apply(x_train,2,max)
      x_train = sweep(x_train,2,mins,'-')
      x_train = sweep(x_train,2,maxs-mins,'/')
      x_test = sweep(x_test,2,mins,'-')
      x_test = sweep(x_test,2,maxs-mins,'/')
      x_exptest = sweep(x_exptest,2,mins,'-')
      x_test = sweep(x_test,2,mins,'-')
      
  }
  class_train = as.factor(traindatas[[j]][,1])
  class_val = as.factor(valdatas[[j]][,1])
  class_test = as.factor(testdatas[[j]][,1])
  class_exptest = as.factor(c(class_test,class_val))
  levels(class_train) = class_levels
  levels(class_test) = class_levels
  levels(class_val) = class_levels
  levels(class_exptest) = class_levels


  cat("\n---------------\nCut",as.character(j),"...\n---------------\n")

      x_train[class_levels[1]] = class_train == class_levels[1]
      x_train[class_levels[2]] = class_train == class_levels[2]
 
      r1 = cv_rbm(x_train,class_train,hyper_params,pos_class)
      
      # calculate response probabilities
      response_train = predict_rbm(r1,x_train,class_levels)
      response_test = predict_rbm(r1,x_test,class_levels)
      response_val = predict_rbm(r1,x_val,class_levels)
      response_exptest = predict_rbm(r1,x_exptest,class_levels)
      response_train[is.na(response_train)] = 0.5
      response_test[is.na(response_test)] = 0.5
      response_val[is.na(response_val)] = 0.5
      response_exptest[is.na(response_exptest)] = 0.5

      response_train_list[[n]] = response_train
      response_test_list[[n]] = response_test
      response_val_list[[n]] = response_val
      response_exptest_list[[n]] = response_exptest

      # predict labels 
      pred_train = as.factor(apply(response_train,1,which.max))
      levels(pred_train) = colnames(response_train)
      pred_test = as.factor(apply(response_test,1,which.max))
      levels(pred_test) = colnames(response_test)
      pred_val = as.factor(apply(response_val,1,which.max))
      levels(pred_val) = colnames(response_val)
      pred_exptest = as.factor(apply(response_exptest,1,which.max))
      levels(pred_exptest) = colnames(response_exptest)

      # Class Stats
      cm_train = confusionMatrix(pred_train, class_train, positive = pos_class)
      cm_test = confusionMatrix(pred_test, class_test, positive = pos_class)
      cm_val = confusionMatrix(pred_val, class_val, positive = pos_class)
      cm_exptest = confusionMatrix(pred_exptest, class_exptest, positive = pos_class)
      
      # ROC train and test
      roc.train = auc(class_train,response_train[,pos_class])
      roc.test = auc(class_test,response_test[,pos_class])
      roc.val = auc(class_val,response_val[,pos_class])
      roc.exptest = auc(class_exptest,response_exptest[,pos_class])

      info[n,'frac'] = fracs[[jj]]
      info[n,'method']='RBM'
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

save_vars = list(cm_train_list,cm_test_list,cm_val_list,cm_exptest_list,info,response_train_list,response_test_list,response_val_list,response_exptest_list)
names(save_vars) = c("cm_train_list","cm_test_list","cm_val_list","cm_exptest_list","info","response_train_list","response_test_list","response_val_list","response_exptest_list")
saveRDS(save_vars,paste(save_dir,f,"_frac_", as.character(fracs[[jj]]), "_RBM_save_binary_hidden.RDS",sep=""))
