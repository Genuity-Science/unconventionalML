library(deepnet)
library(feather)
library(HandTill2001)
library(caret)

sig <- function(x) {1/(1+exp(-x))}

softmax <- function(x) {
    x2 = x - max(x) 
    y = exp(x2)/sum(exp(x2))
}
predict_rbm <- function(r1,test_data,classes) {
    energy_mat = matrix(nrow=nrow(test_data),ncol=length(classes))
    for (m in 1 : length(classes)) {
      test_data[classes] = 0
      test_data[classes[m]] = 1
      test_data_tmp = as.matrix(test_data)
      
      h = rbm.up(r1,as.matrix(test_data_tmp))
      e_cross = diag(h %*% r1$W %*% t(test_data_tmp))
      vb_e = test_data_tmp %*% r1$B
      hb_e = h %*% r1$C
      energy_mat[,m] = hb_e + vb_e + e_cross 
    }
 #   prob_df = data.frame(sig(energy_mat[,1] - energy_mat[,2]),sig(energy_mat[,2]-energy_mat[,1])) # only works for two classes
    prob_df = as.data.frame(t(apply(energy_mat/sqrt(nrow(test_data))/2,1,softmax)))
    names(prob_df) = classes
    return(prob_df)
}

probs_list = list()
ii=1
files='6cancer'

save_dir = '/boston_ailab/users/rli/quantum/Results/'
    print(files[[ii]])
    fdir = sprintf('/boston_ailab/users/rli/quantum/Data/%s_bootstrap_resamples/',files[[ii]])
    pos_class = NULL  
    info=data.frame(dataset=character(),method=character(),tr_acc=double(),tst_acc=double(),
                    tr_bacc=double(),tst_bacc=double(),tr_auroc=double(), tst_auroc=double(),
                    tr_prec=double(), tst_prec=double(), tr_recall=double(), tst_recall=double(),
                    tr_F1=double(),tst_F1=double(),stringsAsFactors=FALSE)
    # initialize list of confusion matrices
    cm_train_list = list()
    cm_test_list = list()
    
    # initialize list of response probabilities
    response_train_list = list()
    response_test_list = list()
    for (n in 1 : 100) {
        train_data = read_feather(sprintf('%sresample_%d_train.feather',fdir,n))
        test_data = read_feather(sprintf('%sresample_%d_test.feather',fdir,n))
        train_labels = read.csv(sprintf('%sresample_%d_train_labels.txt',fdir,n),header=F)[[1]]
        test_labels = read.csv(sprintf('%sresample_%d_test_labels.txt',fdir,n),header=F)[[1]]
        
        classes = levels(train_labels)
    
        for (m in 1 : length(classes)) {    
            train_data[classes[m]] = train_labels == classes[m]
        }
        r1<-rbm.train(as.matrix(train_data),30,numepochs=100,batchsize=8,cd=1)

        # calculate response probabilities
        response_train = predict_rbm(r1,train_data,classes)
        response_test = predict_rbm(r1,test_data,classes)
        response_train_list[[n]] = response_train
        response_test_list[[n]] = response_test
        
        pred_train = as.factor(apply(response_train,1,which.max))
        levels(pred_train) = colnames(response_train)
        pred_test = as.factor(apply(response_test,1,which.max))
        levels(pred_test) = colnames(response_test)
        cm_train = confusionMatrix(pred_train, train_labels, positive = pos_class)
        cm_test = confusionMatrix(pred_test, test_labels, positive = pos_class)
  
        # ROC train and test
        auc.train = auc(multcap(response = train_labels, predicted = data.matrix(response_train)))
        auc.test = auc(multcap(response = test_labels, predicted = data.matrix(response_test)))
        info[n,'dataset']=files[[ii]]
        info[n,'method']='RBM'
        info[n,'tr_acc']=cm_train$overall["Accuracy"]
        info[n,'tst_acc']=cm_test$overall["Accuracy"]
        cm_train_list[[n]] = cm_train
        cm_test_list[[n]] = cm_test
  
        info[n,'tr_bacc']=mean(cm_train$byClass[,'Balanced Accuracy'])
        info[n,'tst_bacc']=mean(cm_test$byClass[,'Balanced Accuracy'])
        info[n,'tr_prec']=mean(cm_train$byClass[,'Precision'])
        info[n,'tst_prec']=mean(cm_test$byClass[,'Precision'])
        info[n,'tr_recall']=mean(cm_train$byClass[,'Recall'])
        info[n,'tst_recall']=mean(cm_test$byClass[,'Recall'])
        info[n,'tr_F1']=mean(cm_train$byClass[,'F1'])
        info[n,'tst_F1']=mean(cm_test$byClass[,'F1'])
  
        info[n,'tr_auroc']=auc.train
        info[n,'tst_auroc']=auc.test
  
  
    }

  save_vars = list(cm_train_list,cm_test_list,info,response_train_list,response_test_list)
  names(save_vars) = c("cm_train_list","cm_test_list","info","response_train_list","response_test_list")
  saveRDS(save_vars,paste(save_dir,files[[ii]],"_rbm_save_all.RDS",sep=""))
