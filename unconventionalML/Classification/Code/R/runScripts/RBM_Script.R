library(deepnet)
library(feather)
library(pROC)
library(caret)
library(R.matlab)

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
    prob_df = as.data.frame(t(apply(energy_mat/sqrt(nrow(test_data))*2,1,softmax)))
    names(prob_df) = classes
    return(prob_df)
}

cv_rbm <- function(x,y,hyper_params,pos_class,k=3) {
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
probs_list = list()

positive_classes = c("tumor","Positive","kirc","luad","Luminal_A","Luminal_A","Luminal_A")
files=c("brcaMatchedTN","ERpn","kirckirp","luadlusc","lumAB","lumAB_gene","lumAB_diffexp_gene")
lumAB_classes = c("Luminal_B","Luminal_A")
useMats = c(F,F,F,F,F,T,T)
useCombined = c(F,F,F,F,F,T,F)

save_dir = '/boston_ailab/users/rli/quantum/Results/'
hyper_params =  expand.grid(k=c(1,3,5),hidden_size=seq(20,40,5),batch_size=2^(seq(3,5)))

for (ii in 1:7) {
    dataMat = useMats[[ii]]
    combined = useCombined[[ii]]
    print(files[[ii]])
    fdir = sprintf('/boston_ailab/users/rli/quantum/Data/%s_bootstrap_resamples/',files[[ii]])
    pos_class = positive_classes[[ii]]    
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
    if (dataMat & combined) {
        mat = readMat(sprintf('/boston_ailab/users/rli/quantum/Data/%s_bootstrap_resamples.mat',files[[ii]]))
        traindatas = unlist(mat$traindatas,recursive=FALSE)
        testdatas = unlist(mat$testdatas,recursive=FALSE)
    }   
    for (n in 1 : 100) {
        
        if (dataMat & !combined) {
            mat = readMat(sprintf('%sresample_%d_data.mat',fdir,n))
            train_data = as.data.frame(mat$traindata[,-1])
            test_data = as.data.frame(mat$testdata[,-1])
            train_labels = as.factor(mat$traindata[,1])
            test_labels = as.factor(mat$testdata[,1])
            levels(train_labels) = lumAB_classes
            levels(test_labels) = lumAB_classes
            classes = lumAB_classes
        } else if (dataMat & combined) {
            train_data = as.data.frame(traindatas[[n]][,-1])
            test_data = as.data.frame(testdatas[[n]][,-1])
            train_labels = as.factor(traindatas[[n]][,1])
            test_labels = as.factor(testdatas[[n]][,1])
            levels(train_labels) = lumAB_classes
            levels(test_labels) = lumAB_classes
            classes = lumAB_classes
        } else {
            train_data = read_feather(sprintf('%sresample_%d_train.feather',fdir,n))
            test_data = read_feather(sprintf('%sresample_%d_test.feather',fdir,n))
            train_labels = read.csv(sprintf('%sresample_%d_train_labels.txt',fdir,n),header=F)[[1]]
            test_labels = read.csv(sprintf('%sresample_%d_test_labels.txt',fdir,n),header=F)[[1]]
            
            classes = levels(train_labels)
        }
        if (grepl("gene",files[[ii]])) {
            mins = apply(train_data,2,min)
            maxs = apply(train_data,2,max)
            train_data = sweep(train_data,2,mins,'-')
            train_data = sweep(train_data,2,maxs-mins,'/')
            test_data = sweep(test_data,2,mins,'-')
            test_data = sweep(test_data,2,maxs-mins,'/')
        }
        for (m in 1 : length(classes) {
            train_data[classes[m]] = train_labels == classes[m]
        } 
        r1 = cv_rbm(x_train,train_labels,hyper_params,pos_class)

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
        auc.train = auc(train_labels,response_train[,pos_class])
        auc.test = auc(test_labels,response_test[,pos_class])
        info[n,'dataset']=files[[ii]]
        info[n,'method']='RBM'
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
  
        info[n,'tr_auroc']=auc.train
        info[n,'tst_auroc']=auc.test
  
  
    }

  save_vars = list(cm_train_list,cm_test_list,info,response_train_list,response_test_list)
  names(save_vars) = c("cm_train_list","cm_test_list","info","response_train_list","response_test_list")
  saveRDS(save_vars,paste(save_dir,files[[ii]],"_rbm_save_all_binary_hidden.RDS",sep=""))
}
