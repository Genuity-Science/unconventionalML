rm(list = ls())
gc()
library(feather)
library(caret)
library(stringr)
library(R.matlab)

Sys.setenv(TZ="US/Eastern")

# Load Data -- LUM A vs. LUM B
base_dir = '/Users/rli/OneDrive - NextCODE Health/Projects/quantumML/Classification/quantum/Data/'
pos_class = "Luminal_A"
classes_levels = c("Luminal_A","Luminal_B")

methods = c('glmnet','glmnet','svmLinear2','rf','nb') #('xgbLinear','svmRadial','svmLinear')
method_disp = c('lasso','ridge','svm','rf','nb')
lassoGrid=expand.grid(.alpha=1, .lambda=2^(seq(-10,10)))
ridgeGrid=expand.grid(.alpha=0, .lambda=seq(0, 100, by = 0.1))
rfGrid=expand.grid(.mtry=c(1:50))#, .ntree=c(1000, 1500, 2000))
tuneGrids=list(lassoGrid,ridgeGrid,NULL,rfGrid,NULL)
cvCtrl = trainControl(method = "cv", number = 10,classProbs = TRUE,returnResamp='all')

# Stats Mat
n_pc = 44
n_splits = 50

n=1

mat = readMat(paste(base_dir,"frac_0.4_data_resamples.mat",sep="")) # doesn't matter which frac because want original dataset
traindatas = unlist(mat$traindatas.scaled,recursive=FALSE)
exptestdatas = unlist(mat$exptestdatas.scaled,recursive=FALSE)
all_data = rbind(traindatas[[1]],exptestdatas[[1]])

x_data = all_data[,2:ncol(all_data)]
x_train_mean = apply(x_data, 2, mean)
x_train_sd = apply(x_data, 2, sd)
x_data = sweep(x_data, 2, x_train_mean, "-")
x_data = sweep(x_data, 2, x_train_sd, "/")
colnames(x_data) = paste('PC_', 1:ncol(x_data), sep = '')

y = all_data[,1]
y_classes = as.factor(y)
levels(y_classes) = classes_levels

info = data.frame()
for (m in 1:length(methods)) {
  method = methods[[m]]
  print(method_disp[[m]])
  fit = train(x_data,y_classes,method=method,trControl = cvCtrl,tuneGrid = tuneGrids[[m]])
  info[n,"method"] = method_disp[[m]]
  info[n,"rep"] = 0 #rep means repetition, 0 means the original
  d = fit$resample$Accuracy
  tmp = matrix(d,nrow=10,byrow=T)
  info[n,"mean_cv_acc"] = max(colMeans(tmp))
  n = n + 1
  for (j in 1:n_splits) {
    
    classes_tmp = sample(y_classes,length(y_classes))
    ######### Model Runs
    fit = train(x_data,classes_tmp,method=method,trControl = cvCtrl,tuneGrid = tuneGrids[[m]])
    info[n,"method"] = method_disp[[m]]
    info[n,"rep"] = j
    d = fit$resample$Accuracy
    tmp = matrix(d,nrow=10,byrow=T)
    info[n,"mean_cv_acc"] = max(colMeans(tmp))
    # Class returns the class of the prediction 
    n = n+1 
  }
}
