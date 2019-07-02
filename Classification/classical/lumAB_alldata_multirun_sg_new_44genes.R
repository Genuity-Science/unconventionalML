# Script to run classical classification algorithms on top 44 genes of PC1 for LumA vs. LumB binomial comparison.
rm(list = ls())
gc()
library(feather)
library(caret)
library(stringr)
library(HandTill2001)
library(MLmetrics)
library(PRROC)
library(flashpcaR)
library(dplyr)
library(data.table)


library(parallel)
library(doParallel)

cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)

Sys.setenv(TZ="US/Eastern")
intervalStart <- Sys.time()


# Load Data -- LUM A vs. LUM B
train_data = read_feather('data5_lumAB_train_normalized.feather')
test_data = read_feather('data5_lumAB_test_normalized.feather')
all_data = rbind(train_data,test_data)

# read top 44 gene from PC1 on train data
top44genes <- read.table("pc1_top_44_genes_ids.txt", header=F)
cols2sel <- c(as.character(top44genes$V1))

# block_path = '/Users/ncilfone/Documents/Projects/quantum_stuff/data/classical_comp/results/'
block_path = 'lumAB_multicut_multirun/'
all_labels = all_data[,1]
y_data = all_labels[[1]]
y_data = gsub(" ", "_", y_data, fixed = TRUE)
y_data = as.factor(y_data)
x_data = all_data[,-1]
run_name = "_LUMAB"

# Run name
Sys.chmod(block_path, "777")
time_path = paste(format(Sys.time(), "%Y_%m_%d"), run_name, "_ALL_MULTICUT/", sep="")
save_dir = paste(block_path, time_path, sep="")
dir.create(save_dir)

### Parameters
# Num PCs
n_pc = 44
n_splits = 100
n_rows = nrow(all_data)
n_tr = nrow(train_data)
pos_class = "Luminal_A"

### stuff for cv
methods = c('glmnet','glmnet','svmLinear2','rf','nb') #('xgbLinear','svmRadial','svmLinear')
method_disp = c('lasso','ridge','svm','rf','nb')
lassoGrid=expand.grid(.alpha=1, .lambda=seq(0, 100, by = 0.1))
ridgeGrid=expand.grid(.alpha=0, .lambda=seq(0, 100, by = 0.1))
rfGrid=expand.grid(.mtry=c(1:44))
tuneGrids=list(lassoGrid,ridgeGrid,NULL,rfGrid,NULL)
cvCtrl = trainControl(method = "cv", number = 10,classProbs = TRUE,allowParallel = TRUE)

# Stats Mat
info=data.frame(method=character(),tr_acc=double(),tst_acc=double(),
                tr_bacc=double(),tst_bacc=double(),tr_auroc=double(), tst_auroc=double(),
                tr_prec=double(), tst_prec=double(), tr_recall=double(), tst_recall=double(),
                tr_F1=double(),tst_F1=double(),stringsAsFactors=FALSE)

### Baseline Models
# Remove Zero Variance
data_var = apply(all_data, 2, sd)
no_var_idx = which(data_var == 0.0)
if (length(no_var_idx) > 0) {
  all_data = all_data[, -no_var_idx]
}

### Generate Splits
### Main Model Loop
cat(paste('\n', format(Sys.time(), "%H:%M"), '\n', sep=""))

# # Split them up
# splits = matrix(0,nrow=n_splits,ncol=n_tr)
# for (n in 1:n_splits) {
#   splits[n,] = sort(sample.int(n_rows,n_tr))
# }
# # Save the splits
# cat("\n---------------\nSave Splits...\n---------------\n")
# write.table(splits,"lumab_all_data_multirun_traincuts.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

# splits = as.matrix(read.table("~/Dropbox-Work/Wuxi/Data/lumab_all_data_multirun_traincuts.txt"))
# index for stats
n=1;
gene_list = list()
cm_train_list = list()
cm_test_list = list()
response_train_list = list()
response_test_list = list()
splits = as.matrix(read.table("bootstrap_resamples/lumAB_bootstrap_resamples.txt"))

# Run through # of Splits
pc1_top44genes_100cuts <- list()
for (j in 1:n_splits) {
  cat(paste("\nSplit ", as.character(j) ,"\n", sep=''))
  # Split Labels
  y_train = y_data[splits[j,]]
  y_test = y_data[-splits[j,]]
  # Split Data
  data_train = x_data[splits[j,],]
  data_test = x_data[-splits[j,],]
  # Remove Zero Variance
  train_var = apply(data_train, 2, sd)
  no_var_idx = which(train_var == 0.0)
  if (length(no_var_idx) > 0) {
    data_train = data_train[, -no_var_idx]
    data_test = data_test[, -no_var_idx]
  }
  # Normalize All Data
  cat("\n---------------\nNormalize Data...\n---------------\n")
  train_mean = apply(data_train, 2, mean)
  train_sd = apply(data_train, 2, sd)
  # Train
  x_train = sweep(data_train, 2, train_mean, "-")
  x_train = sweep(x_train, 2, train_sd, "/")
  # Test
  x_test = sweep(data_test, 2, train_mean, "-")
  x_test = sweep(x_test, 2, train_sd, "/")
  
  # get subset for train and test
  x_train_44genes <- x_train %>% select(cols2sel)
  x_train_44genes <- setcolorder(x_train_44genes,cols2sel)
  
  x_test_44genes <- x_test %>% select(cols2sel)
  x_test_44genes <- setcolorder(x_test_44genes,cols2sel)
  
  cat("\n---------------\nNormalize PC Data...\n---------------\n")
  x_train_44genes_mean = apply(x_train_44genes, 2, mean)
  x_train_44genes_sd = apply(x_train_44genes, 2, sd)
  # Train
  z_x_train_44genes = sweep(x_train_44genes, 2, x_train_44genes_mean, "-")
  z_x_train_44genes = sweep(z_x_train_44genes, 2, x_train_44genes_sd, "/")
  # Test
  z_x_test_44genes = sweep(x_test_44genes, 2, x_train_44genes_mean, "-")
  z_x_test_44genes = sweep(z_x_test_44genes, 2, x_train_44genes_sd, "/")
  
  
  # Fix Data Types
  y_train = as.factor(y_train)
  y_test = factor(y_test, levels = levels(y_train))
  
  
  for (m in 1:length(methods)) {
    method=methods[m]
    print(method_disp[[m]])
    
    ######### Model Runs
    fit = train(z_x_train_44genes,y_train,method=method,trControl = cvCtrl,tuneGrid = tuneGrids[[m]])
    
    # save the model
    save_name = paste(save_dir, 'glm_run', '_', as.character(j), '_' , method_disp[[m]], sep = '')
    save(fit, file=save_name)
 
    
    fit_var_imp = varImp(fit)
    fit_var_imp = as.data.frame(fit_var_imp$importance) #[,1])
    colnames(fit_var_imp) <-  paste('features_run', '_', as.character(j), '_' , method_disp[[m]], sep = '')
    save_name = paste(save_dir, 'features_run', '_', as.character(j), '_' , method_disp[[m]], sep = '')
    write.table(fit_var_imp, file = save_name, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = '\t')
    
    # Class returns the class of the prediction
    pred_train = predict(fit,z_x_train_44genes)
    pred_test = predict(fit,z_x_test_44genes)
    
    # Class Stats
    cm_train = confusionMatrix(pred_train, y_train, positive = pos_class)
    cm_test = confusionMatrix(pred_test, y_test, positive = pos_class)
    
    # calculate response probabilities
    response_train = predict(fit,z_x_train_44genes,type="prob")
    response_test = predict(fit,z_x_test_44genes,type="prob")
    response_train_list[[n]] = response_train
    response_test_list[[n]] = response_test
    
    # AUC train and test
    auc.train = ModelMetrics::auc(as.numeric(y_train), predict(fit, newdata = z_x_train_44genes, type = 'prob')[,1])
    auc.test = ModelMetrics::auc(as.numeric(y_test), predict(fit, newdata = z_x_test_44genes, type = 'prob')[,1])
    #roc.train = auc(multcap(response = y_train, predicted = data.matrix(response_train_list[[n]])))
    #roc.test = auc(multcap(response = y_test, predicted = data.matrix(response_test)))
    
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
    # info[n,'tr_auroc']=roc.train
    # info[n,'tst_auroc']=roc.test
    n = n+1
  }
}
save_vars = list(cm_train_list,cm_test_list,info,response_train_list,response_test_list)
names(save_vars) = c("cm_train_list","cm_test_list","info","response_train_list","response_test_list")
saveRDS(save_vars,paste0(save_dir,"/lumab_multicut_multirun_44genes_save.RDS"))
#save(gene_list, file = paste(save_dir, 'final_gene_list', sep = ''))

#stopCluster(cluster)
intervalEnd <- Sys.time()
paste("100 iterations of lumAB data took",intervalEnd - intervalStart,attr(intervalEnd - intervalStart,"units"))

# average scores for variable importance for 100 runs
methods = c('glmnet','glmnet','svmLinear2','rf','nb') #('xgbLinear','svmRadial','svmLinear')
method_disp = c('lasso','ridge','svm','rf','nb')
for (m in 1:length(methods)) {
  method=methods[m]
  print(method_disp[[m]])
 
  files <- list.files("./", pattern=paste0("features_run_[0-9]*_",as.character(method_disp[[m]])))
  dataMerge <- data.frame()
  for(f in files){ 
    ReadInMerge <- read.table(file=f, header=T, na.strings="NULL")
    dataMerge <- transform(merge(dataMerge, ReadInMerge, by=0, all=T), row.names=Row.names, Row.names=NULL)
    
  }
  write.csv(dataMerge,paste0("features_runs_all_100_",as.character(method_disp[[m]]),".csv"), row.names=T, quote = F)
  dataMerge.avg <- as.data.frame(rowMeans(dataMerge))
  dataMerge.avg <- dataMerge.avg[order(-dataMerge.avg$`rowMeans(dataMerge)`), , drop = FALSE]
  
  dataMerge.avg$ranking <- c(1:44)
  colnames(dataMerge.avg) <- c(paste0(method_disp[[m]],"_avg_score_100runs"),paste0(method_disp[[m]],"_ranking"))
  write.csv(dataMerge.avg,paste0("features_runs_avg_100_",as.character(method_disp[[m]]),".csv"), row.names=T, quote = F)
}


# combine rankings for classical and quantum
setwd("./lumAB_multicut_multirun/2019_05_01_LUMAB_ALL_MULTICUT")
cl_files <- list.files("./", pattern="features_runs_avg_100*")
qn_files <- list.files("./", pattern="*ranking*") #quantum and pc1 loadings

files <- c(cl_files,qn_files )
dataMerge <- read.csv("gene_symbol_mapping.csv", header=T, row.names = 1)
for(f in files){ 
  ReadInMerge <- read.csv(file=f, header=T, na.strings="NULL",row.names=1)
  dataMerge <- transform(merge(dataMerge, ReadInMerge, by=0, all=T), row.names=Row.names, Row.names=NULL)
  
}
dataMerge <- dataMerge[order(dataMerge$lasso_ranking), , drop = FALSE]
write.csv(dataMerge,"lumAB_44genes_100runs_classical_quantum_pc1_ranking.csv", row.names=T, quote = F)


               
