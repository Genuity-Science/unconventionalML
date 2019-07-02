# script that creates the serial dilution cut data
# This version chooses the top 44 genes in PC1 on after each cut (i.e. the gene choice changes)
# Author: Omar Gamel - based on older script by Nicholas Cilfone

# Strictly, we should switch definition of test and validation from Nick's. 
# but we left it as is for consistency

rm(list = ls())
gc()

library(glmnet)
library(data.table)
library(methods)
library(feather)
library(caret)
library(flashpcaR)
library(dplyr)
library(stringr)
library(randomForest)

options(warn=-1)

Sys.setenv(TZ="US/Eastern")
# to catch the odd error where on rare splits cross validation doesn't work, use a heuerstic value of lambda=lambda_default in rare case of error
# used ad_hoc method due to lack of time
lambda_default=0.04
overwrite_splits=data.frame(frac=c(0.75,0.75,0.3,0.35,0.35,0.35,0.45,0.45,0.5,0.5,0.5,0.6,0.6,0.65,0.65),j=c(4,15,12,31,41,50,16,45,19,20,23,2,42,11,46))

get_svalue = function(glm){
  # if cross validation return lambda.min, else lambda_default
  svalue = ifelse(class(glm)[[1]] == "cv.glmnet","lambda.min",lambda_default)
  return(svalue)
}

get_glm_results = function(glm, data, y_levels) {
  # to catch the odd error where on rare splits cross validation doesn't work, use a heuerstic value of lambda=lambda_default in rare case of error
  class_pred = predict(glm, data, type="class", s = get_svalue(glm))
  class_pred = factor(class_pred, levels = levels(y_levels))
  return(class_pred)
}

get_stats = function(cm_frame, run, frac, algo, dataset, all_stats, auc_val) {
  stat_vec = c(frac, run, algo, dataset, cm_frame$overall['Accuracy'], cm_frame$byClass['Balanced Accuracy'], cm_frame$byClass['Sensitivity'],
               cm_frame$byClass['Specificity'], cm_frame$byClass['Precision'], cm_frame$byClass['Recall'], cm_frame$byClass['F1'], auc_val)
  names(stat_vec) = c('Frac', 'Run', 'Algorithm', 'Dataset', 'Accuracy', 'Bal.Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'Recall', 'F1', 'AUC')
  all_stats[[1]][[all_stats[[2]]]] = stat_vec
  all_stats[[2]] = all_stats[[2]] + 1
  return(all_stats)
}

#replace %subset% with 'train' or 'test'
datasetpath='~/NextCODE Health/quantumMachineLearning - Documents/rawdata/data5_lumAB_subset_normalized.feather' 
cutpath='~/NextCODE Health/quantumMachineLearning - Documents/Cuts/'
splitsfolder = '2018_07_21_09_17_LUMAB_PC_AND_GENE_MULTICUT/'
splitsfile = 'split_lists_frac_num.txt' #replace num with the fraction
save_dir=paste0(cutpath,splitsfolder)


# initialize output structure
all_stats = list(list(), 1)

### Parameters
n_genes = 44

# Split Frac Ranges
low_fracs = seq(0.18, 0.2, 0.02)
high_fracs = seq(0.30, 0.95, 0.05)
all_fracs = c(rev(high_fracs), rev(low_fracs))
## override
#all_fracs=high_fracs[-(1:8)]

data_train=read_feather(str_replace(datasetpath,"subset","train"))
data_test=read_feather(str_replace(datasetpath,"subset","test"))

X_train=as.matrix(data_train[-1])
y_train=droplevels(as.factor(data_train$cancer))
levels(y_train)=str_replace(levels(y_train)," ","_")
X_test=as.matrix(data_test[-1])
y_test=droplevels(as.factor(data_test$cancer))
levels(y_test)=str_replace(levels(y_test)," ","_")

### Generate Splits
### Main Model Loop

for (frac in all_fracs) {

  # Load the splits
  splits=read.csv(paste0(save_dir,gsub("num",as.character(frac),splitsfile)),header = TRUE,sep="\t")
  n_splits = ncol(splits)

  cat("\n\n---------------\nLoading Splits ...\n---------------")
  
  # Run through # of Splits
  for (j in 1:n_splits ) {
    cat(paste0("\nFrac ", as.character(frac), " - Split "))
    
    # ad-hoc hardcoded solution due to lack of time, keep track of problematic splits (by frac and j) and replace problematic ones with the next one 
    check_splits=data.frame(frac=c(frac),j=c(j))
    while (nrow(merge(check_splits,overwrite_splits))>0) {cat("(skip) "); j=j%%n_splits+1; check_splits$j=c(j)} 
    
    cat(paste0(as.character(j) ,": "))
    
    # Split Labels y
    y_train_j = y_train[splits[[j]]]
    y_test_j = y_train[-splits[[j]]]
    y_valid_j = y_test
    y_exp_test_j = unlist(list(y_test_j, y_valid_j)) # use expanded test for the results
    
    # Split Data X
    X_train_j = X_train[splits[[j]],]
    X_test_j = X_train[-splits[[j]],]
    X_valid_j = X_test 
    # ideally original test set should keep its name, and the part cut from training should be validation .. 
    # but kept inherited naming convention
    
    # Remove Zero Variance
    train_var = apply(X_train_j, 2, sd)
    no_var_idx = which(train_var == 0.0)
    # 5% Zeros Cutoff
    cut_frac = 0.02
    min_smp = ceiling(nrow(X_train_j) * cut_frac)
    nz_sum = colSums(X_train_j != 0)
    nz_idx = which(nz_sum < min_smp)
    no_var_idx = unique(c(no_var_idx, nz_idx))
    if (length(no_var_idx) > 0) {
      X_train_j = X_train_j[, -no_var_idx]
      X_test_j = X_test_j[, -no_var_idx]
      X_valid_j = X_valid_j[, -no_var_idx]
    }
    
    # identify top genes in pc1
    pc=flashpca(X_train_j, ndim = 1, stand = "sd", do_loadings = TRUE)
    ord=order(abs(pc$loadings[,1]),decreasing = TRUE)[1:n_genes] #rankings of top n_genes from PC1
    X_train_j=X_train_j[,ord]
    X_test_j=X_test_j[,ord]
    X_valid_j=X_valid_j[,ord]

    # X_exp_test_j = rbind(X_test_j, X_valid_j)
    
    cat("Normalize ... ")
    X_train_j_mean = apply(as.data.frame(X_train_j), 2, mean)
    X_train_j_sd = apply(X_train_j, 2, sd)
    # Train
    X_train_j_z = sweep(X_train_j, 2, X_train_j_mean, "-")
    X_train_j_z = sweep(X_train_j_z, 2, X_train_j_sd, "/")
    # Test
    X_test_j_z = sweep(X_test_j, 2, X_train_j_mean, "-")
    X_test_j_z = sweep(X_test_j_z, 2, X_train_j_sd, "/")
    # Valid
    X_valid_j_z = sweep(X_valid_j, 2, X_train_j_mean, "-")
    X_valid_j_z = sweep(X_valid_j_z, 2, X_train_j_sd, "/")
    # Exp_test: Combo Valid + Test
    X_exp_test_j_z = rbind(X_test_j_z, X_valid_j_z)
    
    ######### Model Runs
    # LASSO
    cat("lasso ... ")
    # cv can cause strange error, if it does, use simple algorithm withouyt cv
    glm_lasso = tryCatch(cv.glmnet(X_train_j_z, y_train_j, nfolds=10, alpha=1, family='binomial', standardize = FALSE),
                         error= function(err) {return(glmnet(X_train_j_z, y_train_j, alpha=1, family='binomial', standardize = FALSE))})
    lasso_cm_train = caret::confusionMatrix(get_glm_results(glm_lasso, X_train_j_z, y_train_j), y_train_j)
    lasso_cm_test = caret::confusionMatrix(get_glm_results(glm_lasso, X_test_j_z, y_train_j), y_test_j)
    lasso_cm_valid = caret::confusionMatrix(get_glm_results(glm_lasso, X_valid_j_z, y_train_j), y_valid_j)
    lasso_cm_exp_test = caret::confusionMatrix(get_glm_results(glm_lasso, X_exp_test_j_z, y_train_j), y_exp_test_j)
    lasso_auc_train = ModelMetrics::auc(as.numeric(y_train_j), 1-predict(glm_lasso, X_train_j_z, type="response", s = get_svalue(glm_lasso)))
    lasso_auc_test = ModelMetrics::auc(as.numeric(y_test_j), 1-predict(glm_lasso, X_test_j_z, type="response", s = get_svalue(glm_lasso)))
    lasso_auc_valid = ModelMetrics::auc(as.numeric(y_valid_j), 1-predict(glm_lasso, X_valid_j_z, type="response", s = get_svalue(glm_lasso)))
    lasso_auc_exp_test = ModelMetrics::auc(as.numeric(y_exp_test_j), 1-predict(glm_lasso, X_exp_test_j_z, type="response", s = get_svalue(glm_lasso)))
    # Gather All Stats
    all_stats = get_stats(lasso_cm_train, j, frac, 'LASSO', 'Train', all_stats, lasso_auc_train)
    all_stats = get_stats(lasso_cm_test, j, frac, 'LASSO', 'Test', all_stats, lasso_auc_test)
    all_stats = get_stats(lasso_cm_valid, j, frac, 'LASSO', 'Valid', all_stats, lasso_auc_valid)
    all_stats = get_stats(lasso_cm_exp_test, j, frac, 'LASSO', 'ExpandTest', all_stats, lasso_auc_exp_test)
    
    # RIDGE
    cat("ridge ... ")
    glm_ridge = tryCatch(cv.glmnet(X_train_j_z, y_train_j, nfolds=10, alpha=0, family='binomial', standardize = FALSE),
                         error= function(err) {return(glmnet(X_train_j_z, y_train_j, alpha=0, family='binomial', standardize = FALSE))})
    ridge_cm_train = caret::confusionMatrix(get_glm_results(glm_ridge, X_train_j_z, y_train_j), y_train_j)
    ridge_cm_test = caret::confusionMatrix(get_glm_results(glm_ridge, X_test_j_z, y_train_j), y_test_j)
    ridge_cm_valid = caret::confusionMatrix(get_glm_results(glm_ridge, X_valid_j_z, y_train_j), y_valid_j)
    ridge_cm_exp_test = caret::confusionMatrix(get_glm_results(glm_ridge, X_exp_test_j_z, y_train_j), y_exp_test_j)
    ridge_auc_train = ModelMetrics::auc(as.numeric(y_train_j), 1-predict(glm_ridge, X_train_j_z, type="response", s = get_svalue(glm_ridge)))
    ridge_auc_test = ModelMetrics::auc(as.numeric(y_test_j), 1-predict(glm_ridge, X_test_j_z, type="response", s = get_svalue(glm_ridge)))
    ridge_auc_valid = ModelMetrics::auc(as.numeric(y_valid_j), 1-predict(glm_ridge, X_valid_j_z, type="response", s = get_svalue(glm_ridge)))
    ridge_auc_exp_test = ModelMetrics::auc(as.numeric(y_exp_test_j), 1-predict(glm_ridge, X_exp_test_j_z, type="response", s = get_svalue(glm_ridge)))
    # Gather All Stats
    all_stats = get_stats(ridge_cm_train, j, frac, 'RIDGE', 'Train', all_stats, ridge_auc_train)
    all_stats = get_stats(ridge_cm_test, j, frac, 'RIDGE', 'Test', all_stats, ridge_auc_test)
    all_stats = get_stats(ridge_cm_valid, j, frac, 'RIDGE', 'Valid', all_stats, ridge_auc_valid)
    all_stats = get_stats(ridge_cm_exp_test, j, frac, 'RIDGE', 'ExpandTest', all_stats, ridge_auc_exp_test)
    
    # SVM
    cat("svm ... ")
    svm_ctrl = trainControl(method = "cv", number = 10)
    svmGrid = expand.grid(cost= 2^c(0:5))
    svm_mod = tryCatch(train(X_train_j_z, y_train_j, method = 'svmLinear2', trControl = svm_ctrl, tuneGrid = svmGrid, preProc = c("center", "scale"), probability = TRUE),
                       warning=function(wr) {return(train(X_train_j_z, y_train_j, method = 'svmLinear2', probability = TRUE))})
    svm_cm_train = caret::confusionMatrix(predict(svm_mod, newdata = X_train_j_z, type = 'raw'), y_train_j)
    svm_cm_test = caret::confusionMatrix(predict(svm_mod, newdata = X_test_j_z, type = 'raw'), y_test_j)
    svm_cm_valid = caret::confusionMatrix(predict(svm_mod, newdata = X_valid_j_z, type = 'raw'), y_valid_j)
    svm_cm_exp_test = caret::confusionMatrix(predict(svm_mod, newdata = X_exp_test_j_z, type = 'raw'), y_exp_test_j)
    svm_auc_train = ModelMetrics::auc(as.numeric(y_train_j), predict(svm_mod, newdata = X_train_j_z, type = 'prob')[,1])
    svm_auc_test = ModelMetrics::auc(as.numeric(y_test_j), predict(svm_mod, newdata = X_test_j_z, type = 'prob')[,1])
    svm_auc_valid = ModelMetrics::auc(as.numeric(y_valid_j), predict(svm_mod, newdata = X_valid_j_z, type = 'prob')[,1])
    svm_auc_exp_test = ModelMetrics::auc(as.numeric(y_exp_test_j), predict(svm_mod, newdata = X_exp_test_j_z, type = 'prob')[,1])
    # Gather All Stats
    all_stats = get_stats(svm_cm_train, j, frac, 'SVM', 'Train', all_stats, svm_auc_train)
    all_stats = get_stats(svm_cm_test, j, frac, 'SVM', 'Test', all_stats, svm_auc_test)
    all_stats = get_stats(svm_cm_valid, j, frac, 'SVM', 'Valid', all_stats, svm_auc_valid)
    all_stats = get_stats(svm_cm_exp_test, j, frac, 'SVM', 'ExpandTest', all_stats, svm_auc_exp_test)

    # RF
    cat("random forest ... ")
    rf_ctrl = trainControl(method = "cv", number = 10, savePredictions = "final")
    rfGrid = expand.grid(mtry=c(1:15))
    rf_mod = suppressWarnings(train(X_train_j_z, y_train_j, method = 'rf', trControl = rf_ctrl, tuneGrid = rfGrid, preProc = c("center", "scale")))
    # rf_mod = tryCatch(train(X_train_j_z, y_train_j, method = 'rf', trControl = rf_ctrl, tuneGrid = rfGrid, preProc = c("center", "scale")),
    #          warning=function(wr) {return(randomForest(X_train_j_z,y_train_j))})
    rf_cm_train = caret::confusionMatrix(rf_mod$pred[order(rf_mod$pred$rowIndex),2], y_train_j)
    rf_cm_test = caret::confusionMatrix(predict(rf_mod, newdata = X_test_j_z, type = 'raw'), y_test_j)
    rf_cm_valid = caret::confusionMatrix(predict(rf_mod, newdata = X_valid_j_z, type = 'raw'), y_valid_j)
    rf_cm_exp_test = caret::confusionMatrix(predict(rf_mod, newdata = X_exp_test_j_z, type = 'raw'), y_exp_test_j)
    rf_auc_train = ModelMetrics::auc(as.numeric(y_train_j), predict(rf_mod, newdata = X_train_j_z, type = 'prob')[,1])
    rf_auc_test = ModelMetrics::auc(as.numeric(y_test_j), predict(rf_mod, newdata = X_test_j_z, type = 'prob')[,1])
    rf_auc_valid = ModelMetrics::auc(as.numeric(y_valid_j), predict(rf_mod, newdata = X_valid_j_z, type = 'prob')[,1])
    rf_auc_exp_test = ModelMetrics::auc(as.numeric(y_exp_test_j), predict(rf_mod, newdata = X_exp_test_j_z, type = 'prob')[,1])
    # Gather All Stats
    all_stats = get_stats(rf_cm_train, j, frac, 'RF', 'Train', all_stats, rf_auc_train)
    all_stats = get_stats(rf_cm_test, j, frac, 'RF', 'Test', all_stats, rf_auc_test)
    all_stats = get_stats(rf_cm_valid, j, frac, 'RF', 'Valid', all_stats, rf_auc_valid)
    all_stats = get_stats(rf_cm_exp_test, j, frac, 'RF', 'ExpandTest', all_stats, rf_auc_exp_test)

    # NAIVE BAYES
    cat("naive bayes ... ")
    nb_ctrl = trainControl(method = "cv", number = 10)
    nb_mod = suppressWarnings(train(X_train_j_z, y_train_j, method = 'naive_bayes', trControl = nb_ctrl, preProc = c("center", "scale")))
    nb_cm_train = caret::confusionMatrix(predict(nb_mod, newdata = X_train_j_z, type = 'raw'), y_train_j)
    nb_cm_test = caret::confusionMatrix(predict(nb_mod, newdata = X_test_j_z, type = 'raw'), y_test_j)
    nb_cm_valid = caret::confusionMatrix(predict(nb_mod, newdata = X_valid_j_z, type = 'raw'), y_valid_j)
    nb_cm_exp_test = caret::confusionMatrix(predict(nb_mod, newdata = X_exp_test_j_z, type = 'raw'), y_exp_test_j)
    nb_auc_train = ModelMetrics::auc(as.numeric(y_train_j), predict(nb_mod, newdata = X_train_j_z, type = 'prob')[,1])
    nb_auc_test = ModelMetrics::auc(as.numeric(y_test_j), predict(nb_mod, newdata = X_test_j_z, type = 'prob')[,1])
    nb_auc_valid = ModelMetrics::auc(as.numeric(y_valid_j), predict(nb_mod, newdata = X_valid_j_z, type = 'prob')[,1])
    nb_auc_exp_test = ModelMetrics::auc(as.numeric(y_exp_test_j), predict(nb_mod, newdata = X_exp_test_j_z, type = 'prob')[,1])
    # Gather All Stats
    all_stats = get_stats(nb_cm_train, j, frac, 'NB', 'Train', all_stats, nb_auc_train)
    all_stats = get_stats(nb_cm_test, j, frac, 'NB', 'Test', all_stats, nb_auc_test)
    all_stats = get_stats(nb_cm_valid, j, frac, 'NB', 'Valid', all_stats, nb_auc_valid)
    all_stats = get_stats(nb_cm_exp_test, j, frac, 'NB', 'ExpandTest', all_stats, nb_auc_exp_test)
    }
  }

# Generate Output Tables
unroll_stats = all_stats[[1]]
names(unroll_stats) = paste('Iteration_', seq(1,length(unroll_stats)), sep = '')
stats_df = as.data.frame(t(bind_rows(unroll_stats)), stringsAsFactors = FALSE)
stats_df[,5:12] = sapply(stats_df[,5:12], as.numeric)
stats_df[,2] = as.numeric(stats_df[,2])
colnames(stats_df) = c('Frac', 'Run', 'Algorithm', 'Dataset', 'Accuracy', 'Bal.Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'Recall', 'F1', 'AUC')

# Write All Stats
write.table(stats_df, file = paste(save_dir, 'all_gene_stats.txt', sep=''), quote = FALSE, sep = '\t')
mean_stats = aggregate(.~Frac+Algorithm+Dataset, stats_df, mean)
std_stats = aggregate(.~Frac+Algorithm+Dataset, stats_df, sd)

# Write Mean and Std Over Cuts
write.table(mean_stats, file = paste(save_dir, 'mean_gene_stats.txt', sep=''), quote = FALSE, sep = '\t', row.names = FALSE)
write.table(std_stats, file = paste(save_dir, 'std_gene_stats.txt', sep=''), quote = FALSE, sep = '\t', row.names = FALSE)
