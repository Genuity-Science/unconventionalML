# Output .mat file with multinomial cuts for DW. 
# Example usage: Rscript 6cancer_output_serial_dilutions_pca.R frac j
# Did this because R has pretty poor memory management and would run out of 
# space. This way only handle one at a time and can do the pca

library(methods)
library(feather)
library(caret)
library(flashpcaR)
library(dplyr)
library(R.matlab)

Sys.setenv(TZ="US/Eastern")

# function that switches how do PCA, dependent on number of features. Takes care of standard scaling
switch_pca = function(x_train, n_pc) {
  # PC Model
  if (nrow(x_train) >= ((2 * n_pc)+1)) {
    pc_model = flashpca(as.matrix(x_train), ndim = n_pc, stand = "sd", do_loadings = TRUE)
    pc_flag = 'flash'
    cat("\nUsing FlashPCA...\n")
  }
  else {
    pc_model = prcomp(as.matrix(x_train), center = TRUE, scale. = TRUE)
    pc_flag = 'base'
    cat("\nFalling Back to PRCOMP...\n")
  }
  return(list(pc_model, pc_flag))
}

predict_correct_pca = function(pc_res, x, n_pc) {
  pc_model = pc_res[[1]]
  if (pc_res[[2]] == 'flash') {
    x_proj = project(as.matrix(x), loadings = pc_model$loadings,  orig_mean = pc_model$center, orig_sd = pc_model$scale)
    x_pred = x_proj$projection
  }
  else {
    x_proj = predict(pc_model, newdata=as.matrix(x))
    x_pred = x_proj[,1:n_pc]/sqrt(dim(x)[2])
  }
  return(x_pred)
}


args = commandArgs(trailingOnly=TRUE)
i = as.numeric(args[1])
j = as.numeric(args[2])

# load data 
train_data = read_feather("/home/richard/Dropbox-Work/Wuxi/Data/data5_6cancer_train_normalized.feather")
test_data = read_feather("/home/richard/Dropbox-Work/Wuxi/Data/data5_6cancer_test_normalized.feather")

# identify dataset and load 
train_labels = as.factor(gsub(" ","_",train_data$cancer,fixed=TRUE))
test_labels = as.factor(gsub(" ","_",test_data$cancer,fixed=TRUE))
train_data = train_data[,-1]
test_data = test_data[,-1]

# Run name
block_path = '~/Dropbox-Work/Wuxi/Data/6cancer_splits/'
dir.create(block_path)


# n_pc = 44
n_pc = 13
n_splits = 50

#low_fracs = seq(0.02, 0.2, 0.02)
#high_fracs = seq(0.25, 0.95, 0.05)
#all_fracs = c(low_fracs,high_fracs)
all_fracs = seq(0.05,0.95,0.05)

#min_frac = (n_pc + 1) / nrow(train_data)
# # Remove those that eclipse the min threshold
#rm_min_frac = which(all_fracs < min_frac)
#if (length(rm_min_frac) > 0) {
#  all_fracs = all_fracs[-rm_min_frac]
#}

#for (i in 1:length(all_fracs)) {
  cat(paste('\n', format(Sys.time(), "%H:%M"), '\n', sep=""))
  
  # Split them up - if doing for the first time
  # splits = createDataPartition(train_labels, times = n_splits, p = all_fracs[i])
  # # Save the splits
  # cat("\n---------------\nSave Splits...\n---------------\n")
   save_dir = paste(block_path, "split_lists_frac_", as.character(all_fracs[i]),'_dir/',sep="")
   split_name = paste(block_path, "split_lists_frac_", as.character(all_fracs[i]),'.txt',sep="")
  # write.table(split_name,x=as.data.frame(splits),quote=F,row.names=F,sep='\t')
  # 
  splits=read.table(split_name,header=TRUE)
  
  # Run through # of Splits
#  for (j in 1:n_splits) {
    # load(split_name)
    
    cat(paste("\nFrac ", as.character(all_fracs[i]), " - Split ", as.character(j) ,"\n", sep=''))
    # Split Labels
    class_train = train_labels[splits[[j]]]
    class_test = train_labels[-splits[[j]]]
    class_valid = test_labels
    class_exp_test = c(class_test, class_valid)
    # Split Data
    data_train = train_data[splits[[j]],]
    data_test = train_data[-splits[[j]],]
    # Remove Zero Variance
    train_var = apply(data_train, 2, sd)
    no_var_idx = which(train_var == 0.0)
    if (length(no_var_idx) > 0) {
      data_train = data_train[, -no_var_idx]
      data_test = data_test[, -no_var_idx]
      data_valid = test_data[, -no_var_idx]
    }
    # Normalize All Data
    cat("\n---------------\nNormalize Data...\n---------------\n")
    x_train = data_train
    x_test = data_test
    x_valid = data_valid
    # Combo Valid + Test
    x_exp_test = rbind(x_test, x_valid)
    
    cat("\n---------------\nPrincipal Components...\n---------------\n")
    # PC Model
    pc_res = switch_pca(x_train, n_pc = n_pc)
    # Train
    z_pc_train = predict_correct_pca(pc_res, x_train, n_pc)
    # Test
    z_pc_test = predict_correct_pca(pc_res, x_test, n_pc)
    # Valid
    z_pc_valid = predict_correct_pca(pc_res, x_valid, n_pc)
    # Exp Test
    z_pc_exp_test = predict_correct_pca(pc_res, x_exp_test, n_pc)
    
    # Rename Cols
    colnames(z_pc_train) = paste('PC_', seq(1,ncol(z_pc_train)), sep = '')
    colnames(z_pc_test) = paste('PC_', seq(1,ncol(z_pc_test)), sep = '')
    colnames(z_pc_valid) = paste('PC_', seq(1,ncol(z_pc_valid)), sep = '')
    colnames(z_pc_exp_test) = paste('PC_', seq(1,ncol(z_pc_exp_test)), sep = '')
    z_pc_exp_test = cbind(class_exp_test-1,z_pc_exp_test)
    z_pc_train = cbind(as.numeric(class_train)-1,z_pc_train)
    z_pc_test = cbind(as.numeric(class_test)-1,z_pc_test)
    z_pc_valid = cbind(as.numeric(class_valid)-1,z_pc_valid)
    writeMat(paste(save_dir,"resample_", as.character(j), "_data.mat", sep=''),
             testdata=as.matrix(z_pc_test),valdata = data.matrix(z_pc_valid),
             traindata=as.matrix(z_pc_train),exptest=as.matrix(z_pc_exp_test))
    gc()
#  }
#}
