# Output .mat files with resample data. Used for the binomial datasets.
# Outputting for 6cancer dataset took too much memory; see output_6cancer_bootstrap_resamples_pca.R

library(feather)
library(R.matlab)
library(flashpcaR)

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
    x_pred = x_proj[,1:n_pc]
  }
  return(x_pred)
}
# 6cancer took too much memory. See output_6cancer_bootstrap_resamples_for_DW.R
files=c("brcaMatchedTN","ERpn","kirckirp","luadlusc","lumAB")#,"6cancer") "lumAB"
base_dir = '~/Dropbox-Work/Wuxi/Data/'
positive_classes = c("tumor","Positive","kirc","luad","Luminal_A")

for (n in 1:length(files)) {
  gc()
  f = files[n]
  cat(paste("Dataset: ", f ,"\n", sep=''))
  cat("\n---------------\nLoading train data...\n---------------\n")
  train_data = read_feather(paste(base_dir,"data5_",f,"_train_normalized.feather",sep=""))

  cat("\n---------------\nLoading test data...\n---------------\n")
  test_data = read_feather(paste(base_dir,"data5_",f,"_test_normalized.feather",sep=""))

  all_data = rbind(train_data,test_data)
  block_path = '~/Dropbox-Work/Wuxi/Data/'

  # fix factor names
  all_labels = all_data[,1]
  y_data = all_labels[[1]]
  y_data = gsub(" ", "_", y_data, fixed = TRUE)
  y_data = as.factor(y_data)
  x_data = all_data[,-1]
  run_name = paste("_",f,sep="")
  
  # Run name
  Sys.chmod(block_path, "777")
  time_path = paste(f, '_bootstrap_resamples/',sep="")
  save_dir = paste(block_path, time_path, sep="")
  dir.create(save_dir)
  n_pc = 118
  n_splits = 100
  n_rows = nrow(all_data)
  n_cols = ncol(x_data)
  #n_tr = nrow(train_data)
  n_tr = round(0.8*n_rows)
  pos_class = positive_classes[n]
  
  ### Baseline Models
  # Remove Zero Variance
  data_var = apply(x_data, 2, sd)
  no_var_idx = which(data_var == 0.0)
  if (length(no_var_idx) > 0) {
    all_data = all_data[, -no_var_idx]
  }

  ### generate splits for the data, or load pre-existing file
  # splits = matrix(0,nrow=n_splits,ncol=n_tr)
  # for (n in 1:n_splits) {
  #   splits[n,] = sort(sample.int(n_rows,n_tr))
  # }
  # # Save the splits
  # cat("\n---------------\nSave Splits...\n---------------\n")
  # write.table(splits,paste(save_dir,f,"_bootstrap_resamples.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)
  splits = as.matrix(read.table(paste("~/Dropbox-Work/Wuxi/Data/",f,"_bootstrap_resamples.txt",sep="")))

  for (j in 1:n_splits) {
    cat(paste("\nSplit ", as.character(j) ,"\n", sep=''))
    # Split Labels
    class_train = y_data[splits[j,]]
    class_test = y_data[-splits[j,]]
    # Split Data
    data_train = x_data[splits[j,],]
    data_test = x_data[-splits[j,],]
    # Remove Zero Variance
    train_sd = apply(data_train, 2, sd)
    no_var_idx = which(train_sd == 0.0)
    if (length(no_var_idx) > 0) {
      data_train = data_train[, -no_var_idx]
      data_test = data_test[, -no_var_idx]
      train_sd = train_sd[-no_var_idx]
    }
    
    # Normalize All Data
    cat("\n---------------\nNormalize Data...\n---------------\n")

    train_mean = apply(data_train, 2, mean)
    #train_sd = apply(data_train, 2, sd)
    # Train
    x_train = sweep(data_train, 2, train_mean, "-")
    x_train = sweep(x_train, 2, train_sd, "/")
    # Test
    x_test = sweep(data_test, 2, train_mean, "-")
    x_test = sweep(x_test, 2, train_sd, "/")
    
    
    cat("\n---------------\nPrincipal Components...\n---------------\n")
    # PC Model
    pc_res = switch_pca(x_train, n_pc = n_pc)
    # Train
    pc_train = predict_correct_pca(pc_res, x_train, n_pc)
    # Test
    pc_test = predict_correct_pca(pc_res, x_test, n_pc)
    # Output for DW
    dw_train = cbind(class_train==pos_class,pc_train)
    dw_test = cbind(class_test==pos_class,pc_test)
    #dw_train = cbind(as.numeric(class_train)-1,pc_train)
    #dw_test = cbind(as.numeric(class_test)-1,pc_test)
    writeMat(paste(save_dir,"resample_", as.character(j), "_data.mat", sep=''),
             traindata=dw_train,testdata=dw_test)
  }
}
