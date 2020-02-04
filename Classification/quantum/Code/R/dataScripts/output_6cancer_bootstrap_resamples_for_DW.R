# Output .mat file with multinomial cuts for DW. 
# Example usage: Rscript output_6cancer_cuts_for_DW.R j
# Did this because R has pretty poor memory management and would run out of 
# space. This way only handle one at a time and can do the pca

args = commandArgs(trailingOnly=TRUE)

library(feather)
library(R.matlab)
library(flashpcaR)


# identify dataset and load 
f = "6cancer"
all_data=read_feather("/boston_scratch/TCGA_Pan6Cancer/combinedData/data5_6cancer_all.feather")
block_path = '/boston_ailab/users/rli/quantum/Data/'
y_data = all_data$cancer

# fix labels so that don't have space in them
x_data = all_data[,-1]

# Run name
run_name = paste("_",f,sep="")
Sys.chmod(block_path, "777")
time_path = paste(f, '_bootstrap_resamples/',sep="")
save_dir = paste(block_path, time_path, sep="")

# dir.create(save_dir)
# n_pc = 44
n_pc = 13
n_splits = 100
n_rows = nrow(all_data)
n_cols = ncol(x_data)
#n_tr = nrow(train_data)
#n_tr = round(0.8*n_rows)

splits=as.matrix(read.table("/boston_ailab/users/rli/quantum/Data/6cancer_bootstrap_resamples.txt",header=T))
#   splits = createDataPartition(y_data, times = n_splits, p = 0.8)
#split_name = "/boston_ailab/users/rli/quantum/Data/6cancer_bootstrap_resamples.txt"
   # Save the splits
#   cat("\n---------------\nSave Splits...\n---------------\n")
#   write.table(split_name,x=as.data.frame(splits),quote=F,row.names=F,sep='\t')
j = as.numeric(args[1])
cat(paste("\nSplit ", as.character(j) ,"\n", sep=''))
# Split Labels
class_train = y_data[splits[,j]]
class_test = y_data[-splits[,j]]
# Split Data
data_train = x_data[splits[,j],]
data_test = x_data[-splits[,j],]
# Remove Zero Variance
train_var = apply(data_train, 2, sd)
no_var_idx = which(train_var == 0.0)
if (length(no_var_idx) > 0) {
  data_train = data_train[, -no_var_idx]
  data_test = data_test[, -no_var_idx]
}
gc()
# Normalize All Data
cat("\n---------------\nNormalize Data...\n---------------\n")
# train_mean = vector("numeric",n_cols)
# train_sd = vector("numeric",n_cols)

train_mean = apply(data_train, 2, mean)
gc()
train_sd = apply(data_train, 2, sd)
gc()
# Train
x_train = sweep(data_train, 2, train_mean, "-")
x_train = sweep(x_train, 2, train_sd, "/")
# Test
x_test = sweep(data_test, 2, train_mean, "-")
x_test = sweep(x_test, 2, train_sd, "/")
# 
# 
cat("\n---------------\nPrincipal Components...\n---------------\n")
# PC Model
pc_model = flashpca(as.matrix(data_train), ndim = n_pc, stand = "sd", do_loadings = TRUE)
# PCs of All Data
# Train
pc_train = pc_model$projection
# Test
pc_test = project(as.matrix(data_test), loadings = pc_model$loadings,  orig_mean = pc_model$center, orig_sd = pc_model$scale)
pc_test = pc_test$projection
# Output for DW
#dw_train = cbind(class_train==pos_class,pc_train)
#dw_test = cbind(class_test==pos_class,pc_test)
dw_train = cbind(as.numeric(class_train)-1,pc_train)
dw_test = cbind(as.numeric(class_test)-1,pc_test)
writeMat(paste(save_dir,"resample_", as.character(j), "_data.mat", sep=''),
         traindata=dw_train,testdata=dw_test)
gc()
