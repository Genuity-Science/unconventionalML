# Output .mat files with resample data with genes from PC1. 

library(feather)
library(R.matlab)
library(flashpcaR)

block_path = '~/Dropbox-Work/Wuxi/Data/'

f = 'lumAB'
cat(paste("Dataset: ", f ,"\n", sep=''))
cat("\n---------------\nLoading train data...\n---------------\n")
train_data = read_feather(paste(block_path,"data5_",f,"_train_normalized.feather",sep=""))

cat("\n---------------\nLoading test data...\n---------------\n")
test_data = read_feather(paste(block_path,"data5_",f,"_test_normalized.feather",sep=""))

all_data = rbind(train_data,test_data)

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
n_splits = 100
n_rows = nrow(all_data)
n_cols = ncol(x_data)
#n_tr = nrow(train_data)
n_tr = round(0.8*n_rows)
pos_class = "Luminal_A"

### generate splits for the data, or load pre-existing file
splits = as.matrix(read.table(paste("~/Dropbox-Work/Wuxi/Data/",f,"_bootstrap_resamples.txt",sep="")))
gene_ids = fread('lumAB_pc1_top_44_genes_ids.txt',header=F)
coln = colnames(x_data)
idxs = match(t(as.vector(ids)),coln)
for (j in 1:n_splits) {
  cat(paste("\nSplit ", as.character(j) ,"\n", sep=''))
  # Split Labels
  class_train = y_data[splits[j,]]
  class_test = y_data[-splits[j,]]
  # Split Data
  data_train = x_data[splits[j,],idxs]
  data_test = x_data[-splits[j,],idxs]
  
  # Output for DW
  dw_train = cbind(class_train==pos_class,data_train)
  dw_test = cbind(class_test==pos_class,data_test)
  #dw_train = cbind(as.numeric(class_train)-1,pc_train)
  #dw_test = cbind(as.numeric(class_test)-1,pc_test)
  writeMat(paste(save_dir,"gene_resample_", as.character(j), "_data.mat", sep=''),
           traindata=as.matrix(dw_train),testdata=as.matrix(dw_test))
}
