
library(feather)
library(R.matlab)

files=c("brcaMatchedTN","ERpn","kirckirp","luadlusc","lumAB")#,"6cancer") "lumAB"
base_dir = '~/Dropbox-Work/Wuxi/Data/'
save_dir = '~/Dropbox-Work/Wuxi/Data/'
positive_classes = c("tumor","Positive","kirc","luad","Luminal_A")

for (n in 1:length(files)) {
  f = files[n]
  cat(paste("Dataset: ", f ,"\n", sep=''))
  cat("\n---------------\nLoading train data...\n---------------\n")
  train_data = read_feather(paste(base_dir,"data5_",f,"_top44genes_train_normalized.feather",sep=""))

  cat("\n---------------\nLoading test data...\n---------------\n")
  test_data = read_feather(paste(base_dir,"data5_",f,"_top44genes_test_normalized.feather",sep=""))

#  all_data = rbind(train_data,test_data)
  block_path = '~/Dropbox-Work/Wuxi/Data/'

  # fix factor names
  train_labels = train_data[,1][[1]]
  train_labels = as.factor(gsub(" ", "_", train_labels, fixed = TRUE))
  test_labels = test_data[,1][[1]]
  test_labels = as.factor(gsub(" ", "_", test_labels, fixed = TRUE))
  train_x = train_data[,-1]
  test_x = test_data[,-1]
  pos_class = positive_classes[n]

  dw_train = as.matrix(cbind(train_labels==pos_class,train_x))
  dw_test = as.matrix(cbind(test_labels==pos_class,test_x))

  writeMat(paste(save_dir,f, "_top44genes_data.mat", sep=''),
    traindata=dw_train,testdata=dw_test)
}
