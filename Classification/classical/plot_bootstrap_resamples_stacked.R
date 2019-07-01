rm(list=ls())
gc()
options(stringsAsFactors = FALSE)

library(data.table)
library(ggplot2)
library(gridExtra)
require(RColorBrewer)
library(dplyr)
library(stringr)
library(grid)
library(tidyr)



base_dir = '~/OneDrive - NextCODE Health/Omar_Code/quantumMachineLearning - Documents/Cuts/'
sem <- function(x) {sd(x)/sqrt(length(x))}

datasets = c("brcaMatchedTN","ERpn","kirckirp","luadlusc","lumAB")
# Classical info
all_cl_info = data.frame()
for (n in 1:length(datasets)) {
  dataset = datasets[[n]]
  #l = readRDS(paste(base_dir,dataset,"_bootstrap_resamples/", dataset,'_bootstrap_resamples_save.RDS',sep=""))
  l = readRDS(paste(base_dir, dataset,'_multirun_save.RDS',sep=""))
  #l = readRDS(paste(base_dir,"lumAB_singlecut_multirun/", dataset,'_singlecut_multirun_tuneNB_save.RDS',sep=""))
  #l = readRDS(paste(base_dir,"2018_12_14_LUMAB_ALL_MULTICUT/", dataset,'_multicut_multirun_tuneNB_save.RDS',sep=""))
  info = l$info 
  info$dataset = dataset
  all_cl_info = rbind(all_cl_info,info)
}

# read lumAB with auprc metrics
#l = readRDS("lumAB_bootstrap_resamples_save_sg.RDS")
#l.info <- l$info
#l.info <- cbind(dataset = "lumAB", l.info)

# simulated annealing info
all_sa_info = readRDS(paste(base_dir,"latest.files.from.richard/bootstrap_resamples_sa_nsols_20_ntotsols_1000[2].RDS",sep=""))
all_sa_info[is.na(all_sa_info)] <- 0

# DWave info
all_dw_info = readRDS(paste(base_dir,"latest.files.from.richard/bootstrap_resamples_dw_nsols_20_ntotsols_1000.RDS",sep=""))
all_dw_info[is.na(all_dw_info)] <- 0

# Field info
all_field_info = readRDS(paste(base_dir,"latest.files.from.richard/bootstrap_resamples_field.RDS",sep=""))
all_field_info[is.na(all_field_info)] <- 0
all_field_info$method <- "Field"

#Random info
all_rand_info = readRDS(paste(base_dir,"latest.files.from.richard/bootstrap_resamples_rand_nsols20.RDS",sep=""))
all_rand_info$method <- "Random"
all_rand_info[is.na(all_rand_info)] <- 0
all_rand_info <- all_rand_info[,!(names(all_rand_info) %in% c("tst_auprc","tr_auprc"))]


# add empty columns of AUPRC to classical and IBM sim methods
#colvector <- c("tr_auprc", "tst_auprc")
#all_cl_info[ , colvector] <- 0

#all_cl_info.lumAB[, colvector] <- 0

# Combining DW, Classical, SA, Random , Field results in one data frame # IBM sim
#all_info = rbind(all_dw_info,all_cl_info,l.info,all_sa_info,all_ibm_info,all_rand_info,all_fl_info)
all_info = rbind(all_cl_info,all_dw_info,all_sa_info,all_rand_info,all_field_info)
all_info[is.na(all_info)] = 0 # Take care of some datasets having NA for train and test auprc
mean_stats = aggregate(.~method+dataset,all_info,mean)
sem_stats = aggregate(.~method+dataset,all_info,sem)
n_cols = ncol(mean_stats)

melt_test_sem = melt(sem_stats,id=c(1,2),variable.name="metric",value.name="sem",measure=seq(4,n_cols,2)) 
melt_test_mean = melt(mean_stats,id=c(1,2),variable.name="metric",value.name="value",measure=seq(4,n_cols,2))

algo_df = merge(melt_test_mean,melt_test_sem)

#rename metrics
algo_df <- algo_df %>%  mutate(metric = str_replace(metric, "tst_acc", "Accuracy"))
algo_df <- algo_df %>%  mutate(metric = str_replace(metric, "tst_auroc", "AUC"))
algo_df <- algo_df %>% mutate(metric = str_replace(metric, "tst_bacc", "Bal.Accuracy"))
algo_df <- algo_df %>%  mutate(metric = str_replace(metric, "tst_F1", "F1 score"))
algo_df <- algo_df %>%  mutate(metric = str_replace(metric, "tst_prec", "Precision"))
algo_df <- algo_df %>%  mutate(metric = str_replace(metric, "tst_recall", "Recall"))

# rename methods
algo_df <- algo_df %>%  mutate(method = str_replace(method, "lasso", "LASSO"))
algo_df <- algo_df %>%  mutate(method = str_replace(method, "nb", "NB"))
algo_df <- algo_df %>% mutate(method = str_replace(method, "rf", "RF"))
algo_df <- algo_df %>%  mutate(method = str_replace(method, "ridge", "Ridge"))
algo_df <- algo_df %>%  mutate(method = str_replace(method, "svm", "SVM"))
algo_df <- algo_df %>%  mutate(method = str_replace(method, "dw", "DWave"))
algo_df <- algo_df %>%  mutate(method = str_replace(method, "sa", "SA"))

levels(algo_df$metric)=c("Accuracy","Bal.Accuracy","AUC","Precision","Recall","F1 score")
#algo_df$method = factor(algo_df$method,levels=c("dw","lasso","ridge","svm","rf","nb","sa"))
algo_df$method = factor(algo_df$method,levels=c("SVM","LASSO","Ridge","RF","NB","SA","DWave","Field","Random")) #"IBMsim"
#algo_df$method = factor(algo_df$method,levels=unique(algo_df$method)) # can specify a preferred order 
algo_df2 = subset(algo_df,metric %in% c("Bal.Accuracy","Accuracy","AUC","F1 score"))

algo_df2 <- algo_df2 %>% 
  mutate(dataset = str_replace(dataset, "lumAB", "LumA vs. LumB"))
algo_df2 <- algo_df2 %>% 
  mutate(dataset = str_replace(dataset, "brcaMatchedTN", "BRCA vs. Normal"))
algo_df2 <- algo_df2 %>% 
  mutate(dataset = str_replace(dataset, "ERpn", "ERpos vs. ERneg"))
algo_df2 <- algo_df2 %>% 
  mutate(dataset = str_replace(dataset, "kirckirp", "KIRC vs. KIRP"))
algo_df2 <- algo_df2 %>% 
  mutate(dataset = str_replace(dataset, "luadlusc", "LUAD vs. LUSC"))

#algo_df$dataset <- fct_rev(algo_df$dataset) %>% levels
nleg = length(unique(algo_df2$dataset))
nmeth = length(unique(algo_df2$method))
nmetric = length(unique(algo_df2$metric))
final_methods = intersect(levels(algo_df2$method),unique(algo_df2$method))
# ADD IBM BALANCED ACCURACY METRICS
#bacc_ibm <- read.csv("~/OneDrive - NextCODE Health/Omar_Code/Code/QISKit/Qclassify/sgujja/ibm_bacc_stats_5_binomials_df.csv", header=T)
#bacc_algo_df <- as.data.frame(rbind(bacc_algo_df,bacc_ibm))

# plot below is if want to keep the order of the methods fixed
# p = ggplot(algo_df2,aes(x=method,y=value,group=metric)) + geom_errorbar(aes(color=metric,width=0.1,ymin=(value-sem),ymax=(value+sem))) + geom_line(aes(color=metric))
# p = p + facet_wrap(vars(dataset),scales="free",nrow=5) + theme_bw()

tmp_df <- algo_df2 %>% group_by(method)
tmp_df = tmp_df[with(tmp_df,order(dataset,metric,method)),]
bacc_order = tmp_df %>% filter(metric=="Bal.Accuracy") %>% group_by(dataset) %>% mutate(sorto=(order(-value)),newo=order(sorto))
tmp_m = matrix(bacc_order$newo,nrow=nmeth)
group_order = as.vector(tmp_m[,rep(1:nleg,each=nmetric)])
tmp_df$bacc_order = group_order

#ungroup tmp_df??
tmp_df = tmp_df %>% ungroup() %>% group_by(dataset,metric) %>% arrange(dataset,bacc_order) %>%
  unite("dataset_baccorder",dataset,bacc_order,sep="_",remove=FALSE) %>% ungroup()

tmp_df$dataset_baccorder = factor(tmp_df$dataset_baccorder)

p_test = ggplot(tmp_df,aes(x=dataset_baccorder,y=value,group=metric)) + 
  geom_errorbar(aes(color=metric,width=0.1,ymin=(value-sem),ymax=(value+sem))) + 
  geom_line(aes(color=metric)) + 
  facet_wrap(vars(dataset),scales="free",nrow=nleg) + 
  theme_bw() + 
  scale_x_discrete(breaks=tmp_df$dataset_baccorder,labels=tmp_df$method) +
  xlab("Algorithm") + ylab("Value") 

show(p_test)
