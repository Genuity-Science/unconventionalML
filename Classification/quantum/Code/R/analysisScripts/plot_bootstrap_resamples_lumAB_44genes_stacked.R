# Script to plot results from bootstrap resamples with the lumAB, using 44 genes
# from PC1. Compares classical, DW, SA, random, field results.


#rm(list=ls())
#gc()
options(stringsAsFactors = FALSE)

library(data.table)
library(ggplot2)
library(gridExtra)
require(RColorBrewer)
library(dplyr)
library(stringr)
library(grid)
library(tidyr)


base_dir ='/Users/rli/OneDrive - NextCODE Health/Projects/quantumML/Classification/quantum/Results/bootstrap_resamples/'
sem <- function(x) {sd(x)/sqrt(length(x))}

# read lumAB after running 100 runs on 44 genes derived from top pc1 for a single cut
l = readRDS(paste(base_dir,"lumAB_gene_bootstrap_resamples.RDS",sep=""))
l.info <- l$info
l.info <- cbind(dataset = "lumAB", l.info)

datasets = c("lumAB")

# simulated annealing info
all_sa_info = readRDS(paste(base_dir,"lumAB_gene_bootstrap_resamples_sa_nsols20.RDS",sep=""))
all_sa_info[is.na(all_sa_info)] <- 0

# DWave info
all_dw_info = readRDS(paste(base_dir,"lumAB_gene_bootstrap_resamples_dw_nsols20.RDS",sep=""))
all_dw_info[is.na(all_dw_info)] <- 0

# Field info
all_fl_info = readRDS(paste(base_dir,"lumAB_gene_bootstrap_resamples_field.RDS",sep=""))
all_fl_info[is.na(all_fl_info)] <- 0
all_fl_info$method <- "Field"
all_fl_info_lumAB <- filter(all_fl_info, dataset == "lumAB")

#Random info
all_rand_info = readRDS(paste(base_dir,"lumAB_gene_bootstrap_resamples_rand_nsols20.RDS",sep=""))
all_rand_info$method <- "Random"
all_rand_info[is.na(all_rand_info)] <- 0


# Combining DW, Classical, SA, random, field results in one data frame
all_info = rbind(all_dw_info,l.info,all_sa_info,all_rand_info,all_fl_info)

# make dataset name consistent
all_info$dataset <- "lumAB_gene"

# calculate mean and average
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
algo_df$method = factor(algo_df$method,levels=c("SVM","LASSO","Ridge","RF","NB" ,"SA","DWave","IBMsim","Random","Field")) #
algo_df2 = subset(algo_df,metric %in% c("Bal.Accuracy","Accuracy","AUC","F1 score"))

algo_df2 <- algo_df2 %>% 
  mutate(dataset = str_replace(dataset, "lumAB_gene", "LumA.vs.LumB"))

nleg = length(unique(algo_df2$dataset))
nmeth = length(unique(algo_df2$method))
nmetric = length(unique(algo_df2$metric))
final_methods = intersect(levels(algo_df2$method),unique(algo_df2$method))

# plot below is if want to keep the order of the methods fixed
# p = ggplot(algo_df2,aes(x=method,y=value,group=metric)) + geom_errorbar(aes(color=metric,width=0.1,ymin=(value-sem),ymax=(value+sem))) + geom_line(aes(color=metric))
# p = p + facet_wrap(vars(dataset),scales="free",nrow=5) + theme_bw()

tmp_df <- algo_df2 %>% group_by(method)
tmp_df = tmp_df[with(tmp_df,order(dataset,metric,method)),]
bacc_order = tmp_df %>% filter(metric=="Bal.Accuracy") %>% group_by(dataset) %>% mutate(sorto=(order(-value)),newo=order(sorto))
tmp_m = matrix(bacc_order$newo,nrow=nmeth)
group_order = as.vector(tmp_m[,rep(1:nleg,each=nmetric)])
tmp_df$bacc_order = group_order

tmp_df = tmp_df %>% ungroup() %>% group_by(dataset,metric) %>% arrange(dataset,bacc_order) %>%
  unite("dataset_baccorder",dataset,bacc_order,sep="_",remove=FALSE) %>% ungroup()

tmp_df$dataset_baccorder = factor(tmp_df$dataset_baccorder)

p_test = ggplot(tmp_df,aes(x=dataset_baccorder,y=value,group=metric)) + 
  geom_errorbar(aes(color=metric,width=0.1,ymin=(value-sem),ymax=(value+sem))) + 
  geom_line(aes(color=metric)) + 
  facet_wrap(vars(dataset),scales="free",nrow=nleg) + 
  theme_bw() + 
  scale_x_discrete(breaks=tmp_df$dataset_baccorder,labels=tmp_df$method,expand=c(0.03,0)) +
  theme(legend.position="bottom",text=element_text(size=16)) + 
  xlab("Algorithm") + ylab("Value") + labs(color="Metric")
