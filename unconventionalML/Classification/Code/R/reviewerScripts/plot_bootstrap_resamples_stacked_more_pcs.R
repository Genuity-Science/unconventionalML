# Script to plot results from bootstrap resamples for binomial datasets. 
# Compares classical, DW, SA, random, field results.


rm(list=ls())
gc()
options(stringsAsFactors = FALSE)

library(data.table)
library(ggplot2)
#library(gridExtra)
#require(RColorBrewer)
library(dplyr)
library(stringr)
#library(grid)
library(tidyr)

base_dir ='/Users/rli/OneDrive - NextCODE Health/Projects/quantumML/Classification/quantum/Results/bootstrap_resamples_more_pcs/'
sem <- function(x) {sd(x)/sqrt(length(x))}

# define datasets to loop over 
datasets = c("brcaMatchedTN","ERpn","kirckirp","luadlusc","lumAB")

# Classical info
all_cl_info = data.frame()
for (n in 1:length(datasets)) {
  dataset = datasets[[n]]
  l = readRDS(paste(base_dir, dataset,'_max_pcs_multirun_save_new.RDS',sep="")) #_multirun_save.RDS
  info = l$info
  info$dataset = dataset
  all_cl_info = rbind(all_cl_info,info)
}

all_info = all_cl_info
all_info[is.na(all_info)] = 0 # substitute 0 for NA 
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
algo_df <- algo_df %>%  mutate(method = str_replace(method, "dw", "D-Wave"))
algo_df <- algo_df %>%  mutate(method = str_replace(method, "sa", "SA"))

levels(algo_df$metric)=c("Accuracy","Bal.Accuracy","AUC","Precision","Recall","F1 score")
algo_df$method = factor(algo_df$method,levels=c("SVM","LASSO","Ridge","RF","NB"))
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

nleg = length(unique(algo_df2$dataset))
nmeth = length(unique(algo_df2$method))
nmetric = length(unique(algo_df2$metric))
final_methods = intersect(levels(algo_df2$method),unique(algo_df2$method))

tmp_df <- algo_df2 %>% group_by(method)
tmp_df = tmp_df[with(tmp_df,order(dataset,metric,method)),]
bacc_order = tmp_df %>% filter(metric=="Bal.Accuracy") %>% group_by(dataset) %>% mutate(sorto=(order(-value)),newo=order(sorto))
tmp_m = matrix(bacc_order$newo,nrow=nmeth)
group_order = as.vector(tmp_m[,rep(1:nleg,each=nmetric)])
tmp_df$bacc_order = group_order

tmp_df = tmp_df %>% ungroup() %>% group_by(dataset,metric) %>% arrange(dataset,bacc_order) %>%
  unite("dataset_baccorder",dataset,bacc_order,sep="_",remove=FALSE) %>% ungroup() %>% rename(Metric=metric)

tmp_df$dataset_baccorder = factor(tmp_df$dataset_baccorder)

p = ggplot(tmp_df,aes(x=dataset_baccorder,y=value,group=Metric)) + 
  geom_errorbar(aes(color=Metric,width=0.1,ymin=(value-sem),ymax=(value+sem))) + 
  geom_line(aes(color=Metric)) +
  facet_wrap(vars(dataset),scales="free",nrow=nleg) + 
  theme_bw() + 
  scale_x_discrete(breaks=tmp_df$dataset_baccorder,labels=tmp_df$method,expand=c(0.017,0)) +
  xlab("Algorithm") + ylab("Value") + theme(legend.position="bottom",text=element_text(size=15),axis.text.x=element_text(hjust=0.6),legend.background=element_rect(size=0.25,linetype="solid",color="black"))
#+ theme(legend.position=c(0.22,0.905),legend.background=element_rect(size=0.25,linetype="solid",color="black"),legend.title.align=0.5,text=element_text(size=15)) + guides(color=guide_legend(ncol=2))

