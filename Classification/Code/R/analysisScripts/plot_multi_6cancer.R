# Script to plot results for 6-cancer multinomial data. 
# from PC1. Compares classical, DW, SA, random, field results.


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
# library(ggpubr)
library(tidyr)

# method_levels = rev(c("SVM","SA","Ridge","RF","Random","NB","LASSO","Field","D-Wave"))

base_dir = '~/OneDrive - Genuity Science/Projects/quantumML/Classification/output/'
sem <- function(x) {sd(x)/sqrt(length(x))}

# simulated annealing info
all_sa_info = readRDS(paste(base_dir,"6cancer_bootstrap_resamples_sa_nsols20_ntotsols1000.RDS",sep=""))

# D-Wave info
all_dw_info = readRDS(paste(base_dir,"6cancer_bootstrap_resamples_dw_nsols20_ntotsols1000.RDS",sep=""))

#Random info with 20 solutions
all_rand_info = readRDS(paste(base_dir,"6cancer_bootstrap_resamples_rand_nsols20_ntotsols1000.RDS",sep=""))

all_field_info = readRDS(paste(base_dir,"6cancer_bootstrap_resamples_field.RDS",sep=""))

all_rbm_info = readRDS(paste(base_dir,"6cancer_rbm_save_all.RDS",sep=''))$info

# classical info
all_cl_info = readRDS(paste(base_dir,"6cancer_multirun_save.RDS",sep=""))
all_cl_info.1 <- all_cl_info$info

# Combining DW, Classical, SA, IBM sim results in one data frame
all_info = rbind(all_dw_info,all_cl_info.1,all_sa_info,all_rand_info,all_field_info,all_rbm_info)

mean_stats = aggregate(.~method+dataset,all_info,mean)
sem_stats = aggregate(.~method+dataset,all_info,sem)
sd_stats = aggregate(.~method+dataset,all_info,sd)
n_cols = ncol(mean_stats)

melt_test_sem = melt(sem_stats,id=c(1,2),variable.name="metric",value.name="sem",measure=seq(4,n_cols,2)) 
melt_test_sd = melt(sd_stats,id=c(1,2),variable.name="metric",value.name="sd",measure=seq(4,n_cols,2)) 
melt_test_mean = melt(mean_stats,id=c(1,2),variable.name="metric",value.name="value",measure=seq(4,n_cols,2))

algo_df = merge(melt_test_mean,melt_test_sd)
algo_df = merge(algo_df,melt_test_sem)

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
algo_df <- algo_df %>%  mutate(method = str_replace(method, "rand", "Random"))
algo_df <- algo_df %>%  mutate(method = str_replace(method, "field", "Field"))

levels(algo_df$metric)=c("Accuracy","Bal.Accuracy","AUC","Precision","Recall","F1")
method_levels = unique(algo_df$method)
algo_df$method = factor(algo_df$method,levels=method_levels)
nleg = length(levels(algo_df$method))
bacc_algo_df = subset(algo_df,algo_df$metric %in% c("Bal.Accuracy")) 

algo_df_sub = subset(algo_df,algo_df$dataset %in% c("6cancer"))
algo_df_sub = subset(algo_df_sub,algo_df_sub$metric %in% c("Bal.Accuracy", "Accuracy","AUC", "F1 score"))
algo_df_sub$method = factor(algo_df_sub$method,levels=method_levels)
algo_df_sub$metric <- factor(algo_df_sub$metric)
algo_df_sub$metric2 = factor(algo_df_sub$metric, levels = levels(algo_df_sub$metric)[c(4,3,2,1 )]) ## This did the trick of changing the order

legend_ord <- method_levels
#
algo_df_sub1 = algo_df_sub %>% 
  group_by(metric2) %>% 
  mutate(position = rank(value))

summary_colors = brewer.pal(9,"Paired")

p2 <- ggplot(algo_df_sub1, aes(x=metric2, y=value, fill=reorder(method, value), group=position)) + 
  geom_bar(position=position_dodge(), stat="identity", width=0.6) + #, colour="black") +# adds borders
  geom_errorbar(aes(ymin=value, ymax=value+sd),
                width=.6, size=.3  ,                 # Width of the error bars
                position=position_dodge(.6), color="black") + 
  scale_fill_manual(values=summary_colors,
                      limits = legend_ord, 
                      labels = legend_ord)  + 
  coord_flip()  +  theme_bw() +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = -0.1)) +
  theme(axis.text = element_text(size = 12), text = element_text(size=12),
        legend.title=element_blank(), axis.title.x=element_blank(),
        axis.title.y=element_blank(),panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black")) # + xlab("Classifier")

show(p2)

algo_df2 = subset(algo_df_sub,method!="IBMsim")
nleg = length(unique(algo_df2$dataset))
nmeth = length(unique(algo_df2$method))
nmetric = length(unique(algo_df2$metric))
final_methods = intersect(levels(algo_df2$method),unique(algo_df2$method))

tmp_df <- algo_df2 %>% group_by(method)
tmp_df = tmp_df[with(tmp_df,order(dataset,metric,method)),]
bacc_order = tmp_df %>% filter(metric=="Accuracy") %>% group_by(dataset) %>% mutate(sorto=(order(-value)),newo=order(sorto))
tmp_m = matrix(bacc_order$newo,nrow=nmeth)
group_order = as.vector(tmp_m[,rep(1:nleg,each=nmetric)])
tmp_df$bacc_order = group_order

tmp_df = tmp_df %>% ungroup() %>% group_by(dataset,metric) %>% arrange(dataset,bacc_order) %>%
  unite("dataset_baccorder",dataset,bacc_order,sep="_",remove=FALSE) %>% ungroup()

tmp_df$dataset_baccorder = factor(tmp_df$dataset_baccorder,levels = unique(tmp_df$dataset_baccorder))

p_test = ggplot(tmp_df,aes(x=dataset_baccorder,y=value,group=metric)) + 
  geom_errorbar(aes(color=metric,width=0.1,ymin=(value-sem),ymax=(value+sem))) + 
  geom_line(aes(color=metric)) + 
  facet_wrap(vars(dataset),scales="free",nrow=nleg) + 
  theme_bw() + 
  scale_x_discrete(breaks=tmp_df$dataset_baccorder,labels=tmp_df$method) +
  xlab("Algorithm") + ylab("Value") + labs(color="Metric")  +
  theme(legend.position="bottom",text=element_text(size=16)) 
