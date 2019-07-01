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
library(ggpubr)


base_dir = '~/OneDrive - NextCODE Health/Omar_Code/quantumMachineLearning - Documents/Cuts/'
sem <- function(x) {sd(x)/sqrt(length(x))}


# simulated annealing info
all_sa_info = readRDS(paste(base_dir,"6cancer_bootstrap_resamples_sa_nsols20.RDS",sep=""))

# DWave info
all_dw_info = readRDS(paste(base_dir,"6cancer_bootstrap_resamples_dw_nsols20.RDS",sep=""))

#Random info with 20 solutions
all_rand_info = readRDS(paste(base_dir,"6cancer_bootstrap_resamples_rand_nsols20.RDS",sep=""))

# classical info
all_cl_info = readRDS(paste(base_dir,"6cancer_bootstrap_resamples_save.RDS",sep=""))
all_cl_info.1 <- all_cl_info$info

# Combining DW, Classical, SA, IBM sim results in one data frame
all_info = rbind(all_dw_info,all_cl_info.1,all_sa_info,all_rand_info)

mean_stats = aggregate(.~method+dataset,all_info,mean)
#sem_stats = aggregate(.~method+dataset,all_info,sem)
sd_stats = aggregate(.~method+dataset,all_info,sd)
n_cols = ncol(mean_stats)

#melt_test_sem = melt(sem_stats,id=c(1,2),variable.name="metric",value.name="sem",measure=seq(4,n_cols,2)) 
melt_test_sd = melt(sd_stats,id=c(1,2),variable.name="metric",value.name="sd",measure=seq(4,n_cols,2)) 
melt_test_mean = melt(mean_stats,id=c(1,2),variable.name="metric",value.name="value",measure=seq(4,n_cols,2))

algo_df = merge(melt_test_mean,melt_test_sd)

#rename metrics
algo_df <- algo_df %>%  mutate(metric = str_replace(metric, "tst_acc", "Accuracy"))
algo_df <- algo_df %>%  mutate(metric = str_replace(metric, "tst_auroc", "AUC"))
algo_df <- algo_df %>% mutate(metric = str_replace(metric, "tst_bacc", "Bal.Accuracy"))
algo_df <- algo_df %>%  mutate(metric = str_replace(metric, "tst_F1", "F1"))
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
algo_df <- algo_df %>%  mutate(method = str_replace(method, "rand", "Random"))

levels(algo_df$metric)=c("Accuracy","Bal.Accuracy","AUC","Precision","Recall","F1")
#algo_df$method = factor(algo_df$method,levels=c("dw","lasso","ridge","svm","rf","nb","sa"))
algo_df$method = factor(algo_df$method,levels=c("SVM","SA","Ridge","RF","Random","NB","LASSO","DWave"))
#algo_df$dataset <- fct_rev(algo_df$dataset) %>% levels
nleg = length(levels(algo_df$method))
bacc_algo_df = subset(algo_df,algo_df$metric %in% c("Bal.Accuracy")) 

# ADD IBM BALANCED ACCURACY METRICS
#bacc_ibm <- read.csv("~/OneDrive - NextCODE Health/Omar_Code/Code/QISKit/Qclassify/sgujja/ibm_bacc_stats_5_binomials_df.csv", header=T)
#bacc_algo_df <- as.data.frame(rbind(bacc_algo_df,bacc_ibm))

#p_all = ggplot(data=bacc_algo_df,aes(x=dataset,y=value)) + geom_errorbar(aes(color=method,width=0.1,ymin=(value-sem),ymax=(value+sem))) + 
#  geom_point(aes(color=method,shape=method),size=2)+scale_shape_manual(values=0:(nleg-1)) 

#show(p_all)


# prepping for inset b
# update data frame for lumAB dataset

algo_df <- algo_df %>% 
  mutate(metric = str_replace(metric, "F1", "F1 score"))

algo_df_sub = subset(algo_df,algo_df$dataset %in% c("6cancer"))
algo_df_sub = subset(algo_df_sub,algo_df_sub$metric %in% c("Bal.Accuracy", "Accuracy","AUC", "F1 score"))
#algo_df_sub = subset(algo_df_sub,algo_df_sub$metric %in% c("Accuracy","AUC", "F1 score"))
#algo_df_sub$method <- factor(algo_df_sub$method, level=c("lasso","ridge","svm","rf","nb"))
algo_df_sub$method = factor(algo_df_sub$method,levels=c("SVM","SA","Ridge","RF","Random","NB","LASSO","DWave"))
algo_df_sub$metric <- factor(algo_df_sub$metric)
algo_df_sub$metric2 = factor(algo_df_sub$metric, levels = levels(algo_df_sub$metric)[c(4,3,2,1 )]) ## This did the trick of changing the order

#my_PuBu = brewer.pal(n = 9, "BrBG")[9:2] # remove lighter two hues and reverse the ordering
#my_PuBu = brewer.pal(n = 9, "YlGnBu")[6:2]

legend_ord <- c("SVM","SA","Ridge","RF","Random","NB","LASSO","DWave")
#
algo_df_sub1 = algo_df_sub %>% 
  group_by(metric2) %>% 
  mutate(position = rank(value))


p2 <- ggplot(algo_df_sub1, aes(x=metric2, y=value, fill=reorder(method, value), group=position)) + 
  geom_bar(position=position_dodge(), stat="identity", width=0.6) + #, colour="black") +# adds borders
  geom_errorbar(aes(ymin=value, ymax=value+sd),
                width=.6, size=.3  ,                 # Width of the error bars
                position=position_dodge(.6), color="black") + 
  scale_fill_manual(values = c("DWave" = '#01665E',"LASSO" = '#80CDC1',"NB" = '#C7EAE5',"SA" = '#F5F5F5',"Random" = '#F6E8C3',"Ridge" = '#DFC27D',"RF" = '#BF812D',"SVM" = '#8C510A'), 
                      limits = legend_ord, 
                      guide = guide_legend(reverse=TRUE),
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


#ggtitle("b") +
#grid.arrange(p1, p2, nrow = 1)

show(p2)


#############################################################################################
# generate p-values from unpaired Wilcoxon test and padj from Bonferronni correction with n=8
#############################################################################################

#all_info
d <- data.frame()
for(i in c("6cancer") )
{
  for (j in c("lasso","nb","rf","rand","ridge","sa","svm"))
  {
    for (k in c("tst_acc","tst_bacc","tst_auroc","tst_F1"))
    {
      all_info.sub <- all_info[all_info$dataset %in% c(i) & all_info$method %in% c("dw", j),]
      all_info.sub <-  as.data.frame(all_info.sub[,c("method",k)])
      names(all_info.sub) <- c("group", "metric")
      
      # Wilcoxin Rank Sum Test
      # Compute paired Wilcoxon-test 
      #res <- wilcox.test(metric ~ group, data = qml.sub, paired = TRUE)
      # Compute un-paired Wilcoxon-test 
      res <- wilcox.test(metric ~ group, data = all_info.sub, paired = FALSE)
      #res
      pval <- res$p.value
      padj <-p.adjust(res$p.value, method = "bonferroni", n = 7)
      info <- c(i,k,paste0("dave.vs.",j),pval, padj)
      d <- rbind(d,info)
      names(d) <- c("dataset", "metric","comparison", "pval", "padj")
    }
    saveRDS(d,paste0(i,"_",k,"_wilcox_paired.rds"))
  }
}

write.csv(d,"multinomial_wilcoxon_unpaired.csv", quote = F, row.names = F)
