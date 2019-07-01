library("dplyr")
library("ggpubr")

rm(list=ls())
gc()
options(stringsAsFactors = FALSE)

# write.csv(all_info,"qMl_classic_dw_sa_ibm_rand_new.csv",quote=F,sep=",",row.names=F)
#qml <- read.csv("/Users/sgujja/OneDrive - NextCODE Health/Omar_Code/quantumMachineLearning - Documents/Cuts/qMl_classic_dw_sa_ibm_rand_new.csv", sep=",",header=T)
qml <- all_info # file from running plot_bootstrap_resamples_stacked.R

d <- data.frame()

for(i in c("brcaMatchedTN","ERpn","kirckirp","luadlusc","lumAB") )
{
  for (j in c("Field","lasso","nb","Random","rf","ridge","sa","svm"))
  {
    for (k in c("tst_acc","tst_bacc","tst_auroc","tst_F1"))
    {
      qml.sub <- qml[qml$dataset %in% c(i) & qml$method %in% c("dw", j),]
      qml.sub <-  as.data.frame(qml.sub[,c("method",k)])
      names(qml.sub) <- c("group", "metric")
      
      #Compute summary statistics by groups:
      # group_by(qml.sub, group) %>%
      #   summarise(
      #     count = n(),
      #     median = median(metric, na.rm = TRUE),
      #     IQR = IQR(metric, na.rm = TRUE)
      #   )
      
      # Plot weight by group and color by group
      # show(ggboxplot(qml.sub, x = "group", y = "metric", 
      #           color = "group", palette = c("#00AFBB", "#E7B800"),
      #           order = c("dw", i),
      #           ylab = "Metric", xlab = "Groups"))
      
      
      # Wilcoxin Rank Sum Test
      # Compute paired Wilcoxon-test 
      #res <- wilcox.test(metric ~ group, data = qml.sub, paired = TRUE)
      # Compute unpaired Wilcoxon-test 
      res <- wilcox.test(metric ~ group, data = qml.sub, paired = FALSE)
      #res
      pval <- res$p.value
      padj <-p.adjust(res$p.value, method = "bonferroni", n = 8)
      info <- c(i,k,paste0("dave.vs.",j),pval, padj)
      d <- rbind(d,info)
      names(d) <- c("dataset", "metric","comparison", "pval", "padj")
  }
  saveRDS(d,paste0(i,"_",k,"_wilcox_paired.rds"))
}
}
write.csv(d,"binomial_wilcoxon_unpaired.csv", quote = F, row.names = F)



###########################
# Wilcoxon Rank sum Test for reduction in training data sizes for lasso and DW at 20% and 95%
############################

reduc.train <- readRDS("training_size_summary.RDS")
reduc.train <- as.data.frame(reduc.train$dw_lasso_population_data)

print(table(reduc.train$frac))
print(table(reduc.train$method))

lasso.sub <- reduc.train[reduc.train$method %in% c("lasso"),]
print(table( lasso.sub$frac))
#0.2 0.95 
#50   50 
lasso.sub <-  as.data.frame(lasso.sub[,c("frac","val_bacc")])
names(lasso.sub) <- c("group", "metric")


dw.sub <- reduc.train[reduc.train$method %in% c("D-Wave"),]
print(table(dw.sub$frac))
#0.2 0.95 
#50   50 
dw.sub <-  as.data.frame(dw.sub[,c("frac","val_bacc")])
names(dw.sub) <- c("group", "metric")

# Compute unpaired Wilcoxon-test  for lasso at 0.95 and 0.2 fractions
res <- wilcox.test(metric ~ group, data = lasso.sub, paired = FALSE)
#res
pval <- res$p.value
#
padj <-p.adjust(res$p.value, method = "bonferroni", n =135)
print(padj)


# Compute unpaired Wilcoxon-test for D-Wave at 0.95 and 0.2 fractions
res <- wilcox.test(metric ~ group, data = dw.sub, paired = FALSE)
#res
pval <- res$p.value
padj <-p.adjust(res$p.value, method = "bonferroni", n =135)
print(padj)



##################################################
# calculate mean and sem values
###############################################
sem <- function(x) {sd(x)/sqrt(length(x))}

lasso.sub$method <- "lasso"
dw.sub$method <- "D-Wave"

algo_df <- rbind(lasso.sub, dw.sub)
mean_stats = aggregate(.~group+method,algo_df,mean)
sem_stats = aggregate(.~group+method,algo_df,sem)
n_cols = ncol(mean_stats)
melt_test_sem = melt(sem_stats,id=c(1,2),variable.name="val_bacc",value.name="sem",measure=3) 
melt_test_mean = melt(mean_stats,id=c(1,2),variable.name="val_bacc",value.name="mean",measure=3)
algo_df_stats = merge(melt_test_mean,melt_test_sem)


print(algo_df_stats)


#################################################
#unpaired wilcoxon rank sum test for D-Wave vs. SVM at 20%; adjusting p-values for multiple correction with n=8 tests:
##################################################

all_data <- readRDS("all_bacc_pop_data.RDS")

all_data_svm_dw_20 <- all_data[all_data$method %in% c("svm", "D-Wave") & all_data$frac %in% c("0.2"),]
all_data_svm_dw_20<-  as.data.frame(all_data_svm_dw_20[,c("method","val_bacc")])


all_data_svm_dw_30 <- all_data[all_data$method %in% c("svm", "D-Wave") & all_data$frac %in% c("0.3"),]
all_data_svm_dw_30<-  as.data.frame(all_data_svm_dw_30[,c("method","val_bacc")])

# this section of the code is run for each fraction from 0.2 to 0.5
all_data_svm_dw_frac <- all_data[all_data$method %in% c("svm", "D-Wave") & all_data$frac %in% c("0.5"),]
all_data_svm_dw_frac<-  as.data.frame(all_data_svm_dw_frac[,c("method","val_bacc")])

# Compute unpaired Wilcoxon-test for D-Wave vs. SVM at  0.2 fraction
res <- wilcox.test(val_bacc ~ method, data = all_data_svm_dw_frac, paired = FALSE)
#res
pval <- res$p.value
padj <-p.adjust(res$p.value, method = "bonferroni", n =8)
print(padj)

##################################################
# calculate mean and sem values for D-Wave and SVM for 0.2 to 0.5
###############################################
sem <- function(x) {sd(x)/sqrt(length(x))}

svm.sub.20.50 <- all_data[all_data$method %in% c("svm") & all_data$frac %in% c("0.2","0.25","0.3","0.35","0.4","0.45","0.5"),]
dw.sub.20.50 <- all_data[all_data$method %in% c("D-Wave") & all_data$frac %in% c("0.2","0.25","0.3","0.35","0.4","0.45","0.5"),]

algo_df <- rbind(svm.sub.20.50, dw.sub.20.50)
mean_stats = aggregate(.~frac+method,algo_df,mean)
sem_stats = aggregate(.~frac+method,algo_df,sem)
n_cols = ncol(mean_stats)
melt_test_sem = melt(sem_stats,id=c(1,2),variable.name="val_bacc",value.name="sem",measure=3) 
melt_test_mean = melt(mean_stats,id=c(1,2),variable.name="val_bacc",value.name="mean",measure=3)
algo_df_stats = merge(melt_test_mean,melt_test_sem)

print(algo_df_stats)

#################################################
#unpaired wilcoxon rank sum test for RF balanced accuracy at 44 genes and 44 pc level for lumAB ; adjusting p-values for multiple correction with n=8 tests:
##################################################

# read lumAB after running 100 runs on 44 genes derived from top pc1 for a single cut
l.gene = readRDS("lumAB_multicut_multirun/2019_04_09_LUMAB_ALL_MULTICUT/lumab_multicut_multirun_44genes_save.RDS")
l.gene.info <- l.gene$info
l.gene.info <- cbind(dataset = "lumAB_gene", l.gene.info)
l.gene.info <-  as.data.frame(l.gene.info[l.gene.info$method %in% c("rf"),c("dataset","tst_bacc")])
print(dim(l.gene.info))

# read lumAB  after running 100 runs on 44 pcs
l.pc = readRDS('lumAB_multirun_save.RDS')
l.pc.info <- l.pc$info
l.pc.info <- cbind(dataset = "lumAB_pc", l.pc.info)
l.pc.info <-  as.data.frame(l.pc.info[l.pc.info$method %in% c("rf"),c("dataset","tst_bacc")])
print(dim(l.pc.info))

l.gene.pc <- as.data.frame(rbind(l.gene.info,l.pc.info))
names(l.gene.pc) <- c("group", "metric")

# Compute unpaired Wilcoxon-test  for lasso at 0.95 and 0.2 fractions
res <- wilcox.test(metric ~ group, data = l.gene.pc, paired = FALSE)
#res
pval <- res$p.value
#
padj <-p.adjust(res$p.value, method = "bonferroni", n = 9)
print(padj)
