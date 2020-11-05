# Script to run Nonparametric Wilcoxon signed-rank tests 
# Author: Sharvari Gujja

library("dplyr")
library("ggpubr")

# rm(list=ls())
# gc()
options(stringsAsFactors = FALSE)

# new data sets from Richard Feb 3rd 2020
#################################################

# 
# 1. lumAB (serial dilution) #17 fractions 
# 25% for D-Wave versus SVM ; n = 8 tests (9 methods) 
# D-Wave at 95% versus 20% ; n = 144 tests ((17-1) * 9) 
# Lasso at 95% versus 20% ; n = 144 tests ((17-1) * 9) 
#   
# 2. ERpn (serial dilution) # 19 fractions
# dw at 95% versus 10% ; n = 162 tests (18*9)
# Lasso at 95% versus 10%; n = 162 tests (18*9)
#   
# 3. 6cancer (serial dilution) # 19 fractions
# 95% to 5%  dw ; n = 162 tests (18*9)
# 95% to 5% lasso ; n = 162 tests (18*9)
# 
# 4. lumAB_bootstrap_resamples
# RF PC vs gene-level (from pc1)

##################################################

all_data <- readRDS("~/OneDrive - NextCODE Health/Projects/quantumML/Classification/output/lumAB_bacc_pop.RDS")

#25% for D-Wave versus SVM
all_data_svm_dw_25 <- all_data[all_data$method %in% c("svm", "D-Wave") & all_data$frac %in% c("0.25"),]
all_data_svm_dw_25<-  as.data.frame(all_data_svm_dw_25[,c("method","val_bacc")])
# Compute unpaired Wilcoxon-test 
res <- wilcox.test(val_bacc ~ method, data = all_data_svm_dw_25, paired = FALSE)
#res
pval <- res$p.value
padj <-p.adjust(res$p.value, method = "bonferroni", n = 36)
print(padj)

# D-Wave at 95% versus 20%
all_data_dw_95_20 <- all_data[all_data$method %in% c("D-Wave") & all_data$frac %in% c("0.95","0.2"),]
all_data_dw_95_20 <-  as.data.frame(all_data_dw_95_20[,c("frac","val_bacc")])
# Compute unpaired Wilcoxon-test
res <- wilcox.test(val_bacc ~ frac, data = all_data_dw_95_20, paired = FALSE)
#res
pval <- res$p.value
padj <-p.adjust(res$p.value, method = "bonferroni", n = 144)
print(padj)

# Lasso at 95% versus 20%
all_data_lasso_95_20 <- all_data[all_data$method %in% c("lasso") & all_data$frac %in% c("0.95","0.2"),]
all_data_lasso_95_20 <-  as.data.frame(all_data_lasso_95_20[,c("frac","val_bacc")])
# Compute unpaired Wilcoxon-test 
res <- wilcox.test(val_bacc ~ frac, data = all_data_lasso_95_20, paired = FALSE)
#res
pval <- res$p.value
padj <-p.adjust(res$p.value, method = "bonferroni", n = 144)
print(padj)


################################################################################

all_data <- readRDS("~/OneDrive - NextCODE Health/Projects/quantumML/Classification/output/ERpn_bacc_pop.RDS")

# D-Wave at 95% versus 10%
all_data_dw_95_10 <- all_data[all_data$method %in% c("D-Wave") & all_data$frac %in% c("0.95","0.1"),]
all_data_dw_95_10 <-  as.data.frame(all_data_dw_95_10[,c("frac","val_bacc")])
# Compute unpaired Wilcoxon-test 
res <- wilcox.test(val_bacc ~ frac, data = all_data_dw_95_10, paired = FALSE)
#res
pval <- res$p.value
padj <-p.adjust(res$p.value, method = "bonferroni", n = 162)
print(padj)

# Lasso at 95% versus 10%
all_data_lasso_95_10 <- all_data[all_data$method %in% c("lasso") & all_data$frac %in% c("0.95","0.1"),]
all_data_lasso_95_10 <-  as.data.frame(all_data_lasso_95_10[,c("frac","val_bacc")])
# Compute unpaired Wilcoxon-test
res <- wilcox.test(val_bacc ~ frac, data = all_data_lasso_95_10, paired = FALSE)
#res
pval <- res$p.value
padj <-p.adjust(res$p.value, method = "bonferroni", n = 162)
print(padj)

################################################################################

all_data <- readRDS("~/OneDrive - NextCODE Health/Projects/quantumML/Classification/output/6cancer_bacc_pop.RDS")

# D-Wave at 95% versus 5%
all_data_dw_95_5 <- all_data[all_data$method %in% c("D-Wave") & all_data$frac %in% c("0.95","0.05"),]
all_data_dw_95_5 <-  as.data.frame(all_data_dw_95_5[,c("frac","val_bacc")])
# Compute unpaired Wilcoxon-test 
res <- wilcox.test(val_bacc ~ frac, data = all_data_dw_95_5, paired = FALSE)
#res
pval <- res$p.value
padj <-p.adjust(res$p.value, method = "bonferroni", n = 162)
print(padj)

# Lasso at 95% versus 5%
all_data_lasso_95_5 <- all_data[all_data$method %in% c("lasso") & all_data$frac %in% c("0.95","0.05"),]
all_data_lasso_95_5 <-  as.data.frame(all_data_lasso_95_5[,c("frac","val_bacc")])
# Compute unpaired Wilcoxon-test 
res <- wilcox.test(val_bacc ~ frac, data = all_data_lasso_95_5, paired = FALSE)
#res
pval <- res$p.value
padj <-p.adjust(res$p.value, method = "bonferroni", n = 162)
print(padj)

################################################################################

all_data <- readRDS("~/OneDrive - NextCODE Health/Projects/quantumML/Classification/output/lumAB_pca_vs_gene_rf_bacc_pop.RDS")
# RF PC vs. gene-level

# Compute unpaired Wilcoxon-test 
res <- wilcox.test(tst_bacc ~ dataset, data = all_data, paired = FALSE)
#res
pval <- res$p.value
padj <-p.adjust(res$p.value, method = "bonferroni", n = 9)
print(padj)
