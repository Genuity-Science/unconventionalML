all_bacc_pop = readRDS('~/OneDrive - Genuity Science/Projects/quantumML/Classification/output/lumAB_sd_pca_bacc_pop.RDS')

methods = unique(all_bacc_pop$method)
pair_methods = combn(methods,2)

data_025_bacc = all_bacc_pop %>% subset(frac==0.25)

p_values = data.frame()
for (n in 1 : dim(pair_methods)[2]) {
  p_values[n,1:2] = pair_methods[,n]
  tmp_df = data_025_bacc %>% filter(method %in% pair_methods[,n])
  res <- wilcox.test(val_bacc ~ method, data = tmp_df, paired = FALSE)
  p_values[n,3] = res$p.value
}
colnames(p_values) = c("Alg1","Alg2", "orig")
p_values$bonferroni = p.adjust(p_values$orig,method='bonferroni')
p_values$BH = p.adjust(p_values$orig,method='BH')
p_values$Frac = 0.25

#svm_dw_df = p_values %>% subset(Alg1 %in% c("svm","D-Wave") | Alg2 %in% c("svm","D-Wave"))
svm_df = p_values %>% subset(Alg1 == "svm" | Alg2 == "svm")
p.adjust(svm_df$orig,method='BH')
dw_df = p_values %>% subset(Alg1=='D-Wave' | Alg2 == 'D-Wave')
p.adjust(dw_df$orig,method='BH')

##############################################

# do for all fractions, for all methods, relative to....95%? 
frac_95_bacc = all_bacc_pop %>% subset(frac==0.95)
all_fracs = sort(unique(all_bacc_pop$frac))
all_fracs = all_fracs[2:16]
all_frac_p_values = data.frame()
count = 1
for (n in 1:length(all_fracs)) {
  frac = all_fracs[n]
  for (m in 1:length(methods)) {
    tmp_df = all_bacc_pop[all_bacc_pop$frac %in% c(frac,0.95) & all_bacc_pop$method==methods[m],]
    res <- wilcox.test(val_bacc ~ frac, data = tmp_df, paired = FALSE)
    all_frac_p_values[count,1] = frac
    all_frac_p_values[count,2] = methods[m]
    all_frac_p_values[count,3] = res$p.value
    count = count + 1 
  }
}
colnames(all_frac_p_values) = c("Frac","method","orig")
all_frac_p_values$bonferroni = p.adjust(all_frac_p_values$orig,method='bonferroni')
all_frac_p_values$BH = p.adjust(all_frac_p_values$orig,method='BH')


frac_95_20_p_values = all_frac_p_values[all_frac_p_values$Frac==0.2,1:3]
frac_95_20_p_values$BH = p.adjust(frac_95_20_p_values$orig,method='BH')


#####################

