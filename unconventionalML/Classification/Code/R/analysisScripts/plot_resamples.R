sem <- function(x) {
  sd(x,na.rm=T)/sqrt(length(x))
}

save_names = c("RIN_gte5.5_adjCounts_noRace_nodups_noXY_all_resamples_save.RDS",
               "RIN_gte5.5_adjCounts_noRace_nodups_noXY_DE_all_resamples_save.RDS",
               "RIN_gte5.5_adjCounts_noRace_nodups_noXY_mostVarGenes10pct_all_resamples_save.RDS",
               "all_samples_megena_adjCounts_noRace_nodups_noXY_all_resamples_save_new.RDS",
               "all_samples_nGOseq_adjCounts_noRace_nodups_noXY_all_resamples_save.RDS")

print_names = c("RIN >=5.5 all genes","RIN >=5.5 DE genes", "RIN >=5.5 MostVar genes","MEGENA","nGOseq")
file_dir = '~/OneDrive - NextCODE Health/Projects/BI/Results/Classification/FibrosisStage/'
out_dir = '~/OneDrive - NextCODE Health/Projects/BI/Results/Classification/figs/'

all_df = data.frame(name=character(),Method=character(),Metric=character(),
                    Mean.value=double(),sem=double())
out_name = "RIN_gte5.5_noRace_nodups_noXY_comparison"
topGenesInAll = list()
topGenes1000 = list() 
reCalc = F
for (ii in 5 : length(save_names)) {
  vars = readRDS(paste(file_dir,save_names[[ii]],sep=''))
  info = vars$info
  count = 2
  if (reCalc) {
    for (n in 13 : 100) {
      for (m in 1 : 4) {
        true_labels = droplevels(vars$response.lists$test[[m]][[n]]$True.label)
        tmp = vars$response.lists$test[[m]][[n]][,1:6]
        pred_labels = factor(colnames(tmp)[apply(tmp,1,which.max)],levels=levels(true_labels))
        cmt = confusionMatrix(pred_labels,true_labels)
        info[count,"Bal.Acc"] = mean(cmt$byClass[,"Balanced Accuracy"])
        count = count + 2 
      }
    }
  }
  
  varImps = vars$varImps
  mean_info = aggregate(info[,4:7],list(info$datapart,info$method),mean)
  std_info = aggregate(info[,4:7],list(info$datapart,info$method),sem)
  
  for (n in seq(1,8,2)) {
    cat(mean_info[n,"Group.2"],'\t')
    for (m in 3:5) {
      cat(sprintf('%.4f +/- %.4f\t',mean_info[n,m],std_info[n,m]))
    }
    cat('\n')
  }
  tmp_df = melt(mean_info[seq(1,8,2),])[,-1]
  names(tmp_df) = c("Method","Metric","Mean.value")
  tmp_df$sem = melt(std_info[seq(1,8,2),])$value
  tmp_df$name = print_names[[ii]]
  all_df = rbind(all_df,tmp_df)
  
  gene_df = as.data.frame(lapply(varImps,function(x) {cbind(x[[1]])}))
  gene_df2 = as.data.frame(lapply(varImps,function(x) {cbind(x[[3]])}))
  # all_gene_df = cbind(gene_df,gene_df2)
  all_gene_df = gene_df2
  all_genes = unique(unlist(all_gene_df))
  nGenes = length(all_genes)
  nResamples = dim(all_gene_df)[2]
  names(all_gene_df) = paste("Resample",1:nResamples,sep=".")
  gene_ranks = matrix(0,nrow=nGenes,ncol=nResamples)
  rownames(gene_ranks) = all_genes
  for (n in 1 : nGenes) {
    for (m in 1:nResamples) {
      id = which(all_gene_df[,m]==all_genes[[n]])
      gene_ranks[n,m] = ifelse(length(id) == 0, NA, id)
    }
  }
  avg_rank = rowMeans(gene_ranks,na.rm=T)
  avg_rank_nans = rowMeans(gene_ranks)
  sort_avg_rank = sort(avg_rank)
  sort_avg_rank_nans = sort(avg_rank_nans[!is.na(avg_rank_nans)])
  topGenes1000[[ii]] = list(names(sort_avg_rank[1:1000]),sort_avg_rank[1:1000])
  topGenesInAll[[ii]] = list(names(sort_avg_rank_nans),sort_avg_rank_nans)
}
# ggplot(df,aes(fill=Reference,y=value,x=Prediction)) + geom_bar(position='stack',stat='identity'
plot_df <- all_df %>% subset(Metric!="F1")

p = ggplot(plot_df,aes(x=name, y=Mean.value,fill=Method)) + geom_bar(stat="identity", position="dodge") +
  geom_errorbar( aes(ymin=Mean.value-sem, ymax=Mean.value+sem), width=0.4, size=0.4,position=position_dodge(0.9)) +
  labs(x='Dataset',y="Mean value") + facet_wrap(vars(Metric),scales="free",nrow=3) 
# ggsave(paste(out_dir,out_name,'.pdf',sep=''),plot=p,device='pdf',width=8,height=10)

p2 <- ggplot(plot_df,aes(x=Method, y=Mean.value,fill=name)) + geom_bar(stat="identity", position="dodge") +
  geom_errorbar( aes(ymin=Mean.value-sem, ymax=Mean.value+sem), width=0.4, size=0.4,position=position_dodge(0.9)) +
  labs(x='Method',y="Mean value") + facet_wrap(vars(Metric),scales="free",nrow=3)
# ggsave(paste(out_dir,out_name,'_vs_method.pdf',sep=''),plot=p2,device='pdf',width=8,height=10)



