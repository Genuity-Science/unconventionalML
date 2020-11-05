# Script for generating plots for serial dilutions. 

rm(list = ls())
gc()

library(glmnet)
library(data.table)
library(methods)
library(feather)
library(caret)
library(RColorBrewer)
library(dplyr)

# standard error of the mean
std_err_mean = function(x) sd(x)/sqrt(length(x))

# raw data for annealing-type method output as list of matrics for each fraction
# combine ans extract relevant columns
process_rds_file = function(dataraw,meth) {
  #initialize empty dataframe with same columns as dataraw
  data=dataraw[[1]][0,]; colnames(data)[1]="Frac"
  sem=data
  for (i in 1:length(dataraw)){
    # if (class(dataraw[[i]]) == 'data.frame') dataraw[[i]] = as.matrix(dataraw[[i]])
    data[i,]=colMeans(dataraw[[i]],na.rm=T)
    sem[i,]=apply(dataraw[[i]],2,function(x){std_err_mean(x[!is.na(x)])})
    sem[i,1]=dataraw[[i]]$frac[1]
  }
  
  melt_mean=melt(data,id=1,measure=2:ncol(data))
  melt_sem=melt(sem,id=1,measure=2:ncol(data))
  # remove "tst" and "exptest" but keep "tr" and "val"
  melt_mean=melt_mean[grepl("tr_",melt_mean$variable) | grepl("val_",melt_mean$variable) ,] 
  melt_sem=melt_sem[grepl("tr_",melt_sem$variable) | grepl("val_",melt_sem$variable) ,] 
  rownames(melt_mean)=rownames(melt_sem)=c();
  colnames(melt_sem)[colnames(melt_sem)=='value']='sem'
  df=merge(melt_mean,melt_sem,by=c("Frac","variable"))
  
  # split columns like "tr_acc" to two columns
  extra_cols = data.frame(do.call('rbind', strsplit(as.character(df$variable),'_',fixed=TRUE)))
  colnames(extra_cols)=c("Dataset","variable")
  extra_cols$Dataset = as.factor(extra_cols$Dataset)
  extra_cols$variable = as.factor(extra_cols$variable);
  levels(extra_cols$Dataset)=c("Train","Valid") # alphabetical tr then val
  if (nlevels(extra_cols$variable) ==7) {
    levels(extra_cols$variable)=c("Accuracy","AUPRC","AUC","Bal.Accuracy","F1","Precision","Recall")
  }else if (nlevels(extra_cols$variable) == 6) {
    levels(extra_cols$variable)=c("Accuracy","AUC","Bal.Accuracy","F1","Precision","Recall")
  }
  df$variable=NULL
  df=cbind(df,extra_cols)
  df=df[(df$variable %in% c("Accuracy","AUC","Bal.Accuracy","F1")),] #keep only four variables
  df$Algorithm=meth
  
  return(df)
}

# function to read in data from classical file. Uses info as list from an output
# script (e.g. bootstrap_resamples_output.rds)
process_cl_file = function(info){
  new_info = melt(as.data.table(info), id=c("frac","method"), 
                  measure=patterns(Accuracy="_acc$",Bal.Accuracy="_bacc$",AUC="_auroc$",
                                  F1="_F1$"),
                  value.factor=F)
  
  names(new_info)[1:3] = c("Frac","Algorithm","Dataset")
  levels(new_info$Dataset) = c("Train","Test","Valid","ExpandTest")
  new_info = new_info %>% filter(Dataset=="Train" | Dataset == "Valid")
  
  mean_stats = aggregate(.~Frac+Algorithm+Dataset,new_info,mean,na.action=na.pass)
  mean_stats$Algorithm = toupper(mean_stats$Algorithm)
  mean_stats = melt(mean_stats,id=1:3,measure=4:ncol(mean_stats))
  rownames(mean_stats) = c()
  
  sem_stats = aggregate(.~Frac+Algorithm+Dataset,new_info,std_err_mean,na.action=na.pass)
  sem_stats$Algorithm = toupper(sem_stats$Algorithm)
  sem_stats = melt(sem_stats,id=1:3,measure=4:ncol(sem_stats))
  rownames(sem_stats) = c()
  
  cl_df = merge(mean_stats,sem_stats,by=c("Frac","Algorithm","Dataset","variable"))
  colnames(cl_df)[5:6] = c("value","sem")
  
  cl_df$Algorithm[cl_df$Algorithm=="RIDGE"] = "Ridge"
  
  return(cl_df)
}

d = 'ERpn'
type = 'pca'
# where annealing-type .RDS files are saved
base_dir = paste('/Users/rli/OneDrive - Genuity Science/Projects/quantumML/Classification/output/',sep='')

# where classical .RDS files are saved
load_path = paste('/Users/rli/OneDrive - Genuity Science/Projects/quantumML/Classification/output/',sep='')

# colors for plots
summary_plotcolors = brewer.pal(10,"Paired")

nruns = 50

# Read in Data -- LUMAB -- GENE
# load classical machine learning algorithms

clfile = paste(d,'_sd_', type,'_cl_all_info.RDS',sep='')
cldataraw = readRDS(paste(load_path,clfile,sep=''))
cl_df = process_cl_file(cldataraw)

# SA Results
sfx = '_nsols20_ntotsols1000'
safile=paste(d,'_',type,"_split_performance_sa",sfx,".RDS",sep="")
sadataraw=readRDS(paste(base_dir,safile,sep='')) 
sa_df = process_rds_file(sadataraw,"SA")

# rand Results
randfile=paste(d,'_',type,"_split_performance_rand",sfx,".RDS",sep="")
randdataraw=readRDS(paste(base_dir,randfile,sep='')) 
rand_df = process_rds_file(randdataraw,"Random")

# field results
fieldfile=paste(d,'_',type,"_split_performance_field.RDS",sep='')
fielddataraw=readRDS(paste(base_dir,fieldfile,sep='')) 
field_df = process_rds_file(fielddataraw,"Field")

all_cl_df = rbind(cl_df,sa_df,field_df,rand_df)

# D-Wave Results
dw_sfx = '_cinit8_nsols20_ntotsols1000'
dwfile=paste(d,'_',type,"_split_performance_dw",dw_sfx,".RDS",sep='')
dwdataraw=readRDS(paste(base_dir,dwfile,sep='')) 
dw_df = process_rds_file(dwdataraw,"D-Wave")

# RBM Results
rbmfile=paste(d,'_',type,"_split_performance_RBM.RDS",sep='')
rbmdataraw=readRDS(paste(base_dir,rbmfile,sep='')) 
rbmdataraw = lapply(rbmdataraw,function(x){x[-2]}) # remove method
rbm_df = process_rds_file(rbmdataraw,"RBM")


# -------------------- #
dw_df['Algorithm']="D-Wave"
all_data_df=rbind(all_cl_df,rbm_df,dw_df)

drop_fracs = 0.18 #c(0.05,0.15)
all_data_df$Frac = as.factor(all_data_df$Frac)
all_data_df = all_data_df %>% subset(!(Frac %in% drop_fracs))

frac = all_data_df$Frac
all_data_df$Frac = as.numeric(levels(frac))[frac]

#choose only Validation balanced accuracy, then remove columns that no longer vary
baccdata = all_data_df %>% 
  subset(variable == "Bal.Accuracy" & Dataset == "Valid", select=-c(variable,Dataset)) %>%
  arrange(Algorithm,Frac)
baccdata$Algorithm = factor(baccdata$Algorithm,levels=sort(unique(baccdata$Algorithm)))
p_bacc = ggplot(data = baccdata, aes(x=Frac, y=value)) + 
    geom_point(aes(colour = Algorithm,shape=Algorithm),size=3) + scale_shape_manual(values=c(5:9,15:19)) + 
    geom_line(aes(colour = Algorithm),size=0.75) + 
    geom_errorbar(aes(colour = Algorithm,width=0.01,
                      ymin=(value-sem),ymax=(value+sem)),size=0.75) + 
    scale_color_manual(values=summary_plotcolors) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text=element_text(size=14),
          text=element_text(size=14)) +
    labs(x='Fraction of Training Data', y='Balanced Accuracy')

#p_bacc + theme(legend.position="bottom",legend.direction = "horizontal") + scale_x_continuous(expand = c(0,0.02)) + theme(text=element_text(size=18),legend.title=element_text(size=16))

# ggsave(paste(load_path, d,'_',type,'_serial_dilutions_withRBM_v2.pdf', sep=''), plot = p_bacc, device = 'pdf', width = 8, height = 5)


all_bacc_stats_df = all_data_df %>% filter(Dataset=="Valid",variable=="Bal.Accuracy")
cl_bacc_pop = cldataraw[,c("frac","method","val_bacc","tr_bacc")]
cl_bacc_pop$method = toupper(cl_bacc_pop$method)
cl_bacc_pop$method[cl_bacc_pop$method=="RIDGE"] = "Ridge"
dw_bacc_pop = data.frame()
sa_bacc_pop = data.frame()
field_bacc_pop = data.frame()
rand_bacc_pop = data.frame()
rbm_bacc_pop = data.frame()
for (n in 1 : length(dwdataraw)) {
  dw_bacc_pop = rbind(dw_bacc_pop,dwdataraw[[n]][,c("frac","val_bacc","tr_bacc")])
  sa_bacc_pop = rbind(sa_bacc_pop,sadataraw[[n]][,c("frac","val_bacc","tr_bacc")])
  field_bacc_pop = rbind(field_bacc_pop,fielddataraw[[n]][,c("frac","val_bacc","tr_bacc")])
  rand_bacc_pop = rbind(rand_bacc_pop,randdataraw[[n]][,c("frac","val_bacc","tr_bacc")])
  rbm_bacc_pop = rbind(rbm_bacc_pop,rbmdataraw[[n]][,c("frac","val_bacc","tr_bacc")])
}
dw_bacc_pop$method = 'D-Wave'
sa_bacc_pop$method = 'SA'
rand_bacc_pop$method = 'Random'
field_bacc_pop$method = 'Field'
rbm_bacc_pop$method = 'RBM'
all_bacc_pop = rbind(cl_bacc_pop,dw_bacc_pop,sa_bacc_pop,field_bacc_pop,rand_bacc_pop,rbm_bacc_pop)
#saveRDS(all_bacc_pop,paste(load_path,d,'_all_bacc_pop_data',dw_sfx,'.RDS',sep=''))
#saveRDS(all_bacc_stats_df,paste(load_path,d,'_bacc_stats',dw_sfx,'.RDS',sep=''))

# -------------------- #
# overfitting plot
# no test data for Frac = 1

# algorithmlist = unique(all_cl_df$Algorithm)
# nalgorithms = length(algorithmlist)
# ndw = nrow(dw_df)
# dw_fulldf = dw_df[rep(1:ndw,each=nalgorithms),] # multiple copies of dw_df, with algorithms inserted
# rownames(dw_fulldf) = c()
# dw_fulldf = cbind(dw_fulldf,Algorithm=algorithmlist)
# dw_fulldf = dw_fulldf[colnames(all_cl_df)]

# cl_overf= all_cl_df %>%
#   group_by(Frac, Algorithm,variable) %>%
#   summarise(over=value[Dataset=="Train"]-value[Dataset=="Valid"],sem=sqrt(sum(sem^2)),Dataset="Tr-Val") %>%
#   as.data.frame()
# dw_overf= dw_fulldf %>%
#   group_by(Frac, Algorithm,variable) %>%
#   summarise(over=value[Dataset=="Train.DW"]-value[Dataset=="Valid.DW"],sem=sqrt(sum(sem^2)),Dataset="Tr-Val.DW") %>%
#   as.data.frame()
# 
# full_overf = rbind(cl_overf,dw_overf)
# 
# p_overf = ggplot(data = full_overf, aes(x=Frac, y=over)) + 
#   geom_point(aes(colour = Dataset)) + 
#   geom_line(aes(colour = Dataset)) + 
#   geom_errorbar(aes(colour = Dataset,linetype='dashed',width=0.015,ymin=(over-sem),ymax=(over+sem))) +
#   scale_color_manual(values=overf_plotcolors) + 
#   facet_wrap(Algorithm~variable, nrow = nrowplot, ncol = ncolplot) + 
#   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   labs(x='Fraction of Training Data (80%)', y='Metric Overfit (Train - Test)')
#ggsave(paste(base_dir, 'figs/', 'all_classic_dw_overfit_', as.character(name), '.pdf', sep=''), plot = p_overf, device = 'pdf', width = 25, height = 15)

all_data_bacc = all_data_df %>% 
  filter(variable=="Bal.Accuracy") %>% 
  group_by(Frac,Algorithm)
all_data_bacc[all_data_bacc$Dataset=="Valid","Dataset"] = "Test" 


tmp_df = all_bacc_pop %>% subset(val_bacc > 0.7) %>% group_by(method,frac)
# keep so can find p-values 

tmp_df2 = tmp_df %>% subset(val_bacc < tr_bacc) %>% 
  summarise(mean_tr=mean(tr_bacc), 
            mean_tst = mean(val_bacc),
            sem_tr = sd(tr_bacc)/sqrt(length(tr_bacc)), 
            sem_tst = sd(val_bacc)/sqrt(length(val_bacc))
            ) 

mean_diff_df = all_data_bacc %>% group_by(Algorithm) %>% 
  summarise(overf=mean(value[Dataset=="Train"] - value[Dataset=="Test"]))
lvs = mean_diff_df$Algorithm[order(mean_diff_df$overf)]
  
# mean_diff_df = tmp_df2 %>%
#   summarise(overf=mean(mean_tr-mean_tst))
# lvs = mean_diff_df$method[order(mean_diff_df$overf)]
  
all_data_bacc$Algorithm = factor(all_data_bacc$Algorithm,levels=lvs)

p_overf = ggplot(data = all_data_bacc,aes(x=Frac,y=value)) + 
  geom_point(aes(color=Dataset)) + #,size=3) + 
  geom_line(aes(color=Dataset)) + #,size=1) + 
  geom_errorbar(aes(color=Dataset,linetype='dashed',width=0.025,ymin=(value-sem),ymax=value+sem),size=0.75) + 
  scale_color_manual(values=c("blue","red")) + #scale_y_continuous(breaks=c(0.5,0.75,1)) + scale_x_continuous(breaks=c(0.2,0.6,1.0),labels=c("0.2","0.6","1")) + 
  facet_wrap(~Algorithm,nrow=3) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x="Fraction of Training Data",y="Balanced Accuracy") +
  guides(linetype=FALSE,size=FALSE) + theme(text=element_text(size=18), panel.spacing.x=unit(2.5,"lines")) #+ 
  #theme(legend.position="bottom",legend.direction = "horizontal")
show(p_overf)
# ggsave(paste(load_path,'all_bacc_overfit_',as.character(d),'_',type,'_withRBM_v2.pdf',sep=''),plot=p_overf,device='pdf',width=12,height=6)


#ERpn gene: compare NB to D-Wave at 0.02
#lumAB gene: compare SVM and RF to D-Wave at 0.05 
#lumAB: compare SVM to D-Wave at 0.2 
#ERpn: compare SVM to D-Wave at 0.1
p_df = all_bacc_pop %>% subset(frac==0.2 & method %in% c("RBM","SVM")) %>% mutate(diff = tr_bacc - val_bacc)
res <- wilcox.test(tr_bacc - val_bacc ~ method, data = p_df, paired = FALSE)
print(res$p.value)
overf_stats = p_df %>% group_by(method) %>% summarise(mean=mean(diff),sem=sd(diff)/sqrt(length(diff)))
print(overf_stats)

res <- wilcox.test(val_bacc ~ method,data = p_df,paired=F)
print(res$p.value)
st = p_df %>% group_by(method) %>% summarise(mean=mean(val_bacc),sem=sd(val_bacc)/sqrt(length(val_bacc)))
print(st)
# p_df = tmp_df %>% subset(frac==0.05 & method %in% c("D-Wave","RF"))
# res <- wilcox.test(tr_bacc - val_bacc ~ method, data = p_df, paired = FALSE)
# print(res$p.value)
