# File for making serial dilution plots in a grid of algorithm (rows) metric (columns)
# supersedes lum_ab_plots.R
# make sure the working directory is the home directory of this script
# Authors: Omar Gamel, built on one by Nicholas Cilfone (modified by Richard Li)
# 2018

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

# raw data for annealing type method output as list of matrics for each fraction
# combine ans extract relevant columns
process_rds_file = function(dataraw,meth) {
  #initialize empty dataframe with same columns as dataraw
  data=dataraw[[1]][0,]; colnames(data)[1]="Frac"
  sem=data
  for (i in 1:length(dataraw)){
    data[i,]=colMeans(dataraw[[i]],na.rm=T)
    sem[i,]=apply(dataraw[[i]],2,function(x){std_err_mean(x[!is.na(x)])})#,na.rm=T)
    sem[i,1]=dataraw[[i]]$frac[1]
  }
  
  melt_mean=melt(data,id=1,measure=2:ncol(data))
  melt_sem=melt(sem,id=1,measure=2:ncol(data))
  # remove "tst" and "exptest" but keep "tr" and "val"
  melt_mean=melt_mean[grepl("tr_",melt_mean$variable) | grepl("val_",melt_mean$variable) ,] 
  #melt=melt[! grepl("exptst",melt$variable) & ! grepl("val",melt$variable) ,] 
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

d = '6cancer'
# where annealing-type .RDS files are saved
base_dir = paste('/home/richard/Dropbox-Work/Wuxi/Results/', d, '_splits/preds_for_R/',sep='')

# where classical .RDS files are saved
load_path = paste('/home/richard/Dropbox-Work/Wuxi/Results/', d, '_splits/',sep='')

# colors for plots
full_plotcolors=c("lightblue","pink","blue","red")
overf_plotcolors=full_plotcolors[3:4]
summary_plotcolors = brewer.pal(9,"Paired")

nruns = 50

# Read in Data -- LUMAB -- GENE
# load classical machine learning algorithms
clfile = paste(d,'_pca_cl_info_all.RDS',sep='')
cldataraw = readRDS(paste(load_path,clfile,sep=''))
cl_df = process_cl_file(cldataraw)

# SA Results
sfx = '_nsols20_ntotsols1000'
safile=paste(d,"_split_pca_performance_sa",sfx,"_2.RDS",sep="") # for 6cancer, haave the _2.RDS
#safile=paste(d,"_split_pca_performance_sa",sfx,".RDS",sep="")
sadataraw=readRDS(paste(base_dir,safile,sep='')) 
sa_df = process_rds_file(sadataraw,"SA")

# rand Results
randfile=paste(d,"_split_pca_performance_rand",sfx,"_2.RDS",sep="")
randdataraw=readRDS(paste(base_dir,randfile,sep='')) 
rand_df = process_rds_file(randdataraw,"Random")

# field results
fieldfile=paste(d,"_split_pca_performance_field_2.RDS",sep='')
fielddataraw=readRDS(paste(base_dir,fieldfile,sep='')) 
field_df = process_rds_file(fielddataraw,"Field")

all_cl_df = rbind(cl_df,sa_df,field_df,rand_df)
all_cl_df = all_cl_df %>% filter(Frac!=0.18)

# D-Wave Results
dw_sfx = '_cinit8_nsols20_ntotsols1000'
dwfile=paste(d,"_split_pca_performance_dw",dw_sfx,"_2.RDS",sep='')
dwdataraw=readRDS(paste(base_dir,dwfile,sep='')) 
dw_df = process_rds_file(dwdataraw,"D-Wave")
#dw_df = dw_df %>% filter(Frac!=0.18)

# duplicate the dw rows 
algorithmlist = unique(all_cl_df$Algorithm)
nalgorithms = length(algorithmlist)
ndw = nrow(dw_df)
dw_fulldf = dw_df[rep(1:ndw,each=nalgorithms),] # multiple copies of dw_df, with algorithms inserted
rownames(dw_fulldf) = c()
dw_fulldf$Algorithm=NULL
dw_fulldf = cbind(Algorithm=algorithmlist,dw_fulldf)
dw_fulldf = dw_fulldf[colnames(all_cl_df)]

# merge
fulldata = rbind(all_cl_df,dw_fulldf)
# fulldata=cl_df; full_plotcolors=c("lightblue","blue") # override: plot classic only
fulldata$variable=droplevels(fulldata$variable)

# -------------------- #
# train and test metrics plot
nrowplot = length(unique(fulldata$Algorithm)); ncolplot = nlevels(fulldata$variable)
p_all = ggplot(data = fulldata, aes(x=Frac, y=value)) + 
    geom_point(aes(colour = Dataset,shape=Dataset)) + 
    geom_line(aes(colour = Dataset)) +
    geom_errorbar(aes(colour = Dataset,linetype='dashed',width=0.015,
                      ymin=(value-sem),ymax=(value+sem))) +
    scale_color_manual(values=full_plotcolors) + 
    facet_wrap(Algorithm~variable, nrow = nrowplot, ncol = ncolplot) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    labs(x='Fraction of Training Data (80%)', y='Metric')

ggsave(paste(load_path, 'figs/', 'all_classic_dw_sem_valid_', as.character(d), '.pdf', sep=''), 
       plot = p_all, device = 'pdf', width = 25, height = 15)

# -------------------- #
# balanced accuracy test -summary plot
dw_df['Algorithm']="D-Wave"
all_data_df=rbind(all_cl_df,dw_df)

drop_fracs = c(0.05,0.15)
#drop_fracs = 0.06
#drop_fracs = 0.18
all_data_df$Frac = as.factor(all_data_df$Frac)
all_data_df = all_data_df %>% subset(!(Frac %in% drop_fracs))
frac = all_data_df$Frac
all_data_df$Frac = as.numeric(levels(frac))[frac]

#choose only Validation balanced accuracy, then remove columns that no longer vary
#didn't build it from fulldata since it has duplicates
baccdata = all_data_df %>% 
  subset(variable == "Bal.Accuracy" & Dataset == "Valid", select=-c(variable,Dataset)) %>%
  arrange(Algorithm,Frac)
baccdata$Algorithm = factor(baccdata$Algorithm,levels=sort(unique(baccdata$Algorithm)))
p_bacc = ggplot(data = baccdata, aes(x=Frac, y=value)) + 
    geom_point(aes(colour = Algorithm,shape=Algorithm),size=3) + scale_shape_manual(values=c(5:8,15:19)) + 
    geom_line(aes(colour = Algorithm),size=0.75) + 
    geom_errorbar(aes(colour = Algorithm,width=0.01,
                      ymin=(value-sem),ymax=(value+sem)),size=0.75) + 
    scale_color_manual(values=summary_plotcolors) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text=element_text(size=14),
          text=element_text(size=14)) +
    labs(x='Fraction of Training Data', y='Balanced Accuracy')

ggsave(paste(load_path, 'figs/', 'all_classic_dw_bacc_', as.character(d),dw_sfx, '.pdf', sep=''), plot = p_bacc, device = 'pdf', width = 12, height = 7)


all_bacc_stats_df = all_data_df %>% filter(Dataset=="Valid",variable=="Bal.Accuracy")
cl_bacc_pop = cldataraw[,c(1,2,17)]
dw_bacc_pop = data.frame()
for (n in 1 : length(dwdataraw)) {
  dw_bacc_pop = rbind(dw_bacc_pop,dwdataraw[[n]][,c(1,16)])
}
dw_bacc_pop$method = 'D-Wave'
all_bacc_pop = rbind(cl_bacc_pop,dw_bacc_pop)

saveRDS(all_bacc_pop,paste(load_path,d,'_all_bacc_pop_data',dw_sfx,'.RDS',sep=''))
saveRDS(all_bacc_stats_df,paste(load_path,d,'_bacc_stats',dw_sfx,'.RDS',sep=''))

# -------------------- #
# overfitting plot
cl_overf= all_cl_df %>% 
  filter(Frac != 1) %>% 
  group_by(Frac, Algorithm,variable) %>% 
  summarise(over=value[Dataset=="Train"]-value[Dataset=="Valid"],
            sem=sqrt(sum(sem^2)),Dataset="Tr-Val") %>% 
  as.data.frame()

dw_overf= dw_fulldf %>% 
  filter(Frac != 1) %>% 
  group_by(Frac, Algorithm,variable) %>% 
  summarise(over=value[Dataset=="Train"]-value[Dataset=="Valid"],
            sem=sqrt(sum(sem^2)),Dataset="Tr-Val.DW") %>% 
  as.data.frame()

full_overf = rbind(cl_overf,dw_overf)
nrowplot = length(unique(cl_overf$Algorithm))
ncolplot = nlevels(fulldata$variable)

p_overf = ggplot(data = full_overf, aes(x=Frac, y=over)) + geom_point(aes(colour = Dataset)) + 
  geom_line(aes(colour = Dataset)) + 
  geom_errorbar(aes(colour = Dataset,
                    linetype='dashed',
                    width=0.015,
                    ymin=(over-sem),
                    ymax=(over+sem))) + 
  scale_color_manual(values=overf_plotcolors) + 
  facet_wrap(Algorithm~variable, nrow = nrowplot, ncol = ncolplot) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  labs(x='Fraction of Training Data (80%)', y='Metric Overfit (Train - Test)')

ggsave(paste(load_path, 'figs/', 'all_classic_dw_overfit_', as.character(d), '.pdf', sep=''),
       plot = p_overf, device = 'pdf', width = 25, height = 15)

# -------------------- #
# balanced accuracy overfitting - summary plot
#choose only test balanced accuracy, then remove columns that no longer vary
dw_overfbacc = dw_overf %>% 
  subset(variable=="Bal.Accuracy" & Algorithm == "LASSO",select=-c(variable,Dataset))
dw_overfbacc['Algorithm']="D-Wave"

cl_overfbacc = cl_overf %>% 
  subset(variable=="Bal.Accuracy",select=-c(variable,Dataset))

overfbaccdata=rbind(cl_overfbacc,dw_overfbacc)
overfbaccdata=overfbaccdata[order(overfbaccdata$Algorithm,overfbaccdata$Frac),c(2,1,3,4)]
rownames(overfbaccdata)=NULL

p_overfbacc = ggplot(data = overfbaccdata, aes(x=Frac, y=over)) + 
  geom_point(aes(colour = Algorithm)) + 
  geom_line(aes(colour = Algorithm)) +
  geom_errorbar(aes(colour = Algorithm, linetype='dashed',
                    width=0.015,ymin=(over-sem),ymax=(over+sem))) +
  scale_color_manual(values=summary_plotcolors) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x='Fraction of Training Data (80%)', y='Balanced Accuracy Train-Test (Overfitting)')
#ggsave(paste(base_dir, 'figs/', 'all_classic_dw_overfbacc_', as.character(d), '.pdf', sep=''), plot = p_overfbacc, device = 'pdf', width = 25, height = 15)

