# Script that builds the latest heatmaps as of April 2019, based on clustering the rows (i.e. samples) of lumAB
# Try multiple strategies to make the heatmap. Strategy 1 is the chosen one, where two heatmaps are made based on 
# clustering done in python - one uses quantum heirarchical clustering (qHCl) and the other uses a classical 
# "custom Ward" HCl that is identical to the qHCl except for the minimization part
# Original script Nov 2018
# Author: Omar Gamel
# Edited: Sharvari Gujja

rm(list=ls())
gc()

library(RColorBrewer)
library(pheatmap)
library(dynamicTreeCut)
library(ape)
library(flashpcaR)
library(dendextend)
library(feather)
library(WGCNA)
library(glmnet)
library(caret)
library(seriation)
library(gplots)
library(ggplot2)

n_pcs=44 # number of principal components
figpath="/Users/sgujja/OneDrive - NextCODE Health/Omar_Code/Code/QISKit/HClustering_reproducibility_lumAB_44genes/figs/rows_heatmaps/"
hmcols = rev(colorRampPalette(brewer.pal(11,"RdBu"))(11))
blueyellow = function(n){colorpanel(n,"blue","white","yellow")}
barcolors=c("#00FFFFFF","#FF0000FF")

#### Functions
mAnnotatedHeatmap = function(X_in,y_in,clust_info=NULL,title="data",labelsize=0.2,colorHeight=1,dist_method="Pearson", linkage="average", 
                             key=TRUE,keysize=1.5, cap=3.5, key.par=list(cex=1.3)){
  # distance between rows (i.e. samples), pearson correlation distance  or else Euclidean
  dist_matrix = if (dist_method == "Pearson") as.dist((1-cor(t(X_in)))/2) else dist(X_in)
  
  clust = hclust(dist_matrix, method=linkage)
  # if clust_info is passed, use its properties to override those of the calculated hclust
  if (!is.null(clust_info)){
    for (property in names(clust_info)){
      clust[[property]]=clust_info[[property]]
    } 
  }
  
  #reorder clust to keep same labels in y_in adjacent
  ylabels=y_in[,1] #vector
  ycolors=ylabels; levels(ycolors)=barcolors
  n_labelclasses=nlevels(ylabels)
  labelnums = ylabels
  levels(labelnums)=1:n_labelclasses # change labels to sequential integers numbers
  label_dist = dist(labelnums) #distance within each label class is zero
  clust = reorder(clust,label_dist,method="OLO") #"OLO" or "GW"
  
  ### PLOT ###
  X_plot = X_in
  if (is.numeric(cap)){ #cap is the absolute value limit, use it to cap
    X_plot[X_plot> cap]= cap; X_plot[X_plot< -cap]= -cap
  }
  
  pdf(paste0(figpath, title, '.pdf'),width = 25, height = 15)
  hmap=heatmap.2(as.matrix(t(X_plot)), col=redblue, symbreak=TRUE, trace='none', cexRow=0.3, ColSideColors=as.character(ycolors), Colv=as.dendrogram(clust), key=key, keysize=keysize,key.par=key.par,labRow = FALSE, labCol = FALSE)
  dev.off()
  
  return(list(clust,hmap))
}

importClustOverride = function(algorithm="c"){
  
  path = "temp/%s_%s.feather"
  path = sprintf(path,"%s",algorithm)
  
  clust_info = vector(mode="list",length=5)
  names(clust_info) = c("merge","height","order","method","call")
  
  clust_info[[1]]=unname(as.matrix(read_feather(sprintf(path,"merge"))))
  clust_info[[2]]=read_feather(sprintf(path,"height"))[[1]]
  clust_info[[3]]=read_feather(sprintf(path,"order"))[[1]]
  
  clust_info[[4]]=ifelse(algorithm=="c","custom ward","quantum HCl")
  clust_info[[5]]=ifelse(algorithm=="c","python classical","python Durr Hoyer")
  
  return(clust_info)
}

########### 
# load data

datasetname = "lumAB"; rawdata=as.data.frame(read_feather("/Users/sgujja/OneDrive - NextCODE Health/Omar_Code/Code/data/all/data5_lumAB_train_normalized.feather"))

X = rawdata[2:ncol(rawdata)]
y = rawdata[1]
ylabels=y$cancer

pc=flashpca(as.matrix(X), ndim = n_pcs, stand = "sd", do_loadings = TRUE)


# Strategy 1: top n_genes genes from PC1
s=1
n_genes=44
ord=order(abs(pc$loadings[,1]),decreasing = TRUE)[1:n_genes] #rankings of top n_genes from PC1
X_pc1 = X[,ord]
hmap_pc1=mAnnotatedHeatmap(X_pc1,y,title=paste(paste0(datasetname,s),"pc1",n_genes,"genes",sep="_"))
saveRDS(cbind(y,X_pc1),file=paste0(datasetname,"_pc1_",n_genes,".RDS"))
write_feather(cbind(y,X_pc1),path=paste0(datasetname,"_pc1_",n_genes,".feather"))

# Strategy 1Q: top n_genes genes from PC1 --> compared custom classic and quantum
hmap_c=mAnnotatedHeatmap(X_pc1,y,clust_info=importClustOverride("c"),title=paste(paste0(datasetname,s),"pc1",n_genes,"genes","customWardHCl",sep="_"))
hmap_q=mAnnotatedHeatmap(X_pc1,y,clust_info=importClustOverride("q"),title=paste(paste0(datasetname,s),"pc1",n_genes,"genes","quantumHCl",sep="_"),key=FALSE)
