library(dplyr)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(patchwork)
library(ade4)
library(clusterSim)
data<-read.table("new_consented_table_for_analyses.txt", 
                 header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".")


var<-read.table("metadata_8.24.22.txt", header = TRUE,
                sep = '\t')
data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)
data_norm<-as.data.frame(data_norm)

data_norm<-data_norm %>%
  filter(var$Consent_AUG_22 == "Y")
data_norm<-as.matrix(data_norm)

var<-var %>%
  filter(var$Consent_AUG_22 == "Y")

data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)
data_norm <-data_norm %>%
  rownames_to_column("Sample") %>%
  filter(data_norm$mean > 0.0001) %>%
  subset(., select = -c(mean))

data_norm = setNames(data.frame(t(data_norm[,-1])),data_norm[,1])
data_norm<-as.data.frame(t(as.matrix(data_norm)))

data<-data_norm
data<-data[-1,]
dist.JSD <- function(inMatrix, pseudocount=0.00000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

data.dist<-dist.JSD(data)

pam.clustering<-function(x,k) {
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

data.cluster<-pam.clustering(data.dist, k=3)

nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids")

nclusters=NULL

for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index")

data.cluster<-pam.clustering(data.dist, k=2)


obs.silhouette<- mean(silhouette(data.cluster, data.dist)[, 3])
cat(obs.silhouette)

obs.pcoa<-dudi.pco(data.dist, scannf=F, nf=2)
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F)

cluster_list<-cbind(data.cluster,obs.pcoa$li)

source("ANCOM-2.1/scripts/ancom_v2.1.R")

for_ancom<-cbind(var,cluster_list)

data<-read.table("new_consented_table_for_analyses.txt", 
                 header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".")

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Consent_AUG_22 == "Y")
data_norm<-as.matrix(data_norm)


data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0001) #%>%

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Consent_AUG_22 == "Y")
data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

for_ancom$data.cluster<-as.factor(for_ancom$data.cluster)
prepro<-feature_table_pre_process(feature_table = data_filtered,meta_data = for_ancom,sample_var = "xSampleID", group_var = NULL,
                                  out_cut = 0.05, zero_cut = 0.95, lib_cut = 1, neg_lb = FALSE)

feature_table = prepro$feature_table 
meta_data = prepro$meta_data 
struc_zero = prepro$structure_zeros 

ancom_results<-ANCOM(feature_table = feature_table, meta_data = meta_data,
                     struc_zero = struc_zero, main_var = "data.cluster",
                     p_adj_method = "BH", alpha = 0.05)
ancom_results

data<-read.table("new_consented_table_for_analyses.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

data<- data[complete.cases(data), ]
data<-data+1
column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)
data_norm<-as.data.frame(data_norm)

data_norm<-data_norm %>%
  filter(var$Consent_AUG_22 == "Y")
data_norm<-as.matrix(data_norm)

data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)
data_norm <-data_norm %>%
  rownames_to_column("Sample") %>%
  filter(data_norm$mean > 0) %>%
  subset(., select = -c(mean))

data_norm = setNames(data.frame(t(data_norm[,-1])),data_norm[,1])

rownames(ancom_results)<-ancom_results$taxa_id
plotting_data<-feature_table %>%
  filter(ancom_results$detected_0.8 == "TRUE" )

plotting_data<-as.data.frame(t(as.matrix(plotting_data)))
plotting_data<-plotting_data+1
plotting_data<-as.matrix(plotting_data)
melted<-melt(plotting_data)
for_plotting<-cbind(melted,for_ancom)

p1<-ggplot(for_plotting, aes(A1, A2, color = data.cluster)) + 
  geom_point(color= "black", size = 1.5)+ 
  geom_point(aes(color=data.cluster), size=1) + 
  scale_color_manual(values = c("#7570B3", "#E7298A"), name = "PAM Cluster")

p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p4<-p3 + stat_ellipse()
p5_mds<-p4 + theme(axis.title.x = element_text(face = "bold", size = 14),
                   legend.text = element_text(face = "bold", size = 12),
                   axis.title.y = element_text(face = "bold", size = 14), 
                   legend.title  = element_text(face = "bold", size = 14),
                   axis.text = element_text(face = "bold", color = "black")) +
  xlab("PCoA 1") + ylab("PCoA 2")

p5_mds

p1<-ggplot(for_plotting, aes(x = Var2, log(value), fill = data.cluster)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#7570B3", "#E7298A"), name = "PAM Cluster")
p1
p2<-p1 + theme_bw()
p3<-p2 + theme(
  panel.grid.minor = element_blank())
p4<-p3 + labs(x = "ASV", y = "Log10 of ASV Reads")  +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.text.x = element_blank(),
        axis.title = element_text(face = "bold", size = 16),
        axis.ticks.x = element_blank())
p5<-p4 + facet_wrap(~Var2, 
ncol = 5, 
       scales = "free_x", labeller = label_wrap_gen())
p6<-p5 + theme(panel.spacing.x = unit(0.3, "lines"), #panel.border = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 8),
               legend.title = element_text(face = "bold", size = 14),
               legend.text = element_text(face = "bold", size = 12),
               legend.position = "NONE") 
p6

pcombined<-p5_mds+p6
pcombined
pcombined + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", color = "black"))
