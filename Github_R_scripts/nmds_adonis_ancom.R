library(ggplot2)
library(dplyr)
library(tidyverse)
library(vegan)
library(dunn.test)
library(rcartocolor)
library(RColorBrewer)
library(pals)
library(patchwork)
library(ggplotify)

data<-read.table("new_rarefied_table_for_analyses.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

var<-read.table("metadata_8.24.22.txt", header = TRUE, row.names = 1, sep ='\t')

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)
data_norm<-as.data.frame(data_norm)

data_norm<-data_norm %>%
  filter(var$Paired == "Y")
data_norm<-as.matrix(data_norm)

var<-var %>%
  filter(var$Paired == "Y")

data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)
data_norm <-data_norm %>%
  rownames_to_column("Sample") %>%
  filter(data_norm$mean > 0) %>%
  subset(., select = -c(mean))

data_norm = setNames(data.frame(t(data_norm[,-1])),data_norm[,1])

data_nmds <- metaMDS(data_norm, distance = 'bray', autotransform = F, 
                     k = 2, noshare = F, trymax = 100, parallel = 6)
nmds_scores<-scores(data_nmds)
nmds_df<-as.data.frame(nmds_scores)

detach(var)
attach(var)

for_plotting<-cbind(nmds_df, var)

p1<-ggplot(for_plotting, aes(NMDS1, NMDS2, color = Participant)) + 
  geom_point(color= "black", size = 1.5)+ 
  geom_point(aes(color=Participant), size=1) + 
  scale_color_brewer(palette = "Dark2")

p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p4<-p3 + stat_ellipse()
p5_mds<-p4 + theme(axis.title.x = element_text(face = "bold", size = 14),
                   legend.text = element_text(face = "bold", size = 12),
                   axis.title.y = element_text(face = "bold", size = 14), 
                   legend.title  = element_text(face = "bold", size = 14))
p5_mds

bray<- vegdist(data_norm, method = 'bray')

ad<-adonis(data_norm ~ Participant + Pair, 
           data = var, parallel = 4, 
           method = "bray")
ad

shannon<-diversity(data_norm, index = "shannon")

for_plotting<-cbind(nmds_df,var,shannon)

shapiro.test(shannon)
hist(shannon)
qqnorm(shannon)
shannon_krusk<-kruskal.test(shannon~Participant, data = for_plotting)
shannon_krusk
shannon_krusk<-kruskal.test(shannon~Pair, data = for_plotting)
shannon_krusk

###ANCOM###
source("ANCOM-2.1/scripts/ancom_v2.1.R")

data<-read.table("new_consented_table_for_analyses.txt", 
                 header = TRUE, sep = '\t',row.names = 1, quote = "", dec = ".")


var<-read.table("metadata_8.24.22.txt", header = TRUE,
                sep = '\t')
data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Paired == "Y")
data_norm<-as.matrix(data_norm)


data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0001) #%>%

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Paired == "Y")
data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

var<-var %>%
  filter(var$Paired == "Y")
rownames(var)<-var$xSampleID

prepro<-feature_table_pre_process(feature_table = data_filtered,meta_data = var,sample_var = "xSampleID", group_var = NULL,
                                  out_cut = 0.05, zero_cut = 0.95, lib_cut = 1, neg_lb = FALSE)

feature_table = prepro$feature_table 
meta_data = prepro$meta_data 
struc_zero = prepro$structure_zeros 

ancom_results<-ANCOM(feature_table = feature_table, meta_data = meta_data,
                    struc_zero = struc_zero, main_var = "Participant",
                   p_adj_method = "BH", alpha = 0.05)
ancom_results

data<-read.table("new_consented_table_for_analyses.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

var<-read.table("metadata_8.24.22.txt", header = TRUE, row.names = 1, sep ='\t')

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)
data_norm<-as.data.frame(data_norm)

data_norm<-data_norm %>%
  filter(var$Consent_AUG_22 == "Y")
data_norm<-as.matrix(data_norm)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Paired == "Y")
data_norm<-as.matrix(data_norm)

var<-var %>%
 filter(var$Paired == "Y")

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
  filter(ancom_results$detected_0.9 == "TRUE" )

plotting_data<-as.data.frame(t(as.matrix(plotting_data)))
plotting_data<-plotting_data+1
plotting_data<-as.matrix(plotting_data)
melted<-melt(plotting_data)
for_plotting<-cbind(melted,var)

p1<-ggplot(for_plotting, aes(x = Var2, log(value), fill = Participant)) + 
  geom_boxplot() + 
  scale_fill_brewer(palette = "Dark2")
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
p6<-p5 + theme(panel.spacing.x = unit(0.3, "lines"), 
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 8),
               legend.title = element_text(face = "bold", size = 14),
               legend.text = element_text(face = "bold", size = 12),
               legend.position = "NONE") #+

p6

pcombined<-p5_mds+p6
pcombined
pcombined + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", color = "black"))
