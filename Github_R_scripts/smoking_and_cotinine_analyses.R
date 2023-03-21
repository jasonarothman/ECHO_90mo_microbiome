library(ggplot2)
library(dplyr)
library(tidyverse)
library(vegan)
library(rcartocolor)
library(RColorBrewer)
library(lmerTest)
library(patchwork)
library(ggplotify)
library(Maaslin2)
library(reshape2)

data<-read.table("new_rarefied_table_for_analyses.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

var<-read.table("metadata_8.24.22.txt", header = TRUE, row.names = 1, sep ='\t')

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)
data_norm<-as.data.frame(data_norm)

data_norm<-data_norm %>%
  filter(var$Cotinine_concentration_categorical != "NA")
data_norm<-as.matrix(data_norm)

var<-var %>%
  filter(var$Cotinine_concentration_categorical != "NA")

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
brewer.pal(8,"Dark2")

p1<-ggplot(for_plotting, aes(NMDS1, NMDS2, color = Cotinine_concentration_categorical)) + 
  geom_point(color= "black", size = 1.5)+ 
  geom_point(aes(color=Cotinine_concentration_categorical), size=1) + 
  scale_color_manual(values = c("#66A61E", "#E6AB02"), name = "Salivary Cotinine")

p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p4<-p3 + stat_ellipse()
p5_mds_child<-p4 + theme(axis.title.x = element_text(face = "bold", size = 14),
                   legend.text = element_text(face = "bold", size = 12),
                   axis.title.y = element_text(face = "bold", size = 14), 
                   legend.title  = element_text(face = "bold", size = 14),
                   axis.text = element_text(face = "bold", color = "black"))
p5_mds_child

bray<- vegdist(data_norm, method = 'bray')

ad<-adonis(data_norm ~ var$Cotinine_concentration_categorical, 
           data = var, parallel = 4, 
           method = "bray")
ad

shannon<-diversity(data_norm, index = "shannon")

for_plotting<-cbind(nmds_df,var,shannon)

shapiro.test(shannon)
hist(shannon)
qqnorm(shannon)
shannon_krusk<-kruskal.test(shannon~Cotinine_concentration_categorical, data = for_plotting)
shannon_krusk

###ANCOM child smoking###
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
  filter(var$Cotinine_concentration_categorical != "NA")
data_norm<-as.matrix(data_norm)


data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0001) #%>%

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Cotinine_concentration_categorical != "NA")
data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

var<-var %>%
  filter(var$Cotinine_concentration_categorical != "NA")
rownames(var)<-var$xSampleID

prepro<-feature_table_pre_process(feature_table = data_filtered,meta_data = var,sample_var = "xSampleID", group_var = NULL,
                                  out_cut = 0.05, zero_cut = 0.95, lib_cut = 1, neg_lb = FALSE)

feature_table = prepro$feature_table 
meta_data = prepro$meta_data 
struc_zero = prepro$structure_zeros 

ancom_results<-ANCOM(feature_table = feature_table, meta_data = meta_data,
                     struc_zero = struc_zero, main_var = "Cotinine_concentration_categorical",
                     p_adj_method = "BH", alpha = 0.05)
ancom_results

data<-read.table("new_rarefied_table_for_analyses.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1)

var<-read.table("metadata_8.24.22.txt", header = TRUE, row.names = 1, sep ='\t')

data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)
data_norm<-as.data.frame(data_norm)

data_norm<-data_norm %>%
  filter(var$Participant == "Caregiver", var$Smoker != "NA")
data_norm<-as.matrix(data_norm)

var<-var %>%
  filter(var$Participant == "Caregiver", var$Smoker != "NA")

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
brewer.pal(8,"Dark2")

p1<-ggplot(for_plotting, aes(NMDS1, NMDS2, color = Smoker)) + 
  geom_point(color= "black", size = 1.5)+ 
  geom_point(aes(color=Smoker), size=1) + 
  scale_color_manual(values = c("#66A61E", "#E6AB02"), name = "Smoker")

p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p4<-p3 + stat_ellipse()
p5_mds_parent<-p4 + theme(axis.title.x = element_text(face = "bold", size = 14),
                   legend.text = element_text(face = "bold", size = 12),
                   axis.title.y = element_text(face = "bold", size = 14), 
                   legend.title  = element_text(face = "bold", size = 14),
                   axis.text = element_text(face = "bold", color = "black"))
p5_mds_parent

bray<- vegdist(data_norm, method = 'bray')

ad<-adonis(data_norm ~ var$Smoker, 
           data = var, parallel = 4, 
           method = "bray")
ad

shannon<-diversity(data_norm, index = "shannon")

for_plotting<-cbind(nmds_df,var,shannon)

shapiro.test(shannon)
hist(shannon)
qqnorm(shannon)

shannon_krusk<-kruskal.test(shannon~Smoker, data = for_plotting)
shannon_krusk

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
  filter(var$Participant == "Caregiver", var$Smoker != "NA")
data_norm<-as.matrix(data_norm)


data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0001) 

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Participant == "Caregiver", var$Smoker != "NA")
data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

var<-var %>%
  filter(var$Participant == "Caregiver", var$Smoker != "NA")
rownames(var)<-var$xSampleID

prepro<-feature_table_pre_process(feature_table = data_filtered,meta_data = var,sample_var = "xSampleID", group_var = NULL,
                                  out_cut = 0.05, zero_cut = 0.95, lib_cut = 1, neg_lb = FALSE)

feature_table = prepro$feature_table 
meta_data = prepro$meta_data 
struc_zero = prepro$structure_zeros 

ancom_results<-ANCOM(feature_table = feature_table, meta_data = meta_data,
                     struc_zero = struc_zero, main_var = "Smoker",
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
  filter(var$Participant == "Caregiver", var$Smoker != "NA")
data_norm<-as.matrix(data_norm)

var<-var %>%
  filter(var$Participant == "Caregiver", var$Smoker != "NA")

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

brewer.pal(8, "Dark2")

p1<-ggplot(for_plotting, aes(x = Var2, log(value), fill = Smoker)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#66A61E", "#E6AB02"), name = "Smoker")
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
                    ncol = 2, 
                    scales = "free_x", labeller = label_wrap_gen())
p6<-p5 + theme(panel.spacing.x = unit(0.3, "lines"), 
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 10),
               legend.title = element_text(face = "bold", size = 14),
               legend.text = element_text(face = "bold", size = 12),
               legend.position = "NONE") #+

p6_parent_ancom

###CHILD COTININE###

data<-read.table("new_consented_table_for_analyses.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1,
                 stringsAsFactors = FALSE)

var<-read.table("metadata_8.24.22.txt", header = TRUE, 
                row.names = 1, sep ='\t', stringsAsFactors = FALSE)


data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Participant == "Child", var$Cotinine_pseudo > 0)

data_norm<-as.matrix(data_norm)


data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0005) 

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Participant == "Child", var$Cotinine_pseudo > 0)
data_filtered<-as.data.frame((as.matrix(data_filtered)))
var<-var %>%
  filter(var$Participant == "Child", var$Cotinine_pseudo > 0)

fit_data <- Maaslin2(input_data = data_filtered,
                     input_metadata = var,
                     output = "Child_cotinine",
                     min_prevalence = 0.0000001,
                     fixed_effects = c("Cotinine_pseudo"),
                     max_significance = 0.05)

###PLOT CHILD MAASLIN2###
data<-read.table("Child_cotinine/significant_results.tsv", header = TRUE, sep = '\t')

p1<-ggplot(data, aes(reorder(feature,coef), coef,
                     fill = coef > 0)) +
  geom_bar(stat = "identity", color = "Black", size = 1) + 
  coord_flip(ylim = c(-0.1,1)) + 
  scale_fill_manual(values = c("Darkgreen"))
p1
p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major.x = element_line(size = 0.25,
                                                 color = "black"),
               panel.grid.minor = element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               axis.ticks.y = element_blank()) +
  geom_hline(yintercept = 0:0.3, size = 1)
p3
p4<- p3 + theme(axis.text = element_text(size = 12, face = "bold",
                                         color = "black"),
                legend.position = "none",
                axis.title = element_text(size = 16, face = "bold",
                                          color = "black"))+
  labs(x = "ASV", y = "Cotinine Linear Model Coefficient")
p4
p4_child_cotinine<-p4

###PARENT COTININE###

data<-read.table("new_consented_table_for_analyses.txt", sep = '\t', quote = "", 
                 header = TRUE, row.names = 1,
                 stringsAsFactors = FALSE)

var<-read.table("metadata_8.24.22.txt", header = TRUE, 
                row.names = 1, sep ='\t', stringsAsFactors = FALSE)


data<- data[complete.cases(data), ]

column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

data_norm<-as.data.frame(data_norm)
data_norm<-data_norm %>%
  filter(var$Participant == "Caregiver", var$Cotinine_pseudo > 0)

data_norm<-as.matrix(data_norm)


data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0005) 

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Participant == "Caregiver", var$Cotinine_pseudo > 0)
data_filtered<-as.data.frame((as.matrix(data_filtered)))
var<-var %>%
  filter(var$Participant == "Caregiver", var$Cotinine_pseudo > 0)

fit_data <- Maaslin2(input_data = data_filtered,
                     input_metadata = var,
                     output = "Parent_cotinine",
                     min_prevalence = 0.0000001,
                     fixed_effects = c("Cotinine_pseudo"),
                     max_significance = 0.05)

###PLOTTING PARENT COTININE###
data<-read.table("Parent_cotinine/significant_results.tsv", header = TRUE, sep = '\t')

p1<-ggplot(data, aes(reorder(feature,coef), coef,
                     fill = coef > 0)) +
  geom_bar(stat = "identity", color = "Black", size = 1) + 
  coord_flip(ylim = c(-1,2)) + 
  scale_fill_manual(values = c("Darkred","Darkgreen"))
p1
p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major.x = element_line(size = 0.25,
                                                 color = "black"),
               panel.grid.minor = element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               axis.ticks.y = element_blank()) +
  geom_hline(yintercept = 0:0.3, size = 1)
p3
p4<- p3 + theme(axis.text = element_text(size = 12, face = "bold",
                                         color = "black"),
                legend.position = "none",
                axis.title = element_text(size = 16, face = "bold",
                                          color = "black"))+
  labs(x = "ASV", y = "Cotinine Linear Model Coefficient")
p4
p4_parent_cotinine<-p4


pcombined<-(p5_mds_child/p5_mds_parent)|(p4_child_cotinine/p6_parent_ancom/p4_parent_cotinine)
pcombined
pcombined + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", color = "black"))
