library(Maaslin2)
library(dplyr)
library(ggplot2)
library(patchwork)

###Adiponectin###
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
  filter(var$Adiponectin_pseudo != "NA")

data_norm<-as.matrix(data_norm)


data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0005) 

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Adiponectin_pseudo != "NA")
data_filtered<-as.data.frame((as.matrix(data_filtered)))
var<-var %>%
  filter(var$Adiponectin_pseudo != "NA")

fit_data <- Maaslin2(input_data = data_filtered,
                     input_metadata = var,
                     output = "Adiponectin",
                     min_prevalence = 0.0000001,
                     fixed_effects = c("Adiponectin_pseudo"),
                     max_significance = 0.05)

###CRP###

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
  filter(var$CRP_pseudo != "NA")

data_norm<-as.matrix(data_norm)


data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0005) 

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$CRP_pseudo != "NA")
data_filtered<-as.data.frame((as.matrix(data_filtered)))
var<-var %>%
  filter(var$CRP_pseudo != "NA")

fit_data <- Maaslin2(input_data = data_filtered,
                     input_metadata = var,
                     output = "CRP",
                     min_prevalence = 0.0000001,
                     fixed_effects = c("CRP_pseudo"),
                     max_significance = 0.05)


###URIC ACID###
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
  filter(var$Uric_acid_pseudo != "NA")

data_norm<-as.matrix(data_norm)


data_norm<-t(data_norm)
data_norm<-as.data.frame(data_norm)
data_norm$mean<-rowMeans(data_norm)

data_filtered<-data %>%
  filter(data_norm$mean > 0.0005) 

data_filtered<-as.data.frame(t(as.matrix(data_filtered)))

data_filtered<-data_filtered %>%
  filter(var$Uric_acid_pseudo != "NA")
data_filtered<-as.data.frame((as.matrix(data_filtered)))
var<-var %>%
  filter(var$Uric_acid_pseudo != "NA")

fit_data <- Maaslin2(input_data = data_filtered,
                     input_metadata = var,
                     output = "Uric_acid",
                     min_prevalence = 0.0000001,
                     fixed_effects = c("Uric_acid_pseudo"),
                     max_significance = 0.05)

###MAASLIN PLOTS###
###ADIPONECTIN###
data<-read.table("Adiponectin/significant_results.tsv", header = TRUE, sep = '\t')

p1<-ggplot(data, aes(reorder(feature,coef), coef,
                     fill = coef > 0)) +
  geom_bar(stat = "identity", color = "Black", size = 1) + 
  coord_flip(ylim = c(-1,1)) + 
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
  labs(x = "ASV", y = "Adiponectin Linear Model Coefficient")
p4
p4_adiponectin<-p4

###URIC ACID###
data<-read.table("Uric_acid/significant_results.tsv", header = TRUE, sep = '\t')

p1<-ggplot(data, aes(reorder(feature,coef), coef,
                     fill = coef > 0)) +
  geom_bar(stat = "identity", color = "Black", size = 1) + 
  coord_flip(ylim = c(-1.2,0.1)) + 
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
  labs(x = "ASV", y = "Uric Acid Linear Model Coefficient")
p4
p4_uric<-p4

###DBRDA###
library(vegan)
library(ggplot2)
library(dplyr)

table<-read.table("new_rarefied_table_for_analyses.txt", header = TRUE, sep = '\t', row.names = 1)
table<-as.data.frame(t(as.matrix(table)))

var<-read.table("metadata_8.24.22.txt", sep = '\t', header = TRUE, row.names = 1)

table<-table %>%
  filter(var$Participant == "Child", var$Uric_acid_pseudo > 0, var$CRP_pseudo > 0,
         var$Adiponectin_pseudo > 0, var$Cotinine_pseudo > 0)

var<-var %>%
  filter(var$Participant == "Child", var$Uric_acid_pseudo > 0, var$CRP_pseudo > 0,
         var$Adiponectin_pseudo > 0, var$Cotinine_pseudo > 0)
colnames(var)[5]<-"Uric.acid"
colnames(var)[7]<-"CRP"
colnames(var)[9]<-"Adiponectin"
colnames(var)[14]<-"Cotinine"
rda<-dbrda(table ~ Uric.acid	+
             CRP +
             Adiponectin +
             Cotinine, 
           var, distance = "bray", metaMDS = TRUE, sqrt.dist = FALSE)

plot(rda)
anova(rda, permu = 999)
anova(rda, permu = 999, by="term")
summary(rda)

newlab<-list('Adiponectin' = "Adiponectin",
             'Cotinine' = "Cotinine",
             'CRP' = "CRP",
             'Uric.acid' = "Uric acid")

library(ggord)
library(ggplotify)
library(patchwork)

p1<-ggord(rda, size = 1.5, color = "black", 
          vec_ext = 0.8, veclsz = 0.7, repel = TRUE, txt = 5)
p1
p2<-p1 + theme_bw() + theme(panel.grid = element_blank())
p3_dbrda<-p2 + theme(axis.title = element_text(face = "bold", color = "black",
                                               size = 12),
                     axis.text = element_text(face = "bold", color = "black",
                                              size = 10))
p3_dbrda
p3_dbrda<-as.ggplot(p3_dbrda)
pcombined<-p3_dbrda + (p4_adiponectin/p4_uric)
pcombined
pcombined + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", color = "black"))
