library(dplyr)
library(tidyverse)
library(rcartocolor)
library(RColorBrewer)
library(pals)
library(reshape2)
data<-read.table("top_ten_bacterial_families.txt", header = TRUE, row.names = 1,
                 sep= '\t', check.names = FALSE)

var<-read.table("metadata_8.24.22.txt", header = TRUE, row.names = 1, sep ='\t')
melted<-melt(data)
for_plotting<-cbind(melted,var)

safe_pal<-carto_pal(11,"Safe")
p1<-ggplot(for_plotting, aes(x = for_plotting$Pair_participant, value, fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") + 
  scale_fill_manual(values = safe_pal)
p2<-p1 + theme_bw()
p3<-p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p4<-p3 + labs(x = "Sample", y = "Relative Abundance", fill = "Bacterial Family") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
          axis.text.x = element_text(angle = 90, h = 1, v = .5, color = "Black", size=6),
        axis.title = element_text(face = "bold", size = 14))
p5<-p4 + facet_wrap(~for_plotting$Plot_dummy, 
                    nrow = 4, 
                    scales = "free_x")
p6<-p5 + theme(panel.spacing.x = unit(0.3, "lines"), panel.border = element_blank(),
               strip.background = element_blank(),
               strip.text = element_blank(),
               legend.title = element_text(face = "bold", size = 14),
               legend.text = element_text(face = "bold", size = 12),
               legend.position = "bottom") +
  scale_y_continuous(expand = c(0,0))
p6