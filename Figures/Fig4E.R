#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)

### Multiome: PT DEGs ###
data <- data.frame(Direction=rep(c("Up", "Down"), each=4),
                   Group=rep(c("F-KO vs F-WT","M-KO vs M-WT","F-WT vs M-WT","F-WT vs M-KO"),2),
                   No_Gene=c(21,454,1035,638,23,418,736,500)) ### cutoff on common-mean: 0.01
data$Group <- factor(data$Group, levels = rev(c("F-KO vs F-WT","M-KO vs M-WT","F-WT vs M-WT","F-WT vs M-KO")))

ggplot(data=data, aes(x=No_Gene, y=Group, color=Direction)) +
  geom_bar(stat="identity", position=position_dodge(), fill="white") + 
  ylab("") + xlab("No. of DEGs") + xlim(0,1200) +
  scale_color_manual(values=c("black","gray")) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
