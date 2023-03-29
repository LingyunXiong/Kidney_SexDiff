#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)

### Compare Liver sex-specific genes across studies ###
data <- data.frame(Sex=rep(c("F", "M"), each=3),
                   Group=rep(c("In-house", "GSE112947", "GSE174535"),2),
                   No_Gene=c(981,698,600,701,446,391))
data$Group <- factor(data$Group, levels = c("In-house", "GSE112947", "GSE174535"))

ggplot(data=data, aes(x=Group, y=No_Gene, fill=Sex)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("No. of DEGs") + xlab("") +
  scale_fill_manual(values=c("#AE0428","#2657AF")) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
