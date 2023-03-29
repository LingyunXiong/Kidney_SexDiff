#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)

### Bar Plot ###
data <- data.frame(Sex=rep(c("F", "M"), each=4),
                   Group=rep(c("PT-S1","PT-S1/2","PT-S2","PT-S2/3"),2),
                   No_Gene=c(31,90,44,70,10,87,80,92))
data$Group <- factor(data$Group, levels = c("PT-S1","PT-S1/2","PT-S2","PT-S2/3"))

ggplot(data=data, aes(x=Group, y=No_Gene, fill=Sex)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("No. of DEGs") + xlab("") + 
  scale_fill_manual(values=c("#FF8596","#7CADED")) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
