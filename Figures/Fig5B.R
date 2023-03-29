#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)

### PT-DARs #
data <- data.frame(Sex=rep(c("F", "M"), each=3),
                   Group=rep(c("PT-S1","PT-S2","PT-S3"),2),
                   No_Gene=c(1357,1776,8765,1531,7166,13493))
data$Group <- factor(data$Group, levels = c("PT-S1","PT-S2","PT-S3"))

ggplot(data=data, aes(x=Group, y=No_Gene, fill=Sex)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("No. of DARs") + xlab("") + ylim(0,15000) +
  scale_fill_manual(values=c("#FF8596","#7CADED")) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### DARs near sex-biased genes #
data <- data.frame(Sex=rep(c("F", "M"), each=3),
                   Group=rep(c("PT-S1","PT-S2","PT-S3"),2),
                   No_Gene=c(94,207,1264,74,578,1164))
data$Group <- factor(data$Group, levels = c("PT-S1","PT-S2","PT-S3"))

ggplot(data=data, aes(x=Group, y=No_Gene, fill=Sex)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("No. of DARs") + xlab("") + ylim(0,1500) +
  scale_fill_manual(values=c("#FF8596","#7CADED")) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
