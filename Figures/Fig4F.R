#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)

# mKO vs. mWT (highlight single-nuclear sex-biased genes) #
data <- data.frame(Sex=rep(c("F", "M"), each=3),
                   Group=rep(c("PT-S1","PT-S2","PT-S3"),2),
                   No_Gene=c(8,84,180,8,113,153))
data$Group <- factor(data$Group, levels = c("PT-S1","PT-S2","PT-S3"))

ggplot(data=data, aes(x=Group, y=No_Gene, fill=Sex)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("No. of DEGs") + xlab("") + ylim(0,200) +
  scale_fill_manual(values=c("#FF8596","#7CADED")) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
 