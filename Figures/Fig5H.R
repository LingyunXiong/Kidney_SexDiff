#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)

### Bar plot for TF binding + Jaspar motifs (Mouse) ###
data1 <- data.frame(Sex=rep(c("M-WT", "M-KO"), each=2),
                   Group=rep(c("AR-B", "Hnf4a-B"),2),
                   Percent_Gene=c(13.2,6.14,7.87,5.62))
data2 <- data.frame(Sex=rep(c("M-WT", "M-KO"), each=2),
                   Group=rep(c("AR-M", "Hnf4a-M"),2),
                   Percent_Gene=c(5.7,3.1,1.69,2.8))
data <- rbind(data1,data2)
data$Group <- factor(data$Group, levels = c("AR-B", "Hnf4a-B","AR-M", "Hnf4a-M"))
data$Sex <- factor(data$Sex, levels = c("M-WT", "M-KO"))

ggplot(data=data, aes(x=Group, y=Percent_Gene, color=Sex)) +
  geom_bar(stat="identity", position=position_dodge(), fill="white") + 
  ylab("Percentage of Peaks (%)") + xlab("") + ylim(0,16) +
  scale_color_manual(values=c("#3E74D1","#E22146")) + theme_bw()
 