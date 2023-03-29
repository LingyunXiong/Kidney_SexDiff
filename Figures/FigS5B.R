#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)

### PT-DARs (12 pairwise comparison) #
data <- data.frame(Group=rep(c("F-KO vs F-WT","M-KO vs M-WT","F-WT vs M-WT","F-WT vs M-KO"), each=3),
                   PT=rep(c("PT-S1","PT-S2","PT-S3"),4),
                   No_Gene=c(2,4,3,4,189,2027,2929,12837,24399,4377,3627,9309))
data$PT <- factor(data$PT, levels = c("PT-S1","PT-S2","PT-S3"))
data$Group <- factor(data$Group, levels = rev(c("F-KO vs F-WT","M-KO vs M-WT","F-WT vs M-WT","F-WT vs M-KO")))

ggplot(data=data, aes(x=Group, y=No_Gene, fill=PT)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("No. of DARs") + xlab("") + ylim(0,25000) + coord_flip() +
  scale_fill_manual(values=c("#EEFF33","#FFB733","#FF6933")) + 
  theme_bw() + theme(
    axis.text.x = element_text(size =8, face= "bold"),
    axis.text.y = element_text(size =10, face= "bold"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())
