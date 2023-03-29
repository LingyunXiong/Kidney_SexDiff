#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)

### TF putative targets from ChEA3 ###
data <- data.frame(Sex=rep(c("F", "M"), each=6),
                   TF=rep(c("Hnf4a","Ar","Hnf4g","Ppara","Nr1h4","Tbx10"),2),
                   value=c(86/214,39/214,37/214,14/214,28/214,19/214,
                           155/337,61/337,67/337,36/337,37/337,30/337))
data$TF <- factor(data$TF, levels = rev(c("Hnf4a","Hnf4g","Ar","Ppara","Nr1h4","Tbx10")))

ggplot(data=data, aes(x=TF, y=value*100, color=Sex, fill=Sex)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("") + xlab("") + coord_flip() + ylim(0,60) +
  scale_color_manual(values=c("#FF8596","#7CADED")) + 
  scale_fill_manual(values=c("#FF8596","#7CADED")) +
  theme_bw() + theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size =10, face= "bold"),
    axis.text.y = element_text(size =10, face= "bold"),
    panel.grid.major = element_blank(), 
