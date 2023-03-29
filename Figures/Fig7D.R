#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)

### Liver Ar KO DEGs ###
data <- data.frame(Sex=rep(c("F", "M"), each=4),
                   Treatment=rep(c("Alb-Ar-KO","Sox2-Ar-KO","GSE112947-CM","GSE112947-OF"),2),
                   value=c(15,35,32,0.5,16,34,42,0.4))
data$Treatment <- factor(data$Treatment, levels = rev(c("Alb-Ar-KO","Sox2-Ar-KO","GSE112947-CM","GSE112947-OF")))

ggplot(data=data, aes(x=Treatment, y=value, color=Sex, fill=Sex)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("Percentage of the genes perturbed (%)") + xlab("") + coord_flip() +
  scale_color_manual(values=c("#AE0428","#2657AF")) + 
  scale_fill_manual(values=c("#AE0428","#2657AF")) +
  ylim(c(0,100)) +
  theme_bw() + theme(
    axis.text.x = element_text(size =10, face= "bold"),
    axis.text.y = element_text(size =10, face= "bold"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())
