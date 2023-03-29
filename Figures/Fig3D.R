#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)

### Treatment DEGs and core sex genes ###
# Male #
data <- data.frame(Sex=rep(c("F", "M"), each=4),
                   Treatment=rep(c("CM","CM+Tes","Ar KO","ERa KO"),2),
                   value=c(77.6,39.3,55.6,0.5,85.5,39.2,82.5,NA))
data$Treatment <- factor(data$Treatment, levels = rev(c("CM","CM+Tes","Ar KO","ERa KO")))

data <- data.frame(Sex=rep(c("F", "M"), each=4),
                   Treatment=rep(c("CM","CM+Tes","Six2-Ar-KO","Sox2-Ar-KO"),2),
                   value=c(77.6,39.3,55.6,57.9,85.5,39.2,81.6,82.2))
data$Treatment <- factor(data$Treatment, levels = rev(c("CM","CM+Tes","Six2-Ar-KO","Sox2-Ar-KO")))

# Female #
data <- data.frame(Sex=rep(c("F", "M"), each=4),
                   Treatment=rep(c("OF","OF+Tes","Ar KO","ERa KO"),2),
                   value=c(7.5,26.6,NA,1.4,3.6,24,NA,NA))
data$Treatment <- factor(data$Treatment, levels = rev(c("OF","OF+Tes","Ar KO","ERa KO")))

ggplot(data=data, aes(x=Treatment, y=value, color=Sex, fill=Sex)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("Percentage of the core set perturbed (%)") + xlab("") + coord_flip() +
  scale_color_manual(values=c("#FF8596","#7CADED")) + 
  scale_fill_manual(values=c("#FF8596","#7CADED")) +
  ylim(c(0,100)) +
  theme_bw() + theme(
    axis.text.x = element_text(size =10, face= "bold"),
    axis.text.y = element_text(size =10, face= "bold"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())