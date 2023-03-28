#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

### For readings in individual samples ###
df_out <- fread("../Data/TMM_Normalized_Counts_13448_genes.csv")
gene2plot <- args[1]

# 79-week samples: breeder and virgin plotted together #
gene_exp <- as.numeric(filter(df_out, symbol==gene2plot)[,c(2:29,54)])
df_gene <- data.frame(Sex=rep("F",29),level=gene_exp,seq=c(rep(1,4),rep(2,5),rep(3,5),rep(4,5),rep(5,10)))
gene_exp_m <- as.numeric(filter(df_out, symbol==gene2plot)[,c(30:53,55:58)])
df_gene_m <- data.frame(Sex=rep("M",28),level=gene_exp_m,seq=c(rep(1,4),rep(2,5),rep(3,5),rep(4,5),rep(5,9)))
df_gene <- rbind(df_gene,df_gene_m)

### For readings of averaged tracing ###
df_data <- fread("../Data/TMM_Normalized_Counts_13448_genes_Mean_Delta.csv")

GE_trace <- as.numeric(filter(df_data, symbol==gene2plot & Sex=="F")[,2:6])
df_trace <- data.frame(Sex=rep("F",5),level=GE_trace,seq=c(1:5))
GE_trace_m <- as.numeric(filter(df_data, symbol==gene2plot & Sex=="M")[,2:6])
df_trace_m <- data.frame(Sex=rep("M",5),level=GE_trace_m,seq=c(1:5))
df_trace <- rbind(df_trace,df_trace_m)

### box plot and path ###
p <- ggplot(df_gene,aes(x=as.factor(seq),y=level,fill=Sex)) + 
  geom_boxplot(varwidth = TRUE, alpha=0.5) + 
  geom_path(data = df_trace, aes(x=as.factor(seq),y=level,group=Sex,color=Sex), size=2) +
  ylab(paste0("Normalized mRNA level (",gene2plot,")")) + xlab("Age (Weeks)") +
  scale_x_discrete(labels=c("1" = "0", "2" = "2", "3" = "4", "4" = "8", "5" = "79")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(file = paste0("Normalized_Exp_",gene2plot,"_BoxLine.pdf"), width = 4.6, height = 3)
p
dev.off()
