#### Load libraries ####
library(data.table)
library(openxlsx)
library(dplyr)
library(ggplot2)

### Read input file ###
df_m <- fread("../Data/Liver_Alb_Sox2_ArKO_Male_Log2FC.csv")
df_sig_Sox2 <- read.xlsx("../Data/SI_Table4.xlsx",sheet = "Table S4.2")
df_sig_Alb <- read.xlsx("../Data/SI_Table4.xlsx",sheet = "Table S4.3")

### scatter plot of log2FC ###
ggplot(df_m, aes(x=`Sox2-Ar-KO`, y=`Alb-Ar-KO`, color=Sex, label=symbol)) +
  geom_point() + scale_color_manual(values=c("#AE0428","#2657AF")) +
  xlab("log2 FC (Sox2-Ar-KO)") + ylab("log2 FC (Alb-Ar-KO)") + 
  geom_hline(yintercept=0, linetype="dotted", color = "black") +
  geom_vline(xintercept=0, linetype="dotted", color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", color = "gray") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### only highlight perturbed liver sex genes ###
vec_sig <- intersect(df_sig_Alb$symbol, df_sig_Sox2$symbol)
df_m_1 <- filter(df_m, !symbol %in% vec_sig)
df_m_2 <- filter(df_m, symbol %in% vec_sig)

ggplot() + geom_point(data=df_m_1,aes(x=`Sox2-Ar-KO`, y=`Alb-Ar-KO`), color="#F0F0F0") +
  geom_point(data=df_m_2,aes(x=`Sox2-Ar-KO`, y=`Alb-Ar-KO`, color=Sex)) + 
  scale_color_manual(values=c("#AE0428","#2657AF")) +
  xlab("log2 FC (Sox2-Ar-KO)") + ylab("log2 FC (Alb-Ar-KO)") + 
  geom_hline(yintercept=0, linetype="dotted", color = "black") +
  geom_vline(xintercept=0, linetype="dotted", color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", color = "gray") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
