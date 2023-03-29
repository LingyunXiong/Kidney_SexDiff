#### Load libraries ####
library(data.table)
library(openxlsx)
library(dplyr)
library(ggplot2)

### Scatter Plot ###
df_human_mouse <- fread("../Data/HumanMouse_CommonSexGenes_PT_log2fc.csv")

ggplot(df_human_mouse, aes(x=log2FC_FvsM_Mouse, y=log2FC_FvsM_Human, color=Sex, label=symbol)) +
  geom_point() + scale_color_manual(values=c("#FF8596","#7CADED")) +
  xlab("log2 FC (F vs M) -- Mouse") + ylab("log2 FC (F vs M) -- Human") +
  geom_hline(yintercept=0, linetype="dotted", color = "black") +
  geom_vline(xintercept=0, linetype="dotted", color = "black") +
  geom_smooth(method=lm,  linetype="dashed", color="gray", fill="lightgray") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
