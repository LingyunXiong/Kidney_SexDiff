#### Load libraries ####
library(data.table)
library(openxlsx)
library(dplyr)
library(ggplot2)

### Read in the DEG lists ###
df_m <- fread("Bulk_snMultiomic_Six2ArKO_consensus_SexGenes.csv")

# plot #
ggplot(df_m, aes(x=ArKOsn_log2FC, y=ArKObulk_log2FC, color=Sex)) +
  geom_point() + scale_color_manual(values=c("#FF8596","#7CADED")) +
  xlab("log2 FC (sn)") + ylab("log2 FC (bulk)") +
  geom_hline(yintercept=0, linetype="dotted", color = "black") +
  geom_vline(xintercept=0, linetype="dotted", color = "black") +
  geom_smooth(method=lm,  linetype="dashed", color="gray", fill="lightgray") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### test correlations ###
cor(df_m$ArKOsn_log2FC, df_m$ArKObulk_log2FC, method = "spearman") # rho=0.80 
