#### Load libraries ####
library(data.table)
library(ggplot2)

### scatter plot of log2FC ###
df_m <- fread("Output_Data/Kidney_Male_Ar_KO_Log2FC_Six2_Sox2_FvsM.txt")

ggplot(df_m, aes(x=`Six2-Ar-KO`, y=`Sox2-Ar-KO`, color=Sex, label=symbol)) +
  geom_point() + scale_color_manual(values=c("#FF8596","#7CADED")) +
  xlab("log2 FC (Six2-Ar-KO)") +ylab("log2 FC (Sox2-Ar-KO)") + 
  geom_text(data=subset(df_m, (`Six2-Ar-KO` > 2.8 | `Six2-Ar-KO` < -2.8 | `Sox2-Ar-KO` < -5)),
            nudge_x = 0, nudge_y = 0.25, check_overlap = T, size=3) +
  geom_hline(yintercept=0, linetype="dotted", color = "black") +
  geom_vline(xintercept=0, linetype="dotted", color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", color = "black") +
  geom_smooth(method=lm,  linetype="dashed", color="gray", fill="lightgray") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
