#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)

# Human Kidney #
df_meta <- fread("../Data/Integrated_Human_UMAP_PFnorm_rPCA_meta.csv")
tbl <- table(df_meta$cell_type,df_meta$sex)
data <- as.data.frame(tbl)
colnames(data)[1:2] <- c("CellType","Group")
data$CellType <- as.character(data$CellType)
data$CellType <- factor(data$CellType,levels = c("Pod","PEC","PT-VCAM1","PT-S1","PT-S1/2","PT-S2","PT-S2/3",
                                                 "LOH","TAL","DCT","CNT","PC","IC-A","IC-B","Vascular","Interstitial","Immune"))

ggplot(data, aes(x = CellType, y = Freq, fill = Group)) + 
  geom_col(position = "fill", alpha=0.8) +
  ylab("Proportion") + xlab("") +
  scale_fill_manual(values = c("#FF8596","#7CADED")) +
  theme_bw(base_size = 10) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
