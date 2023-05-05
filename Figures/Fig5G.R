#### Load libraries ####
library(data.table)
library(openxlsx)
library(dplyr)
library(ggplot2)

#### Read in peaks info and DARs ####
df_sn_sex <- fread("../Data/Multiome_PT_fWTvsmWT_DEGs_full.csv")
vec_sn_F <- filter(df_sn_sex,Direction=="F")$symbol
vec_sn_M <- filter(df_sn_sex,Direction=="M")$symbol

df_DARs_anno <- fread("../Data/GlobalPeaks_Anno_Statistics_PT-S3_mWTvsfWT.csv")
df_DARs_s <- filter(df_DARs_anno, nearestGene %in% c(vec_sn_M, vec_sn_F))

df_DAR_mWTvsfWT <- fread("../Data/ArchR_PTpeaks_GlobalPeaks_DARs_PT-S3_mWTvsfWT.csv")
df_DAR_mWTvsfWT_up <- filter(df_DAR_mWTvsfWT, Direction=="UP")
df_DAR_mWTvsfWT_up_s <- filter(df_DAR_mWTvsfWT_up, nearestGene %in% vec_sn_M)
df_DAR_mWTvsfWT_down <- filter(df_DAR_mWTvsfWT, Direction=="DOWN")
df_DAR_mWTvsfWT_down_s <- filter(df_DAR_mWTvsfWT_down, nearestGene %in% vec_sn_F)

df_DAR_mKOvsmWT <- fread("../Data/ArchR_PTpeaks_GlobalPeaks_DARs_PT-S3_mKOvsmWT.csv")
df_DAR_mKOvsmWT_up <- filter(df_DAR_mKOvsmWT, Direction=="UP")
df_DAR_mKOvsmWT_up_s <- filter(df_DAR_mKOvsmWT_up, nearestGene %in% vec_sn_F)
df_DAR_mKOvsmWT_down <- filter(df_DAR_mKOvsmWT, Direction=="DOWN")
df_DAR_mKOvsmWT_down_s <- filter(df_DAR_mKOvsmWT_down, nearestGene %in% vec_sn_M)

vec_peaks_mWT_down_in_mKO <- intersect(df_DAR_mWTvsfWT_up_s$peakID,df_DAR_mKOvsmWT_down_s$peakID)
vec_peaks_fWT_up_in_mKO <- intersect(df_DAR_mWTvsfWT_down_s$peakID,df_DAR_mKOvsmWT_up_s$peakID)

#### Volcano Plot (PT-S3) ####
df_0 <- select(df_DARs_s,peakID,Log2FC,FDR,Direction)
df_2d <- select(df_DAR_mKOvsmWT,peakID,Direction)
df_2e <- select(df_DAR_mWTvsfWT,peakID,Log2FC,FDR)
df_2e <- filter(df_2e, peakID %in% c(vec_peaks_mWT_down_in_mKO,vec_peaks_fWT_up_in_mKO))
df_2 <- merge(df_2e,df_2d,by="peakID",all.x=T)
df_1 <- select(df_DAR_mWTvsfWT,peakID,Log2FC,FDR,Direction)
df_1 <- filter(df_1, peakID %in% c(df_DAR_mWTvsfWT_up_s$peakID,df_DAR_mWTvsfWT_down_s$peakID))
df_1 <- filter(df_1, !peakID %in% c(vec_peaks_mWT_down_in_mKO,vec_peaks_fWT_up_in_mKO))
df_1$SexBias <- unlist(sapply(df_1$Direction, function(x){ifelse(x=="UP","M","F")}))
df_0 <- filter(df_0, FDR>0.05); df_0 <- df_0[sample(1:dim(df_0)[1],5000),]
df_1$Log2FC <- -df_1$Log2FC; df_0$Log2FC <- -df_0$Log2FC
df_2$Log2FC <- -df_2$Log2FC

ggplot() + geom_point(data=df_0,aes(x=Log2FC, y=-log10(FDR)), color="#E0E0E0") +
  geom_point(data=df_1,aes(x=Log2FC, y=-log10(FDR)), color="#EEEEEE") +
  geom_point(data=df_2,aes(x=Log2FC, y=-log10(FDR), col=Direction)) + 
  scale_color_manual(values=c("#3E74D1","#E22146")) + theme_minimal() +
  xlab("log2 FC (F-WT vs. M-WT)") + ylab("- log10(q-value)") + xlim(-7,5) + 
  geom_vline(xintercept = c(-0.25,0.25), linetype="dotted", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype="dotted", color = "black") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
