#### Load libraries ####
library(data.table)
library(openxlsx)
library(dplyr)
library(ggplot2)

#### Read in PT DEG lists ####
#### F-WT vs. M-WT ####
df_pt1_fWTvsmWT <- fread("../Data/ArchRpeaks_WNN_PT-S1_fWTvsmWT_DEGs_DESeq2Norm_SSeq.csv") 
df_pt1_fWTvsmWT <- filter(df_pt1_fWTvsmWT, significant=="TRUE") # 445 genes; 158 up, 287 down
vec_pt1_fWTvsmWT_up <- filter(df_pt1_fWTvsmWT, Direction=="UP")$symbol 
vec_pt1_fWTvsmWT_down <- filter(df_pt1_fWTvsmWT, Direction=="DOWN")$symbol 

df_pt2_fWTvsmWT <- fread("../Data/ArchRpeaks_WNN_PT-S2_fWTvsmWT_DEGs_DESeq2Norm_SSeq.csv") 
df_pt2_fWTvsmWT <- filter(df_pt2_fWTvsmWT, significant=="TRUE") # 785 genes; 385 up, 400 down
vec_pt2_fWTvsmWT_up <- filter(df_pt2_fWTvsmWT, Direction=="UP")$symbol 
vec_pt2_fWTvsmWT_down <- filter(df_pt2_fWTvsmWT, Direction=="DOWN")$symbol 

df_pt3_fWTvsmWT <- fread("../Data/ArchRpeaks_WNN_PT-S3_fWTvsmWT_DEGs_DESeq2Norm_SSeq.csv") 
df_pt3_fWTvsmWT <- filter(df_pt3_fWTvsmWT, significant=="TRUE") # 1171 genes; 739 up, 432 down
vec_pt3_fWTvsmWT_up <- filter(df_pt3_fWTvsmWT, Direction=="UP")$symbol 
vec_pt3_fWTvsmWT_down <- filter(df_pt3_fWTvsmWT, Direction=="DOWN")$symbol 

vec_fWTvsmWT_up <- unique(c(vec_pt1_fWTvsmWT_up,vec_pt2_fWTvsmWT_up,vec_pt3_fWTvsmWT_up)) 
vec_fWTvsmWT_down <- unique(c(vec_pt1_fWTvsmWT_down,vec_pt2_fWTvsmWT_down,vec_pt3_fWTvsmWT_down)) 
vec_fWTvsmWT_deg <- unique(c(vec_fWTvsmWT_up,vec_fWTvsmWT_down)) 

#### M-KO VS. M-WT ####
df_pt1_mKOvsmWT <- fread("../Data/ArchRpeaks_WNN_PT-S1_mKOvsmWT_DEGs_SSeq.csv") 
df_pt1_mKOvsmWT <- filter(df_pt1_mKOvsmWT, significant=="TRUE") # 42 genes; 25 up, 17 down 
vec_pt1_mKOvsmWT_up <- filter(df_pt1_mKOvsmWT, Direction=="UP")$symbol 
vec_pt1_mKOvsmWT_down <- filter(df_pt1_mKOvsmWT, Direction=="DOWN")$symbol 

df_pt2_mKOvsmWT <- fread("../Data/ArchRpeaks_WNN_PT-S2_mKOvsmWT_DEGs_SSeq.csv") 
df_pt2_mKOvsmWT <- filter(df_pt2_mKOvsmWT, significant=="TRUE") # 303 genes; 132 up, 171 down 
vec_pt2_mKOvsmWT_up <- filter(df_pt2_mKOvsmWT, Direction=="UP")$symbol 
vec_pt2_mKOvsmWT_down <- filter(df_pt2_mKOvsmWT, Direction=="DOWN")$symbol 

df_pt3_mKOvsmWT <- fread("../Data/ArchRpeaks_WNN_PT-S3_mKOvsmWT_DEGs_SSeq.csv") 
df_pt3_mKOvsmWT <- filter(df_pt3_mKOvsmWT, significant=="TRUE") # 688 genes; 376 up, 312 down 
vec_pt3_mKOvsmWT_up <- filter(df_pt3_mKOvsmWT, Direction=="UP")$symbol 
vec_pt3_mKOvsmWT_down <- filter(df_pt3_mKOvsmWT, Direction=="DOWN")$symbol 

vec_mKOvsmWT_up <- unique(c(vec_pt1_mKOvsmWT_up,vec_pt2_mKOvsmWT_up,vec_pt3_mKOvsmWT_up)) 
vec_mKOvsmWT_down <- unique(c(vec_pt1_mKOvsmWT_down,vec_pt2_mKOvsmWT_down,vec_pt3_mKOvsmWT_down)) 
vec_mKOvsmWT_deg <- unique(c(vec_mKOvsmWT_up,vec_mKOvsmWT_down)) 

#### Identify common gene sets ####
vec_up_cc <- intersect(vec_fWTvsmWT_up,vec_mKOvsmWT_up) 
vec_down_cc <- intersect(vec_fWTvsmWT_down,vec_mKOvsmWT_down) 

#### Volcano Plot to show mKO effect among fWT-mWT: Fig.4G ####
df_plot <- select(df_pt3_fWTvsmWT,symbol,log2FC,p_adj,Direction) 
#df_plot <- select(df_pt2_fWTvsmWT,symbol,log2FC,p_adj,Direction) # Fig.S4F
#df_plot <- select(df_pt1_fWTvsmWT,symbol,log2FC,p_adj,Direction) # Fig.S4F
df_plot <- filter(df_plot, !symbol %in% c(grep("Rik",df_plot$symbol,value = T), grep("^Gm",df_plot$symbol,value = T)))
df_1 <- filter(df_plot, !symbol %in% c(vec_up_cc,vec_down_cc))
df_2 <- filter(df_plot, symbol %in% c(vec_up_cc,vec_down_cc))
df_0 <- filter(df_pt2_fWTvsmWT_full, significant=="FALSE")
df_0 <- select(df_0,symbol,log2FC,p_adj,Direction)
df_0 <- df_0[sample(1:dim(df_0)[1],5000),]

ggplot() + geom_point(data=df_0,aes(x=log2FC, y=-log10(p_adj)), color="#E0E0E0") +
  geom_point(data=df_1,aes(x=log2FC, y=-log10(p_adj)), color="lightgrey") +
  geom_point(data=df_2,aes(x=log2FC, y=-log10(p_adj), col=Direction)) + 
  scale_color_manual(values=c("#3E74D1","#E22146")) + theme_minimal() +
  xlab("log2 FC (F-WT vs. M-WT)") + ylab("- log10(FDR)") + 
  geom_vline(xintercept = c(-0.25,0.25), linetype="dotted", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype="dotted", color = "black") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#### Scatter Plot of fold changes (PT-S3): Fig.4H ####
df_m1 <- filter(select(df_pt3_fWTvsmWT,symbol,log2FC,Direction), symbol %in% c(vec_up_cc,vec_down_cc))
colnames(df_m1)[2] <- "log2FC_fWT_mWT"
df_m2 <- filter(select(df_pt3_mKOvsmWT,symbol,log2FC), symbol %in% c(vec_up_cc,vec_down_cc))
colnames(df_m2)[2] <- "log2FC_mKO_mWT"
df_m <- merge(df_m1,df_m2,by="symbol")
df_m <- df_m[-grep("Rik",df_m$symbol),]

ggplot(df_m, aes(x=log2FC_fWT_mWT, y=log2FC_mKO_mWT, color=sex, label=symbol)) +
  geom_point() + scale_color_manual(values=c("#FF8596","#7CADED")) +
  xlab("log2 FC (F-WT vs. M-WT)") + ylab("log2 FC (M-KO vs. M-WT)") +
  geom_hline(yintercept=0, linetype="dotted", color = "black") +
  geom_vline(xintercept=0, linetype="dotted", color = "black") +
  geom_smooth(method=lm,  linetype="dashed", color="gray", fill="lightgray") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
