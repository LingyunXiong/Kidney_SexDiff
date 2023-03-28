#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)

#### Read in MEME-SEA outputs ####
df_pt1_f <- fread("../Data/PTs1_DARs_nearSexBiasedGenes_fWTvsmWT_F.tsv")
df_pt1_m <- fread("../Data/PTs1_DARs_nearSexBiasedGenes_fWTvsmWT_M.tsv")
df_pt2_f <- fread("../Data/PTs2_DARs_nearSexBiasedGenes_fWTvsmWT_F.tsv")
df_pt2_m <- fread("../Data/PTs2_DARs_nearSexBiasedGenes_fWTvsmWT_M.tsv")
df_pt3_f <- fread("../Data/PTs3_DARs_nearSexBiasedGenes_fWTvsmWT_F.tsv")
df_pt3_m <- fread("../Data/PTs3_DARs_nearSexBiasedGenes_fWTvsmWT_M.tsv")

vec_common <- sort(intersect(df_pt3_f$ID, df_pt3_m$ID))
vec_f_only <- sort(setdiff(df_pt3_f$ID, df_pt3_m$ID))
vec_m_only <- sort(setdiff(df_pt3_m$ID, df_pt3_f$ID))

#### Subset TFs and rearrange tables ####
vec_motifs <- c("ANDR_MOUSE.H11MO.0.A",
                "RFX3_MOUSE.H11MO.0.C",
                "HNF1A_MOUSE.H11MO.0.A",
                "HNF4A_MOUSE.H11MO.0.A",
                "HNF4G_MOUSE.H11MO.0.C",
                "PPARA_MOUSE.H11MO.0.A",
                "PPARA_MOUSE.H11MO.1.A",
                "ETV5_MOUSE.H11MO.0.D",
                "NR1H4_MOUSE.H11MO.0.A",
                "HNF1B_MOUSE.H11MO.0.A",
                "STA5B_MOUSE.H11MO.0.A",
                "STA5A_MOUSE.H11MO.0.A",
                "TEAD1_MOUSE.H11MO.0.A",
                "CEBPA_MOUSE.H11MO.0.A",
                "NFIA_MOUSE.H11MO.0.C")
TF_order <- c("AR","Rfx3","Hnf1a","Hnf1b","Hnf4g","Hnf4a","Ppara","Etv5","Nr1h4","Stat5b","Stat5a","Tead1","Cebpa","Nfia")

df_pt3_fs <- filter(df_pt3_f, ID %in% vec_motifs)
df_pt3_fs <- select(df_pt3_fs,ID,`TP%`,QVALUE)
df_pt3_fs$TF <- c("Hnf1b","Hnf1a","Hnf4g","Hnf4a","Nr1h4","Stat5b","Ppara","Etv5","Stat5a","Tead1","Nfia","Cebpa")

df_pt3_ms <- filter(df_pt3_m, ID %in% vec_motifs)
df_pt3_ms <- select(df_pt3_ms,ID,`TP%`,QVALUE)
df_pt3_ms$TF <- c("Hnf1a","Ppara","Hnf4g","Etv5","Hnf4a","AR","Rfx3","Nr1h4")

df_pt3_f_fdr <- select(df_pt3_fs,TF,QVALUE); colnames(df_pt3_f_fdr)[2] <- "s3f"
df_pt3_m_fdr <- select(df_pt3_ms,TF,QVALUE); colnames(df_pt3_m_fdr)[2] <- "s3m"
df_pt3_fdr <- merge(df_pt3_f_fdr,df_pt3_m_fdr,by="TF",all=T)

df_pt3_f_prop <- select(df_pt3_fs,TF,`TP%`); colnames(df_pt3_f_prop)[2] <- "s3f"
df_pt3_m_prop <- select(df_pt3_ms,TF,`TP%`); colnames(df_pt3_m_prop)[2] <- "s3m"
df_pt3_prop <- merge(df_pt3_f_prop,df_pt3_m_prop,by="TF",all=T)

#### Other segments ####
df_pt1_fs <- filter(df_pt1_f, ID %in% vec_motifs)
df_pt1_fs <- select(df_pt1_fs,ID,`TP%`,QVALUE)
df_pt1_fs <- filter(df_pt1_fs,QVALUE<0.05)
df_pt1_fs$TF <- c("Hnf1b","Hnf4g")
df_pt1_ms <- filter(df_pt1_m, ID %in% vec_motifs)
df_pt1_ms <- select(df_pt1_ms,ID,`TP%`,QVALUE)
df_pt1_ms <- filter(df_pt1_ms,QVALUE<0.05)

df_pt1_f_fdr <- select(df_pt1_fs,TF,QVALUE); colnames(df_pt1_f_fdr)[2] <- "s1f"
#df_pt1_m_fdr <- select(df_pt1_ms,TF,QVALUE); colnames(df_pt1_m_fdr)[2] <- "s1m"
df_pt1_fdr <- df_pt1_f_fdr

df_pt1_f_prop <- select(df_pt1_fs,TF,`TP%`); colnames(df_pt1_f_prop)[2] <- "s1f"
#df_pt1_m_prop <- select(df_pt1_ms,TF,`TP%`); colnames(df_pt1_m_prop)[2] <- "s1m"
df_pt1_prop <- df_pt1_f_prop

df_pt2_fs <- filter(df_pt2_f, ID %in% vec_motifs)
df_pt2_fs <- select(df_pt2_fs,ID,`TP%`,QVALUE)
df_pt2_fs$TF <- c("Hnf1a","Hnf4g","Hnf4a","Ppara","Hnf1b","Cebpa")
df_pt2_ms <- filter(df_pt2_m, ID %in% vec_motifs)
df_pt2_ms <- select(df_pt2_ms,ID,`TP%`,QVALUE)
df_pt2_ms <- df_pt2_ms[-6,]
df_pt2_ms$TF <- c("Hnf1b","Hnf1a","Hnf4g","Hnf4a","Ppara")

df_pt2_f_fdr <- select(df_pt2_fs,TF,QVALUE); colnames(df_pt2_f_fdr)[2] <- "s2f"
df_pt2_m_fdr <- select(df_pt2_ms,TF,QVALUE); colnames(df_pt2_m_fdr)[2] <- "s2m"
df_pt2_fdr <- merge(df_pt2_f_fdr,df_pt2_m_fdr,by="TF",all=T)

df_pt2_f_prop <- select(df_pt2_fs,TF,`TP%`); colnames(df_pt2_f_prop)[2] <- "s2f"
df_pt2_m_prop <- select(df_pt2_ms,TF,`TP%`); colnames(df_pt2_m_prop)[2] <- "s2m"
df_pt2_prop <- merge(df_pt2_f_prop,df_pt2_m_prop,by="TF",all=T)

#### Prepare data for plots ####
df_pt1_fp_fdr <- select(df_pt1_fdr,TF,s1f); colnames(df_pt1_fp_fdr)[2] <- "FDR"
df_pt1_fp_prop <- select(df_pt1_prop,TF,s1f); colnames(df_pt1_fp_prop)[2] <- "Proportion"
df_pt1_fp <- merge(df_pt1_fp_prop,df_pt1_fp_fdr,by="TF"); df_pt1_fp$group <- rep("PT-S1-fWT",2)

df_pt1_mp <- df_pt1_fp
df_pt1_mp$Proportion <- rep(NA,2)
df_pt1_mp$FDR <- rep(NA,2)
df_pt1_mp$group <- rep("PT-S1-mWT",2)

df_pt2_fp_fdr <- select(df_pt2_fdr,TF,s2f); colnames(df_pt2_fp_fdr)[2] <- "FDR"
df_pt2_fp_prop <- select(df_pt2_prop,TF,s2f); colnames(df_pt2_fp_prop)[2] <- "Proportion"
df_pt2_fp <- merge(df_pt2_fp_prop,df_pt2_fp_fdr,by="TF"); df_pt2_fp$group <- rep("PT-S2-fWT",6)

df_pt2_mp_fdr <- select(df_pt2_fdr,TF,s2m); colnames(df_pt2_mp_fdr)[2] <- "FDR"
df_pt2_mp_prop <- select(df_pt2_prop,TF,s2m); colnames(df_pt2_mp_prop)[2] <- "Proportion"
df_pt2_mp <- merge(df_pt2_mp_prop,df_pt2_mp_fdr,by="TF"); df_pt2_mp$group <- rep("PT-S2-mWT",6)

df_pt3_fp_fdr <- select(df_pt3_fdr,TF,s3f); colnames(df_pt3_fp_fdr)[2] <- "FDR"
df_pt3_fp_prop <- select(df_pt3_prop,TF,s3f); colnames(df_pt3_fp_prop)[2] <- "Proportion"
df_pt3_fp <- merge(df_pt3_fp_prop,df_pt3_fp_fdr,by="TF"); df_pt3_fp$group <- rep("PT-S3-fWT",length(TF_order))

df_pt3_mp_fdr <- select(df_pt3_fdr,TF,s3m); colnames(df_pt3_mp_fdr)[2] <- "FDR"
df_pt3_mp_prop <- select(df_pt3_prop,TF,s3m); colnames(df_pt3_mp_prop)[2] <- "Proportion"
df_pt3_mp <- merge(df_pt3_mp_prop,df_pt3_mp_fdr,by="TF"); df_pt3_mp$group <- rep("PT-S3-mWT",length(TF_order))

#df_plot <- rbind(df_pt3_mp,df_pt3_fp)
df_plot <- rbind(df_pt1_fp,df_pt1_mp,df_pt2_fp,df_pt2_mp,df_pt3_mp,df_pt3_fp)
df_plot$TF <- factor(df_plot$TF,levels = rev(TF_order))
#df_plot$group <- factor(df_plot$group,levels = c("PT-S3-mWT","PT-S3-fWT"))
df_plot$group <- factor(df_plot$group,levels = c("PT-S1-mWT","PT-S2-mWT","PT-S3-mWT","PT-S1-fWT","PT-S2-fWT","PT-S3-fWT"))
df_plot$`-log10(FDR)` <- -log10(df_plot$FDR)
df_plot$Proportion <- df_plot$Proportion/100

#### dot plot ####
ggplot(data = df_plot, aes(x = group, y = TF, color = `-log10(FDR)`, size = Proportion)) + 
  geom_point() + theme_bw() + 
  scale_color_gradient(low = "grey", high = "red") +
  ylab("") + xlab("") + ggtitle("Differential Peaks") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face= "bold"))

