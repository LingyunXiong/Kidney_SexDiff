#### Load libraries ####
library(data.table)
library(openxlsx)
library(dplyr)
library(DESeq2)
library(ggplot2)

### read in the raw read table ###
df_ct <- fread("../Data/WashU_Treatment_raw_read_counts.txt")
df_mt <- fread("../Data/WashU_Treatment_meta.txt")

### filter genes ###
df_sex_core <- read.xlsx("../Data/SI_Table1.xlsx",sheet = "Table S1.2")
df_ct <- filter(df_ct2, symbol %in% df_sex_core$symbol)

### DE analysis using DESeq2 ###
df_data <- as.matrix(df_ct[,c(2:49)]) # WashU samples (150PE)
rownames(df_data) <- df_ct$symbol
dds <- DESeqDataSetFromMatrix(countData = df_data, colData = df_mt, design = ~ Sex) ### choice for PCA

### PCA plot ###
### rlog transformation ###
dds <- estimateSizeFactors(dds) 
vsd <- vst(dds,blind=FALSE,nsub=nrow(dds))

# plot with ggplot2 #
pcadata <- plotPCA(vsd, c("Sex", "gene"), ntop = 551, returnData=TRUE)
vec_ttr <- pcadata$gene
vec_ttr[which(vec_ttr %in% c("Ar","Esr1", "WT_CN"))] <- "WT"
vec_ttr[which(vec_ttr=="CM_TES")] <- "CM+TES"
vec_ttr[which(vec_ttr=="OF_TES")] <- "OF+TES"
vec_ttr[which(vec_ttr=="Six2_Ar")] <- "Ar KO"
vec_ttr[which(vec_ttr=="Six2_Esr1")] <- "ERa KO"
pcadata$Treatment <- factor(vec_ttr, levels = c("WT","CM","CM+TES","Ar KO","OF","OF+TES","ERa KO"))
percentVar <- round(100 * attr(pcadata, "percentVar"))

ggplot(arrange(pcadata, Treatment), aes(PC1, PC2, color=Treatment, shape=Sex)) + 
  geom_point(size=3, alpha=0.75) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c("#A0A0A0","#00BFC4","#99FFFF","#7CAE00","#F8766D","#FF0000","#C77CFF")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
