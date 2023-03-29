#### Load libraries ####
library(data.table)
library(openxlsx)
library(dplyr)
library(edgeR)
library(DESeq2)
library(biomaRt)  
library(ggplot2)

### Pre-process the input data ###
df_ct1 <- fread("../Data/Liver_ArKO_Read_Count_Table.txt"); colnames(df_ct1)[1] <- "symbol"
df_mt1 <- fread("../Data/Liver_ArKO_Meta_Table.txt")

### subset input table for specific KO experiment ###
# (4) For PCA: Ar KO (Liver) #
df_mt <- filter(df_mt1, Genotype %in% c("AlbCRE/+, Ar c/y", "Sox2Cre/+; Ar c/y", "Ar c/y", "Ar c/c"))
df_mt$Sample
df_ct <- dplyr::select(df_ct1,symbol,
                       `1_1101.Liver_AlbCRE_Arcy`,`2_1102.Liver_AlbCRE_Arcy`,`3_1103.Liver_AlbCRE_Arcy`,`4_1104.Liver_AlbCRE_Arcy`,`5_1105.Liver_AlbCRE_Arcy`,
                       `6_1106.Liver_Arcy`,`7_1107.Liver_Arcy`,`8_1108.Liver_Arcy`,`9_1109.Liver_Arcy`,`10_1110.Liver_Arcy`,
                       `1146_Sox2CRE_ARcy_M_liver`,`1147_Sox2CRE_ARcy_M_liver`,`1148_Sox2CRE_ARcy_M_liver`,`1149_Sox2CRE_ARcy_M_liver`,`1150_Sox2CRE_ARcy_M_liver`,
                       `1151_ARccF_liver`,`1152_ARccF_liver`,`1153_ARccF_liver`,`1154_ARccF_liver`)
df_ct <- filter(df_ct, symbol %in% c(vec_M,vec_F))

### edgeR processing ###
mtx_data <- as.matrix(df_ct[,c(2:20)])
rownames(mtx_data) <- df_ct$symbol

d <- DGEList(counts=mtx_data,group=factor(df_mt$Category))

cpm(d)[grep("^Ar$",df_ct$symbol),]

keep <- rowSums(cpm(d)>5) >= 1 
length(which(keep)) 
d <- d[keep,]
dim(d)

df_ct_cpm <- filter(df_ct, symbol %in% rownames(d$counts))

### DE analysis using DESeq2 ###
df_ct <- df_ct_cpm
df_data <- as.matrix(df_ct[,c(2:20)])
rownames(df_data) <- df_ct$symbol
df_mt$Category <- factor(df_mt$Category, levels = c("WT","KO"))
dds <- DESeqDataSetFromMatrix(countData = df_data, colData = df_mt, design = ~ Category) 

### PCA plot ###
### rlog transformation ###
dds <- estimateSizeFactors(dds) 
vsd <- vst(dds,blind=FALSE,nsub=nrow(dds))

# plot with ggplot2 #
pcadata <- plotPCA(vsd, c("Category","Sex"), ntop = 1144, returnData=TRUE)
vec_ttr <- df_mt$Genotype
vec_ttr[which(vec_ttr %in% c("Ar c/c", "Ar c/y"))] <- "WT"
vec_ttr[which(vec_ttr == "AlbCRE/+, Ar c/y")] <- "Alb-Ar-KO"
vec_ttr[which(vec_ttr == "Sox2Cre/+; Ar c/y")] <- "Sox2-Ar-KO"
pcadata$Treatment <- factor(vec_ttr, levels = c("WT","Alb-Ar-KO","Sox2-Ar-KO"))
percentVar <- round(100 * attr(pcadata, "percentVar"))

ggplot(arrange(pcadata, Treatment), aes(PC1, PC2, color=Treatment, shape=Sex)) + 
  geom_point(size=3, alpha=0.75) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
