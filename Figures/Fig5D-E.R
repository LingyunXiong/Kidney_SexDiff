#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(RColorBrewer)
library(gplots)

#### Read in files ####
df <- fread("../Data/Multiome_PT_fWTvsmWT_DEGs_ARresponse.csv")
vec_M <- filter(df,Direction=="M")$symbol
vec_F <- filter(df,Direction=="F")$symbol

df_z <- read.csv(file = "../Data/Kidney_snMultiomic_GAS.csv",stringsAsFactors = TRUE); colnames(df_z)[1] <- "symbol"
df_z <- select(df_z,"symbol","PT.S1.mWT","PT.S2.mWT","PT.S3.mWT","PT.S1.fWT","PT.S2.fWT","PT.S3.fWT",
               "PT.S1.mKO","PT.S2.mKO","PT.S3.mKO","PT.S1.fKO","PT.S2.fKO","PT.S3.fKO")
df_z <- filter(df_z, symbol %in% c(vec_F,vec_M))
mtx <- as.matrix(df_z[, c(2:10)])*1000000
rownames(mtx) <- df_z$symbol

# Filter genes with zero scores #
vec_rowsum <- rowSums(mtx)
mtx <- mtx[-which(is.na(vec_rowsum)),]

#vec_max <- rowMaxs(mtx)
vec_mean <- rowMeans(mtx)
mtx_z <- mtx/vec_mean
mtx_z <- mtx_z[which(rownames(mtx_z) %in% vec_M),]
#mtx_z <- mtx_z[which(rownames(mtx_z) %in% vec_F),]

#### cluster heatmaps ####
hr <- hclust(as.dist(1-cor(t(mtx_z), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height)/1.2)
length(unique(mycl))
mycolhc <- rainbow(length(unique(mycl)), start=0.2, end=0.8)
mycolhc <- mycolhc[as.vector(mycl)]
mycol <- rev(brewer.pal(n=8,name="RdGy"))

### plot heatmap (without clustering columns) ###
tiff(file="Heatmap_MultiomeSexGenes_GeneAccessScore_Norm2Mean_M_ArResponse.tiff", units="in", width=8, height=4, res=150)
heatmap.2(mtx_z, Rowv=as.dendrogram(hr), RowSideColors=mycolhc, 
          Colv=FALSE, col=mycol, scale="row",
          density.info="none", trace="none")
dev.off()

### Plot changes in accessibility ###
# fWT vs. mWT #
fc_s1_fWTvsmWT_F <- mtx_z[,4]/mtx_z[,1]
fc_s2_fWTvsmWT_F <- mtx_z[,5]/mtx_z[,2]
fc_s3_fWTvsmWT_F <- mtx_z[,6]/mtx_z[,3]
df1 <- data.frame(symbol=rownames(mtx_z),FC=fc_s1_fWTvsmWT_F,sex=rep("F",dim(mtx_z)[1]),group=rep("PT-S1",dim(mtx_z)[1]))
df2 <- data.frame(symbol=rownames(mtx_z),FC=fc_s2_fWTvsmWT_F,sex=rep("F",dim(mtx_z)[1]),group=rep("PT-S2",dim(mtx_z)[1]))
df3 <- data.frame(symbol=rownames(mtx_z),FC=fc_s3_fWTvsmWT_F,sex=rep("F",dim(mtx_z)[1]),group=rep("PT-S3",dim(mtx_z)[1]))
df_GAS_FC_F <- rbind(df1,df2,df3)

fc_s1_fWTvsmWT_M <- mtx_z[,4]/mtx_z[,1]
fc_s2_fWTvsmWT_M <- mtx_z[,5]/mtx_z[,2]
fc_s3_fWTvsmWT_M <- mtx_z[,6]/mtx_z[,3]
df1 <- data.frame(symbol=rownames(mtx_z),FC=fc_s1_fWTvsmWT_M,sex=rep("M",dim(mtx_z)[1]),group=rep("PT-S1",dim(mtx_z)[1]))
df2 <- data.frame(symbol=rownames(mtx_z),FC=fc_s2_fWTvsmWT_M,sex=rep("M",dim(mtx_z)[1]),group=rep("PT-S2",dim(mtx_z)[1]))
df3 <- data.frame(symbol=rownames(mtx_z),FC=fc_s3_fWTvsmWT_M,sex=rep("M",dim(mtx_z)[1]),group=rep("PT-S3",dim(mtx_z)[1]))
df_GAS_FC_M <- rbind(df1,df2,df3)

df_plot <- rbind(df_GAS_FC_F,df_GAS_FC_M)
write.csv(df_plot,"GAS/Multiome_GAS_FC_fWTvsmWT_ARresponse.csv",row.names = F, quote = F)

# mKO vs. mWT #
fc_s1_fWTvsmWT_F <- mtx_z[,7]/mtx_z[,1]
fc_s2_fWTvsmWT_F <- mtx_z[,8]/mtx_z[,2]
fc_s3_fWTvsmWT_F <- mtx_z[,9]/mtx_z[,3]
df1 <- data.frame(symbol=rownames(mtx_z),FC=fc_s1_fWTvsmWT_F,sex=rep("F",dim(mtx_z)[1]),group=rep("PT-S1",dim(mtx_z)[1]))
df2 <- data.frame(symbol=rownames(mtx_z),FC=fc_s2_fWTvsmWT_F,sex=rep("F",dim(mtx_z)[1]),group=rep("PT-S2",dim(mtx_z)[1]))
df3 <- data.frame(symbol=rownames(mtx_z),FC=fc_s3_fWTvsmWT_F,sex=rep("F",dim(mtx_z)[1]),group=rep("PT-S3",dim(mtx_z)[1]))
df_GAS_FC_F <- rbind(df1,df2,df3)

fc_s1_fWTvsmWT_M <- mtx_z[,7]/mtx_z[,1]
fc_s2_fWTvsmWT_M <- mtx_z[,8]/mtx_z[,2]
fc_s3_fWTvsmWT_M <- mtx_z[,9]/mtx_z[,3]
df1 <- data.frame(symbol=rownames(mtx_z),FC=fc_s1_fWTvsmWT_M,sex=rep("M",dim(mtx_z)[1]),group=rep("PT-S1",dim(mtx_z)[1]))
df2 <- data.frame(symbol=rownames(mtx_z),FC=fc_s2_fWTvsmWT_M,sex=rep("M",dim(mtx_z)[1]),group=rep("PT-S2",dim(mtx_z)[1]))
df3 <- data.frame(symbol=rownames(mtx_z),FC=fc_s3_fWTvsmWT_M,sex=rep("M",dim(mtx_z)[1]),group=rep("PT-S3",dim(mtx_z)[1]))
df_GAS_FC_M <- rbind(df1,df2,df3)

df_plot <- rbind(df_GAS_FC_F,df_GAS_FC_M)
write.csv(df_plot,"GAS/Multiome_GAS_FC_mKOvsmWT_ARresponse.csv",row.names = F, quote = F)

### Box plot ###
ggplot(df_plot, aes(x=group, y=log2(FC), color=sex)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0, linetype="dotted", color = "black") +
  scale_color_manual(values=c("#FF8596","#7CADED")) +
  xlab("") + ylab("log2 (Fold Change in Acessibility)") + 
  scale_y_continuous(breaks=seq(-3,3,1)) +
  theme_bw() + theme(
    axis.text.x = element_text(size =10, face= "bold"),
    axis.text.y = element_text(size =10, face= "bold"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

