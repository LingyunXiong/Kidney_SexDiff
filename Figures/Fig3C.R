#### Load libraries ####
library(data.table)
library(openxlsx)
library(dplyr)
library(matrixStats)
library(RColorBrewer)
library(pheatmap)
library(gplots)

### read in gene expression matrix ###
df_z <- read.csv(file = "../Data/SexGenes_Core_Treatment_WashU_150PE_z-score.csv",stringsAsFactors = TRUE)
df_list <- read.xlsx("../Data/SI_Table1.xlsx",sheet = "Table S1.2")
#vec_filtered <- filter(df_list,Direction=="F")$symbol # for female-biased genes 
vec_filtered <- filter(df_list,Direction=="M")$symbol # for male-biased genes
df_z <- filter(df_z, symbol %in% vec_filtered)

### Hierarchical clustering ###
mtx_z <- as.matrix(df_z[, c(2:11)]) #  for all time points
rownames(mtx_z) <- df_z$symbol
hr <- hclust(as.dist(1-cor(t(mtx_z), method="pearson")), method="complete")

### cluster heatmaps ###
mycl <- cutree(hr, h=max(hr$height)/1.2)
length(unique(mycl))
df_z$heatmap_cl <- mycl
mycolhc <- rainbow(length(unique(mycl)), start=0.2, end=0.8)
mycolhc <- mycolhc[as.vector(mycl)]
mycol <- rev(brewer.pal(n=8,name="RdGy"))

### plot heatmap ###
tiff(file="Heatmap_SexGene_Treatment_Male.tiff", units="in", width=5, height=5, res=150)
heatmap.2(mtx_z, Rowv=as.dendrogram(hr), RowSideColors=mycolhc, 
          Colv=FALSE, col=mycol, scale="row", 
          density.info="none", trace="none")
dev.off()