### Load libraries ###
library(data.table)
library(dplyr)
library(Signac)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
set.seed(1234)

library(RColorBrewer)
my_cols = brewer.pal(4,"Spectral")

#### Download processed single-nuclear multiomic dataset ####
# Link: https://drive.google.com/file/d/1ETHKM4UduP9mjNzKS6pCV6mVzKWurC_R/view?usp=sharing

### Load merged object ###
Merged_combined <- readRDS("Integrated_ArchRpeaks_LogNorm_rpca_K20_WNN_Typed.rds")
DefaultAssay(Merged_combined) <- "RNA"
Merged_combined$cell_ID <- rownames(Merged_combined[[]])

### UMAP plots ###
new_cols <- c("#E68613","#ABA300","#0CB702","#CD9600","#00BFC4","#F8766D","#00A9FF","#FF61CC",
				"#00B8E7","#00C19A","#ED68ED","#7CAE00","#8494FF","#C77CFF","#FF68A1","#00BE67")

new.cell.order <- sample(rownames(Merged_combined[[]]),dim(Merged_combined[[]])[1])
p1 <- DimPlot(Merged_combined, reduction = "umap.rna", label = TRUE, label.size = 4, repel = TRUE, cols=alpha(new_cols, 0.8), cells = new.cell.order) + ggtitle('RNA') & NoLegend()
p2 <- DimPlot(Merged_combined, reduction = "umap.rna", group.by = 'group', cols=alpha(my_cols, 0.5), cells = new.cell.order) + ggtitle('Color by sex-genotype')
pdf("plots/Integrated_ArchRpeaks_UMAP_RNA_LogNrpca_K20_4Group_Typed_4Color.pdf", width = 11, height = 5)
p1 + p2
dev.off()  

p11 <- DimPlot(Merged_combined, reduction = "wnn.umap", label = TRUE, label.size = 4, repel = TRUE, cols=alpha(new_cols, 0.8), cells = new.cell.order) + ggtitle('WNN') & NoLegend()
p8 <- DimPlot(Merged_combined, reduction = "wnn.umap", group.by = 'group', cols=alpha(my_cols, 0.5), cells = new.cell.order) + ggtitle('Color by sex-genotype')
pdf("plots/Integrated_ArchRpeaks_UMAP_WNN_LogNrpca_K20_4Group_Typed_4Color.pdf", width = 11, height = 5)
p11 + p8
dev.off()

### Bar plot: stacked + percent (Multiome composition) ###
df_meta <- fread("~/Documents/USC/P02-Kidney/Sex_Difference/Multiome/SeuratSignac/outputs/DEGs/Integrated_ArchRpeaks_WNNmeta_recode.csv")
df_meta$group <- paste(df_meta$sex,df_meta$type,sep="-")

df_meta_s <- select(df_meta,cell_ID,dataset,sex,type,group,cell_type)
tbl <- table(df_meta_s$cell_type,df_meta_s$group)
data <- as.data.frame(tbl)
colnames(data)[1:2] <- c("CellType","Group")
data$CellType <- as.character(data$CellType)
data$CellType <- factor(data$CellType,levels = c("Pod","PEC","PT-S1","PT-S2-mWT","PT-S3-mWT","PT-S2-f/mKO","PT-S3-f/mKO",
                                               "LOH","TAL","DCT","PC","IC-A","IC-B","Vascular","Interstitial","Immune"))

ggplot(data, aes(x = CellType, y = Freq, fill = Group)) + 
  geom_col(position = "fill", alpha=0.5) +
  ylab("Proportion") + xlab("") +
  scale_fill_manual(values = brewer.pal(4,"Spectral")) +
  theme_bw(base_size = 10) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#### Dot plot (Fig.S4) ####
vec_markers <- c("Nphs2","Ncam1","Hnf4a","Lrp2","Slc5a12","Slc22a6","Slc22a7","Slc7a12","Slc7a13",
				 "Aqp1","Sptssb","Slc12a1","Slc12a3","Calb1","Hsd11b2","Aqp4","Aqp2",
				 "Atp6v1g3","Slc4a1","Slc26a4","Kdr","Flt1","Pdgfrb","Ptprc") 
				 
p <- DotPlot(Merged_combined, features = rev(vec_markers), cols = c("lightgrey", "red")) + RotatedAxis() 
pdf("plots/Integrated_ArchRpeaks_CellType_Marker_DotPlot_HierCluster_WNNr4_ordered_Typed.pdf",width = 9, height = 5)				 
p + coord_flip()
dev.off() 