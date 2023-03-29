### Load libraries ###
library(data.table)
library(dplyr)
library(Signac)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
set.seed(1234)

#### Download processed single-nuclear RNA-seq dataset ####
# Link: https://drive.google.com/file/d/1Vc34v1H11Ldk096OqMFmpWShSCqMOjxC/view?usp=sharing

### Load merged object ###
Merged_combined <- readRDS("Integrated_Human_PFnorm_rPCA_Typed.rds")

### UMAP plots ###
my_cols = c("#FF8596","#7CADED")

p1 <- DimPlot(Merged_combined, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.01) + ggtitle('') & NoLegend()
new.cell.order <- sample(rownames(Merged_combined[[]]),dim(Merged_combined[[]])[1])
p2 <- DimPlot(Merged_combined, reduction = "umap", group.by = 'sex', cols=my_cols, cells = new.cell.order, , pt.size = 0.01) + ggtitle('')

pdf("Integrated_Human_UMAP_PFnorm_rPCA_K20_r5_CellType.pdf", width = 10, height = 5)
p1 + p2
dev.off() 
