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

p5 <- DimPlot(Merged_combined, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p6 <- DimPlot(Merged_combined, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p7 <- DimPlot(Merged_combined, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

pdf("Integrated_ArchRpeaks_UMAP_WNN_LogNrpca_K20_4Group_r5.pdf", width = 14, height = 5)
p5 + p6 + p7 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()
