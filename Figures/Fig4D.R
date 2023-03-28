#### Load libraries ####
library(data.table)
library(dplyr)
library(Signac)
library(Seurat)
library(ggplot2)
set.seed(1234)

#### Download processed single-nuclear multiomic dataset ####
# Link: https://drive.google.com/file/d/1T0BkBz8OKmRX-POintJC9c-wc3I68Bm2/view?usp=sharing

#### Load R object ####
pt <- readRDS("Integrated_ArchRpeaks_PT.rds")
DefaultAssay(pt) <- "RNA"

#### Feature Plot ####
gene.name <- "Prlr"
pt$sortby <- as.vector(pt[['RNA']][which(rownames(pt[['RNA']])==gene.name),])
new.cell.order <- rownames(pt[[]])[order(pt$sortby, decreasing = FALSE)]
pdf(paste0("plots/Integrated_ArchRpeaks_PT_FeaturePlot_", gene.name,".pdf"), width = 5, height = 5)	
FeaturePlot(pt, features = gene.name, reduction = "wnn.umap", cells = new.cell.order, cols = c("lightgrey", "darkblue")) + ggtitle(gene.name)
dev.off()

#### Violin Plot ####
pdf(paste0("plots/Integrated_ArchRpeaks_PT_ViolinPlot_", gene.name,".pdf"), width = 9, height = 5) 
VlnPlot(pt, features = gene.name, group.by = "cell_type_split", pt.size = 0)
dev.off() 

#### Dot Plot ####
my_levels_full <- c("PT-S1-mWT","PT-S1-fWT","PT-S1-fKO","PT-S1-mKO","PT-S2-mWT","PT-S2-fWT","PT-S2-fKO","PT-S2-mKO","PT-S3-mWT","PT-S3-fWT","PT-S3-fKO","PT-S3-mKO")

pt@active.ident <- factor(pt$cell_type_split, levels = my_levels_full)

vec_markers <- c("Slc7a12","Prlr","Cyp4a14","Kynu","Slc22a8","Slc22a7","Cyp51",
				 "Cyp7b1","Slc7a13","Cyp2j13","Cyp4b1","Slc22a28","Slc5a12","Slc5a2",
				 "Hnf1a","Hnf1b","Hnf4a","Hnf4g","Lrp2")

my_levels <- c("PT-S1-mWT","PT-S1-fWT","PT-S1-mKO","PT-S2-mWT","PT-S2-fWT","PT-S2-mKO","PT-S3-mWT","PT-S3-fWT","PT-S3-mKO")					
p <- DotPlot(pt, features = rev(vec_markers), cols = c("lightgrey", "red"), idents = my_levels) + RotatedAxis()

pdf("plots/Integrated_2SexPT_PTMarker_DotPlot_HierCluster_Split_re1.pdf", width = 5.4, height = 5.8)				 
p + coord_flip()
dev.off() 
