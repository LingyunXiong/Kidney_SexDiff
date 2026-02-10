#### Load libraries ####
library(data.table)
library(dplyr)
library(Signac)
library(Seurat)
library(ggplot2)
set.seed(1234)

#### Download processed single-nuclear multiomic dataset ####
# Link: https://ucla.box.com/s/6f43jpkv5sco1bmcv8ttmp19edz9pnaw

#### Load R object ####
pt <- readRDS("Integrated_ArchRpeaks_PT.rds")

### Coverage Plot ###
DefaultAssay(pt) <- "PTpeaks"
my_levels <- c("PT-S1-mKO","PT-S2-mKO","PT-S3-mKO","PT-S1-fKO","PT-S2-fKO","PT-S3-fKO",
			   "PT-S1-mWT","PT-S2-mWT","PT-S3-mWT","PT-S1-fWT","PT-S2-fWT","PT-S3-fWT") 
pt@active.ident <- factor(pt$cell_type_split, levels = my_levels)

# For Forward-Stranded Genes #
gene.name <- "Atp11a"
gene.name <- "Gsta4"

pdf(paste0("plots/Integrated_2Sex_ArchRPTpeaks_Coverage_", gene.name,"_PT-S3.pdf"), width = 6, height = 3)	
CoveragePlot(pt, region = gene.name, features = gene.name, assay = 'PTpeaks', expression.assay = 'RNA', peaks = TRUE, 
			extend.upstream = 5000, extend.downstream = 0, scale.factor = 1e7, idents=c("PT-S3-mWT","PT-S3-fWT","PT-S3-mKO"))
dev.off()
pdf(paste0("plots/Integrated_2Sex_ArchRPTpeaks_Coverage_", gene.name,"_PT-S2.pdf"), width = 6, height = 3)	
CoveragePlot(pt, region = gene.name, features = gene.name, assay = 'PTpeaks', expression.assay = 'RNA', peaks = TRUE, 
			extend.upstream = 5000, extend.downstream = 0, scale.factor = 1e7, idents=c("PT-S2-mWT","PT-S2-fWT","PT-S2-mKO"))
dev.off()
pdf(paste0("plots/Integrated_2Sex_ArchRPTpeaks_Coverage_", gene.name,"_PT-S1.pdf"), width = 6, height = 3)	
CoveragePlot(pt, region = gene.name, features = gene.name, assay = 'PTpeaks', expression.assay = 'RNA', peaks = TRUE, 
			extend.upstream = 5000, extend.downstream = 0, scale.factor = 1e7, idents=c("PT-S1-mWT","PT-S1-fWT","PT-S1-mKO"))
dev.off()

# For Backward-Stranded Genes #
gene.name <- "Slco1a1"
gene.name <- "Cyp2j13"
gene.name <- "Abcc3"
gene.name <- "Hao2"
gene.name <- "Gsta2"

pdf(paste0("plots/Integrated_2Sex_ArchRPTpeaks_Coverage_", gene.name,"_PT-S3.pdf"), width = 6, height = 3)	
CoveragePlot(pt, region = gene.name, features = gene.name, assay = 'PTpeaks', expression.assay = 'RNA', peaks = TRUE, 
			extend.upstream = 0, extend.downstream = 5000, scale.factor = 1e7, idents=c("PT-S3-mWT","PT-S3-fWT","PT-S3-mKO"))
dev.off()
pdf(paste0("plots/Integrated_2Sex_ArchRPTpeaks_Coverage_", gene.name,"_PT-S2.pdf"), width = 6, height = 3)	
CoveragePlot(pt, region = gene.name, features = gene.name, assay = 'PTpeaks', expression.assay = 'RNA', peaks = TRUE, 
			extend.upstream = 0, extend.downstream = 5000, scale.factor = 1e7, idents=c("PT-S2-mWT","PT-S2-fWT","PT-S2-mKO"))
dev.off()
pdf(paste0("plots/Integrated_2Sex_ArchRPTpeaks_Coverage_", gene.name,"_PT-S1.pdf"), width = 6, height = 3)	
CoveragePlot(pt, region = gene.name, features = gene.name, assay = 'PTpeaks', expression.assay = 'RNA', peaks = TRUE, 
			extend.upstream = 0, extend.downstream = 5000, scale.factor = 1e7, idents=c("PT-S1-mWT","PT-S1-fWT","PT-S1-mKO"))
dev.off()

