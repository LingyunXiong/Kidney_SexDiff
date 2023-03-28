### Load libraries ###
library(data.table)
library(dplyr)
library(Signac)
library(Seurat)
library(ggplot2)
set.seed(1234)

library(RColorBrewer)
my_cols = brewer.pal(4,"Spectral")

#### Download processed single-nuclear multiomic dataset ####
# Link:

### Load merged object ###
Merged_f2 <- readRDS("Integrated_2Sex_ArchRpeaks_Filtered.rds")
DefaultAssay(Merged_f2) <- "RNA"

Merged_f2$dataset <- Merged_f2$orig.ident
Merged_f2$type <- unlist(sapply(Merged_f2$dataset, function(x){ifelse(x %in% c("A1","A2","B1","B3","B5"),"WT","KO")}))
Merged_f2$sex <- unlist(sapply(Merged_f2$dataset, function(x){ifelse(x %in% c("A1","A2","A3","A4"),"F","M")}))

### Split the dataset into a list of two seurat objects (M and F) ###
Merged_list <- SplitObject(Merged_f2, split.by = "sex")
rm(Merged_f2); gc()

### Normalize and identify variable features for individual datasets ###
Merged_list <- lapply(X = Merged_list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

### Select variable features that are shared across dataset ###
features <- SelectIntegrationFeatures(object.list = Merged_list, nfeatures = 2000)
Merged_list <- lapply(X = Merged_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

### Perform integration ###
integration.anchors <- FindIntegrationAnchors(
        object.list = Merged_list,
        normalization.method = "LogNormalize",
        anchor.features = features,
        dims = 1:30,
        reduction = "rpca",
        k.anchor = 20)

rm(Merged_list); gc()
Merged_combined <- IntegrateData(
        anchorset = integration.anchors,
        normalization.method = "LogNormalize",
        new.assay.name = "integratedRNA",
        dims = 1:30)

### Clustering on integrated dataset ###
Merged_combined$group <- paste(Merged_combined$sex,Merged_combined$type,sep="-")
DefaultAssay(Merged_combined) <- "integratedRNA"

# Run the standard workflow for visualization and clustering
Merged_combined <- ScaleData(Merged_combined, verbose = FALSE)
Merged_combined <- RunPCA(Merged_combined, npcs = 30, verbose = FALSE)
Merged_combined <- FindNeighbors(Merged_combined, reduction = "pca", dims = 1:30)
Merged_combined <- FindClusters(Merged_combined, resolution = 0.5)
Merged_combined <- RunUMAP(Merged_combined, reduction = "pca", dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# visualization
p1 <- DimPlot(Merged_combined, reduction = "umap.rna", label = TRUE) + ggtitle('RNA Clustering') & NoLegend()
new.cell.order <- sample(rownames(Merged_combined[[]]),dim(Merged_combined[[]])[1])
p2 <- DimPlot(Merged_combined, reduction = "umap.rna", group.by = 'group', cols=alpha(my_cols, 0.5), cells = new.cell.order) + ggtitle('Color by sex-genotype')

pdf("Integrated_ArchRpeaks_UMAP_RNA_LogNrpca_K20_4Group_r5.pdf", width = 11, height = 5)
p1 + p2
dev.off()

### Integrate ArchR peak data ###
DefaultAssay(Merged_combined) <- "peaks"
Merged_combined <- RunTFIDF(Merged_combined)
Merged_combined <- FindTopFeatures(Merged_combined, min.cutoff = "q25")
Merged_combined <- RunSVD(Merged_combined)
Merged_combined <- RunUMAP(Merged_combined, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

p3 <- DimPlot(Merged_combined, reduction = "umap.atac", group.by = 'group', cols=alpha(my_cols, 0.5), cells = new.cell.order) + ggtitle('Peaks Clustering (Merged)')

#integrate embeddings and output new object to prevent overwriting integrated RNA
Merged_combined_atac <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  new.reduction.name = "integratedLSI",
  reductions = Merged_combined@reductions$lsi
)

#copy integrated LSI from duplicate seurat object to original object
Merged_combined@reductions$integratedLSI <- Merged_combined_atac@reductions$integratedLSI
Merged_combined <- RunUMAP(Merged_combined, reduction = 'integratedLSI', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
rm(Merged_combined_atac); gc()

p4 <- DimPlot(Merged_combined, reduction = "umap.atac", group.by = 'group', cols=alpha(my_cols, 0.5), cells = new.cell.order) + ggtitle('Peaks Clustering (Integrated)')

pdf("Integrated_ArchRpeaks_UMAP_peaks_LogNrpca_K20_4Group_r5.pdf", width = 10, height = 5)
p3 + p4
dev.off()

### WNN ###
Merged_combined <- FindMultiModalNeighbors(Merged_combined, reduction.list = list("pca", "integratedLSI"), dims.list = list(1:30, 2:30))
Merged_combined <- FindClusters(Merged_combined, resolution = 0.5, graph.name = "wsnn", algorithm = 3) ### original resolution: 0.5
Merged_combined <- RunUMAP(Merged_combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

p5 <- DimPlot(Merged_combined, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p6 <- DimPlot(Merged_combined, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p7 <- DimPlot(Merged_combined, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

pdf("Integrated_ArchRpeaks_UMAP_WNN_LogNrpca_K20_4Group_r5.pdf", width = 14, height = 5)
p5 + p6 + p7 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

### Save object ###
saveRDS(Merged_combined, file = "Integrated_ArchRpeaks_LogNorm_rpca_K20.rds")
write.csv(Merged_combined[[]],"Integrated_ArchRpeaks_LogNorm_meta_WNN.csv",quote=F)
