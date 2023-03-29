### Load libraries ###
library(data.table)
library(dplyr)
library(ArchR)
library(dineq)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
library(ggplot2)

#### GAS calculation functions ####
Extend <- function(x, upstream = 0, downstream = 0, from.midpoint = FALSE){
   if (any(strand(x = x) == "*")) {
      warning("'*' ranges were treated as '+'")
   }
   on_plus <- strand(x = x) == "+" | strand(x = x) == "*"
   if (from.midpoint) {
      midpoints <- start(x = x) + (width(x = x) / 2)
      new_start <- midpoints - ifelse(
         test = on_plus, yes = upstream, no = downstream
      )
      new_end <- midpoints + ifelse(
         test = on_plus, yes = downstream, no = upstream
      )
   } else {
      new_start <- start(x = x) - ifelse(
         test = on_plus, yes = upstream, no = downstream
      )
      new_end <- end(x = x) + ifelse(
         test = on_plus, yes = downstream, no = upstream
      )
   }
   ranges(x = x) <- IRanges(start = new_start, end = new_end)
   x <- trim(x = x)
   return(x)
}

find.peak.regions <- function(archrObject, GENE.info, Upstream.Extend) {
   check.dist <- Upstream.Extend
   while (check.dist > 0) {
      TSS_checkregion <- Extend(resize(GENE.info, 1), upstream = check.dist, downstream = -1)
      
      
      overlapping_genes_idx <- findOverlaps(TSS_checkregion, archrObject@geneAnnotation$genes[-gene.idx,])
      if (length(overlapping_genes_idx) != 0) {
         overlapping_gene_region <- archrObject@geneAnnotation$genes[-gene.idx,][subjectHits(overlapping_genes_idx),]
         check.dist <- GenomicRanges::distance(GENE.info, overlapping_gene_region) - 1
      } else {
         check.dist <- 0
      }
   }
   gene.info_extended <- union(GENE.info, TSS_checkregion)
   idx_peaks_in_gene.info_extended <- queryHits(findOverlaps(proj@peakSet, gene.info_extended))
   return(idx_peaks_in_gene.info_extended)
}

accessibility.of.cluster_archr <- function(archrObject, celltype, peak_index){
   # rows are genomic regions, columns are clusters
   sum.access.data <- matrix(data = 0, nrow = length(peak_index), ncol = length(celltype))
   # for each accessible region we are looking at...
   # accessibility of each region within each cell
   access.counts <- assays(peakMatrix)$PeakMatrix[peak_index,, drop =F]
   # access.counts <- access.counts/total_counts_perCell
   for (i in 1:length(celltype)){
      access.counts.s <- access.counts[,colnames(access.counts) %in% archrObject$cellNames[which(archrObject$Clusters_ivy_v4_refined==celltype[i])],drop = F]
      # average over the top numCellsAvg cells
      sum.access.data[,i] <- rowSums(access.counts.s)/totalCounts_perCluster[i]
   }
   return(sum.access.data)
}

z.Gini <- function(region.access){ # Calculate Gini score and z normalize
   gini.scores <- rep(0, nrow(region.access))
   for (i in 1:nrow(region.access)){
      if (!is.na(gini.wtd(region.access[i,]))){gini.scores[i] <- gini.wtd(region.access[i,])}
   }
   # z-standardize
   if (nrow(region.access) > 1){
      gini.mean <- mean(gini.scores)
      gini.sdev <- sd(gini.scores)
      gini.scores <- (gini.scores - gini.mean)/gini.sdev
   } else {
      # if there is only one peak, just set gini score to 0
      gini.scores[1] <- 0
   }
   return(gini.scores)
}

gene.access.score <- function(region.access, final.weights){ # Compute gene accessibility score
   gene.scores <- rep(0, ncol(region.access))
   for (z in 1:ncol(region.access)){
      # using region.access which average of accessibility counts
      # normalizing by number of peaks == length(final.weights)
      gene.scores[z] <- sum(region.access[,z] * final.weights)/length(final.weights)
   }
   return(gene.scores)
}

#### GAS setting up objects ####
peakMatrix <- getMatrixFromProject(
   ArchRProj = proj_ls[[1]],
   useMatrix = "PeakMatrix",
   useSeqnames = NULL,
   verbose = F,
   binarize = FALSE,
   threads = getArchRThreads(),
   logFile = NULL
) 
proj <- proj_ls[[1]]

vec_cluster <- names(table(proj$Clusters_ivy_v4_refined))
sample_list <- c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "B5", "B6")
total_counts_perCell <- colSums(assays(peakMatrix)$PeakMatrix)
totalCounts_perCluster <- c()
for (i in 1:length(vec_cluster)){
   tmp <- total_counts_perCell[names(total_counts_perCell) %in% proj$cellNames[which(proj$Clusters_ivy_v4_refined==vec_cluster[i])]]
   totalCounts_perCluster[i]<- sum(tmp)
}

#### Calculate Jenssens Gene Accessibility Score ####
vec_DEG <- read.csv(file = "Multiome_PT_fWTvsmWT_DEGs_ARresponse.csv", header = T, row.names = 1, stringsAsFactors = F)
vec_genes <- rownames(vec_DEG)
upstream.extend <- 5000
mtx_gene_GAS <- matrix(data = NA, nrow = length(vec_genes), ncol = length(vec_cluster)+1)

# Get gene range and extend it 5k upstream
for (i in 1:length(vec_genes)) {
   gene.name <- vec_genes[i]
   gene.idx <- which(proj@geneAnnotation$genes$symbol == gene.name)
   gene.info <- proj@geneAnnotation$genes[gene.idx,]
   if (length(gene.info) == 0) {print("gene not found in annotation");next}
   gene.info_extended <- Extend(gene.info, upstream = upstream.extend, downstream = -2)
   gene.width <- gene.info@ranges@width
   
   # Find if the extended gene range is overlapping with other genes
   idx_peaks_in_gene.info_extended <- find.peak.regions(proj, gene.info, upstream.extend)
   overlapping_gr <- proj@geneAnnotation$genes[subjectHits(findOverlaps(gene.info_extended, proj@geneAnnotation$genes)),]
   if (length(idx_peaks_in_gene.info_extended) == 0) {print("no peaks in gene region");next}
	   
   # Calculate average counts in each peak # 
   region.access <- accessibility.of.cluster_archr(proj, vec_cluster, idx_peaks_in_gene.info_extended)
   
   # using peaks that overlap with extended gene body region,
   # calculate distance from each peak to the gene TSS 
   overlapping.granges <- proj@peakSet[idx_peaks_in_gene.info_extended,]
   dist.to.TSS <- rep(0, length(idx_peaks_in_gene.info_extended))
   
   for (j in 1:length(overlapping.granges@seqnames)){
      if(as.character(gene.info_extended@strand) == "+") {
         # if start of genomic region upstream of gene, dist = end - tss
         if (overlapping.granges@ranges@start[j] < gene.info@ranges@start){
            # and end of genomic region is still upstream of gene,
            if (overlapping.granges@ranges@start[j]+overlapping.granges@ranges@width[j]-1 < gene.info@ranges@start) {
               dist.to.TSS[j] <- abs(overlapping.granges@ranges@start[j]+overlapping.granges@ranges@width[j]-1 - gene.info@ranges@start)
            } else { # if TSS is within genomic region, find shortest distance
               dist.to.TSS[j] <- 0
            }
         } else { # else genomic region starts downstream of TSS, dist = start - tss
            dist.to.TSS[j] <- abs(gene.info@ranges@start - overlapping.granges@ranges@start[j])
         }
      } else { # for negative strand gene
         # genomic region starts downstream of tss,
         if (overlapping.granges@ranges@start[j] < gene.info@ranges@start+gene.info@ranges@width-1){
            # and ends downstream of tss
            if(overlapping.granges@ranges@start[j]+overlapping.granges@ranges@width[j]-1 < gene.info@ranges@start+gene.info@ranges@width-1) {
               dist.to.TSS[j] <- gene.info@ranges@start+gene.info@ranges@width-1 - overlapping.granges@ranges@start[j]+overlapping.granges@ranges@width[j]-1
            } else { # if TSS is within genomic region, find shortest distance
               dist.to.TSS[j] <- 0
            }
         } else { # else genomic region starts downstream of TSS, dist = start - tss
            dist.to.TSS[j] <- abs(gene.info@ranges@start+gene.info@ranges@width-1 - overlapping.granges@ranges@start[j])
         }
      }
   }
   
   w.d <- exp(-dist.to.TSS/5000) + exp(-1) 
   gini.scores <- z.Gini(region.access)
   w.b <- exp(gini.scores)
   final.weights <- w.d*w.b
   final.gene.scores <- gene.access.score(region.access, final.weights) *10^5
   mtx_gene_GAS[i,] <- c(final.gene.scores, gene.width)
   print(paste0(gene.name, " is done"))
}

rownames(mtx_gene_GAS) <- vec_genes
colnames(mtx_gene_GAS) <- c(vec_cluster, "gene_length")
rm(i, gene.name, gene.idx, gene.info, gene.info_extended, gene.width, idx_peaks_in_gene.info_extended, overlapping_gr, region.access)
rm(overlapping.granges, dist.to.TSS, w.d, gini.scores, w.b, final.weights, final.gene.scores)

#### Write out GAS results ####
write.csv(x = mtx_gene_GAS, file = "Kidney_snMultiomic_GAS.csv")
rm(df, p)


