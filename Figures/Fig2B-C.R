#### Load libraries ####
library(data.table)
library(dplyr)
library(edgeR)
library(DESeq2)
library(ggplot2)

### Read input files ###
df_ct <- fread("../Data/WashU_age_sex_raw_read_counts.txt")
df_mt <- fread("../Data/WashU_age_sex_meta.txt")

vec_age <- gsub("v","",df_mt$age)
vec_age <- gsub("b","",vec_age)
df_mt$Age <- vec_age
df_mt$Sex <- df_mt$gender

### Data processing using edgeR ###
mtx_data <- as.matrix(df_ct[,c(4:60)])
rownames(mtx_data) <- df_ct$symbol

d <- DGEList(counts=mtx_data,group=factor(df_mt$Sex))
keep <- rowSums(cpm(d)>5) >= 2 
length(which(keep)) 
d <- d[keep,]
dim(d)

df_ct_cpm <- filter(df_ct, symbol %in% rownames(d$counts))

### Run DESeq2 ###
df_ct <- df_ct_cpm
df_data <- as.matrix(df_ct[,c(4:60)])
rownames(df_data) <- df_ct$symbol
dds <- DESeqDataSetFromMatrix(countData = df_data, colData = df_mt, design = ~ Sex) 

### Experiment: code age as a continuous variable ###
vec_age <- df_mt$age
vec_age[which(vec_age=="0")] <- 1L
vec_age[which(vec_age=="2")] <- 2L
vec_age[which(vec_age=="4")] <- 3L
vec_age[which(vec_age=="8")] <- 4L
vec_age[which(vec_age %in% c("79b","79v"))] <- 5L
df_mt$time <- as.integer(vec_age)
dds <- DESeqDataSetFromMatrix(countData = df_data, colData = df_mt, design = ~ time)
#################################################

### PCA plot ###
### rlog transformation ###
dds <- estimateSizeFactors(dds) 
rld <- rlog(dds, blind=FALSE)

# Option 1: use built-in plot function in DESeq2
plotPCA(rld, intgroup = c("Age", "Sex")) 

# Option 2: plot with ggplot2 #
pcadata <- plotPCA(rld, intgroup = c("Age", "Sex"), ntop = 500, returnData=TRUE)
pcadata$Age <- factor(pcadata$Age, levels = c("0","2","4","8","79"))
percentVar <- round(100 * attr(pcadata, "percentVar"))

ggplot(pcadata, aes(PC1, PC2, color=Age, shape=Sex)) + 
  geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#################################################

### Differential expression analysis using DESeq2 ###
dds <- DESeq(dds)
res <- results(dds) 
summary(res) ### check direction of comparison
df_res <- as.data.frame(res[order(res$padj),])

## Extract DEGs ###
df_res_sig <- filter(df_res, abs(log2FoldChange)>=0.5 & padj < 0.05) ### adjust thresholds here
df_res_sig$symbol <- rownames(df_res_sig) ### DEG list for use

### Annotate gene functions (optional) ###
mart <- readRDS(file = "../Data/ensembl_mouse_genes.rds")
gene_set <- df_res_sig$symbol
attr = c("external_gene_name","chromosome_name", "description")
gene_anno <- getBM(attributes = attr, filters = "external_gene_name", values = gene_set, mart = mart, useCache = FALSE)        
colnames(gene_anno)[c(1,2)] <- c("symbol","Chr")
df_res_out <- merge(gene_anno,df_res_sig,by="symbol",all.y=T)

### Annotate DE directions ###
df_res_out$Direction <- unlist(sapply(df_res_out$log2FoldChange, function(x){ifelse(x>0,"M","F")}))
table(df_res_out$Direction)

### Write out DEG list ###
write.table(df_res_out,"{OUTPUT_NAME}.txt",row.names = F,quote = F, sep = "\t") 

### Bar plot: sex-specific genes over time ###
data <- data.frame(Sex=rep(c("F", "M"), each=5),
                   Group=rep(c("0", "2", "4", "8", "79"),2),
                   No_Gene=c(2,2,160,869,862,0,4,307,864,693))
data$Group <- factor(data$Group, levels = c("0", "2", "4", "8", "79"))

ggplot(data=data, aes(x=Group, y=No_Gene, fill=Sex)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("No. of DEGs") + xlab("Age (week)") +
  scale_fill_manual(values=c("#FF8596","#7CADED")) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

