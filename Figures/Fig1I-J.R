library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

### ToppCluster ###
df_pathway <- fread("ToppFun/ToppCluster_output_Pathway_FDR_MinGene5.tsv")
df_pathway_s <- filter(df_pathway,`Core-F_logP`>6 | `Core-M_logP`>6 | `All-M_logP`>6 | `All-F_logP`>6)
df_pathway_s <- filter(df_pathway_s, !ID %in% c("132956","131226","PW:0000065","PW:0000485"))

table_process <- function(data,col_p,col_gene){
  data_s <- select(data,`Title (or Source)`,`col_p`,`col_gene`)
  colnames(data_s) <-c("Pathway","-log10(FDR)","GeneSet")
  vec_mol <- as.character(data_s$GeneSet)
  data_s$`Gene Count` <- unlist(sapply(vec_mol, function(x){length(unlist(strsplit(x, ",")))}))
  return(data_s)
}

df_F_full <- table_process(df_pathway_s,"All-F_logP","All-F_GeneSet")
df_F_full$group <- rep("All-F",nrow(df_F_full))
df_F_core <- table_process(df_pathway_s,"Core-F_logP","Core-F_GeneSet")
df_F_core$group <- rep("Core-F",nrow(df_F_core))
df_M_full <- table_process(df_pathway_s,"All-M_logP","All-M_GeneSet")
df_M_full$group <- rep("All-M",nrow(df_M_full))
df_M_core <- table_process(df_pathway_s,"Core-M_logP","Core-M_GeneSet")
df_M_core$group <- rep("Core-M",nrow(df_F_core))

df_plot <- rbind(df_F_core, df_F_full, df_M_core, df_M_full)
vec_go_order <- arrange(df_plot,-`-log10(FDR)`)$Pathway
vec_go_order <- vec_go_order[!duplicated(vec_go_order)]
df_plot$term <- factor(df_plot$Pathway, levels = rev(vec_go_order))
df_plot$group <- factor(df_plot$group, levels = c("All-F","Core-F","All-M","Core-M"))

ggplot(data = df_plot, aes(x = group, y = term, color = `-log10(FDR)`, size = `Gene Count`)) + 
  geom_point() + theme_bw() + scale_size(range = c(0, 8)) +
  scale_color_gradient(low = "grey", high = "red") +
  ylab("") + xlab("") + ggtitle("ToppFun pathway") + 
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face= "bold"))
