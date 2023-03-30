#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(tidyr)

mean_score <- function(data){
  mtx <- as.matrix(data)
  vec_mean <- rowMeans(mtx)
  return(vec_mean)
}

### read in DoRothEA output: scaled activity score ###
df_8wk_core <- fread("../Data/SexGenes_8wk_Core_DoRothEA_TF_activity_score.csv")
df_8wk_core$MeanAS_F <- mean_score(df_8wk_core[,c(2:6)])
df_8wk_core$MeanAS_M <- mean_score(df_8wk_core[,c(7:11)])
df_8wk_core_s <- select(df_8wk_core,symbol,MeanAS_F,MeanAS_M)

df_PT <- fread("../Data/PTgenes_8wk_GE_DoRothEA_TF_activity_score.csv")
df_PT$MeanAS_F <- mean_score(df_PT[,c(2:6)])
df_PT$MeanAS_M <- mean_score(df_PT[,c(7:11)])
df_PT_s <- select(df_PT,symbol,MeanAS_F,MeanAS_M)

df_4wk <- fread("../Data/SexGenes_4wk_DoRothEA_TF_activity_score.csv")
df_4wk$MeanAS_F <- mean_score(df_4wk[,c(2:6)])
df_4wk$MeanAS_M <- mean_score(df_4wk[,c(7:11)])
df_4wk_s <- select(df_4wk,symbol,MeanAS_F,MeanAS_M)

### select TFs ###
vec2 <- sort(intersect(df_8wk_core$symbol,df_PT$symbol))
vec_tf <- c(vec2,"Esr1")

### prepare plot data frame ###
df_8wk_core_t <- filter(df_8wk_core_s, symbol %in% vec_tf)
colnames(df_8wk_core_t)[2:3] <- c("8wk-Core-F","8wk-Core-M")
df_8wk_core_p <- gather(df_8wk_core_t, group, NES, `8wk-Core-F`:`8wk-Core-M`, factor_key=TRUE)

df_PT_t <- filter(df_PT_s, symbol %in% vec_tf)
colnames(df_PT_t)[2:3] <- c("PT-F","PT-M")
df_PT_p <- gather(df_PT_t, group, NES, `PT-F`:`PT-M`, factor_key=TRUE)

df_4wk_t <- filter(df_4wk_s, symbol %in% vec_tf)
colnames(df_4wk_t)[2:3] <- c("4wk-F","4wk-M")
df_4wk_p <- gather(df_4wk_t, group, NES, `4wk-F`:`4wk-M`, factor_key=TRUE)

df_plot <- rbind(df_8wk_core_p,df_PT_p,df_4wk_p)
vec_tf_order <- arrange(df_8wk_core_t,`8wk-Core-F`)$symbol
vec_tf_order <- c(vec_tf_order[1:9],"Esr1",vec_tf_order[10:19])
df_plot$symbol <- factor(df_plot$symbol, levels = vec_tf_order)
df_plot$group <- factor(df_plot$group, levels = c("4wk-F","8wk-Core-F","PT-F","4wk-M","8wk-Core-M","PT-M"))

### tile heatmap ###
ggplot(df_plot, aes(x = symbol, y = group, fill = NES)) +
  geom_tile(color = "white", lwd = 1, linetype = 1) + coord_fixed() +
  scale_fill_gradient2(low = "darkblue", high = "indianred", mid = "white", midpoint = 0) + 
  theme_minimal() + xlab("") + ylab("") + coord_flip() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
