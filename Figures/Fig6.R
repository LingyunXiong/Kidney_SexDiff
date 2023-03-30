#### Load libraries ####
library(data.table)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(tidyr)
library(RColorBrewer)
library(gplots)

#### Read in and scale input data ####
# Fig. 6D #
df_z <- fread("../Data/CandidateGenes_KCE_Expression.csv") 
mtx <- as.matrix(df_z[, c(2:7)])
rownames(mtx) <- df_z$symbol
mtx_scale <- t(scale(t(mtx)))
df_p <- gather(as.data.frame(mtx_scale), group, Level, mS1:fS3, factor_key=TRUE)

# Fig. 6E #
df_z <- fread("../Data/CandidateGenes_RNAscope_Levels.csv")
mtx <- as.matrix(df_z[, c(2:10)])
rownames(mtx) <- df_z$symbol
mtx_scale <- t(scale(t(mtx)))
df_p <- gather(as.data.frame(mtx_scale), group, Level, mS1:mSix2KO_S3, factor_key=TRUE)

# Fig. 6F #
df_z <- fread("../Data/CandidateGenes_Multiome_GAS.csv")
mtx <- as.matrix(df_z[, c(2:10)])
rownames(mtx) <- df_z$symbol
mtx_scale <- t(scale(t(mtx)))
df_p <- gather(as.data.frame(mtx_scale), group, Level, PT.S1.mWT:PT.S3.mKO, factor_key=TRUE)

### tile heatmap ###
df_p$symbol <- factor(df_z$symbol,levels = rev(c("Atp11a","Slco1a1","Cyp2j13","Abcc3","Gsta4","Hao2")))

ggplot(df_p, aes(x = symbol, y = group, fill = Level)) +
  geom_tile(color = "white", lwd = 1, linetype = 1) + coord_fixed() +
  scale_fill_gradient2(low = "grey", high = "red") + 
  theme_minimal() + xlab("") + ylab("") + coord_flip() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

