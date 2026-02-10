### Load libraries ###
library(data.table)
library(dplyr)
library(Signac)
library(Seurat)
library(ggplot2)
set.seed(1234)

#### Download processed single-nuclear multiomic dataset ####
# Link: https://ucla.box.com/s/6f43jpkv5sco1bmcv8ttmp19edz9pnaw

### Functions ###
compute_sseq_params <- function (x, zeta_quantile = 0.995)
{
    params <- list()
    N <- nrow(x)
    G <- ncol(x)
    grand_median <- median(rowSums(x))
    s_ij <- rowSums(x)/grand_median
    params$s_ij <- s_ij
    x_size_norm <- x/s_ij
    mu_g <- colMeans(x_size_norm)
    v_g <- apply(x_size_norm, 2, var)
    use_g <- v_g > 0
    params$mu_g <- mu_g
    params$v_g <- v_g
    params$use_g <- use_g
    phi_mm_g <- pmax(0, (N * v_g - mu_g * sum(1/s_ij))/(mu_g^2 * sum(1/s_ij)))
    params$phi_mm_g <- phi_mm_g
    zeta_hat <- quantile(params$phi_mm_g, zeta_quantile, na.rm = T)
    params$zeta_hat <- zeta_hat
    mean_phi_mm_g <- mean(phi_mm_g[use_g])
    delta <- (sum((phi_mm_g[use_g] - mean_phi_mm_g)^2)/(G - 1))/(sum((phi_mm_g[use_g] - zeta_hat)^2)/(G - 2))
    params$delta <- delta
    phi_g <- rep(NA, G)
    phi_g[use_g] <- (1 - delta) * phi_mm_g[use_g] + delta * zeta_hat
    params$phi_g <- phi_g
    params
}

sseq_differential_expression <- function (x, cond0, cond1, sseq_params, gene_symbols){
    x_a <- x[cond0, , drop = F]
    x_b <- x[cond1, , drop = F]
    G <- ncol(x)
    s_a <- sum(sseq_params$s_ij[cond0])
    s_b <- sum(sseq_params$s_ij[cond1])
    x_ga <- colSums(x_a)
    x_gb <- colSums(x_b)
    p_res <- rep(NA, G)
    p_res[sseq_params$use_g] <- sapply(which(sseq_params$use_g),
                                       function(g) nb_exact_test(x_ga[g], x_gb[g], s_a, s_b,
                                                                 sseq_params$mu_g[g], sseq_params$phi_g[g]))
    p_adj <- rep(1, G)
    p_adj[sseq_params$use_g] <- p.adjust(p_res[sseq_params$use_g], method = "BH")
    data.frame(symbol = gene_symbols, tested = sseq_params$use_g, 
		sum_a = x_ga, sum_b = x_gb,
		common_mean = sseq_params$mu_g, dispersion = sseq_params$phi_g,
		mean_a_sizenorm = x_ga/s_a, mean_b_sizenorm = x_gb/s_b,
		log2FC = log2((1 + x_ga)/(1 + s_a)) - log2((1 + x_gb)/(1 + s_b)), 
		p = p_res, p_adj = p_adj)
}

nb_exact_test <- function (x_a, x_b, s_a, s_b, mu, phi)
{
    all_x_a <- seq(0, x_a + x_b, 1)
    all_x_b <- seq(x_a + x_b, 0, -1)
    .prob <- function(x, s) dnbinom(x, mu = s * mu, size = 1/(phi/s))
    p_obs <- .prob(x_a, s_a) * .prob(x_b, s_b)
    p_all <- .prob(all_x_a, s_a) * .prob(all_x_b, s_b)
    sum(p_all[p_all <= p_obs])/sum(p_all)
}

### Run SSeq ###
pt <- readRDS("Integrated_ArchRpeaks_PT.rds")
DefaultAssay(pt) <- "RNA"
t_mat <- t(as.matrix(pt@assays$RNA@counts))
sseq_params <- compute_sseq_params(t_mat)
fdata <- rownames(pt[["RNA"]])
cluster_labels <- pt$branch 
group0 <- cluster_labels == "KO"
group1 <- cluster_labels == "WT"
de_result <- sseq_differential_expression(t_mat, group0, group1, sseq_params, gene_symbols = fdata)

min_mean <- 0.01; p_cutoff <- 0.05; min_log2fc <- 0.5
de_result$significant <- with(de_result, common_mean >= min_mean & p_adj < p_cutoff & abs(log2FC) >= min_log2fc)
de_result$Direction <- unlist(sapply(de_result$log2FC, function(x){ifelse(x>0,"UP","DOWN")}))
de_result <- de_result[with(de_result, order(Direction, p_adj)),]
write.csv(de_result,"{OUTPUT_NAME}_DEGs_KOvsWT.csv",row.names=F,quote=F) 
