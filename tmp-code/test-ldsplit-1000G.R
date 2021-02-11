library(dplyr)
library(ggplot2)
library(bigsnpr)

obj.bed <- bed(bedfile <- download_1000G(dir = "../bigsnpr/tmp-data/"))
fam2 <- bigreadr::fread2(sub_bed(bedfile, ".fam2"))
THR_R2 <- 0.05

CHR <- 22
pop <- "EAS"

ind.row <- `if`(pop == "all", rows_along(obj.bed), which(fam2$`Super Population` == pop))
ind.col <- which(obj.bed$map$chromosome == CHR)
obj.bigsnp <- snp_attach(snp_readBed2(bedfile, tempfile(), ind.row, ind.col, ncores = 4))
G <- obj.bigsnp$genotypes
ind.col2 <- which(snp_MAF(G, ncores = 4) > 0.05)
map <- obj.bigsnp$map[ind.col2, ]
POS2 <- snp_asGeneticPos(map$chromosome, map$physical.pos, dir = "../bigsnpr/tmp-data/")
corr <- runonce::save_run(
  snp_cor(G, ind.col = ind.col2, infos.pos = POS2, size = 3 / 1000,
          thr_r2 = THR_R2, ncores = 4),
  file = paste0("tmp-data/corr_", pop, "_chr", CHR, ".rds")
)
dim(corr)

# ldetect
url <- paste0(
  "https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7ff3e90be31f278a7b6a61daaba0dc/",
  pop, "/fourier_ls-chr", CHR, ".bed")
if (RCurl::url.exists(url)) {
  ldetect_chr <- bigreadr::fread2()
  n_block0 <- nrow(ldetect_chr)
  all_ind0 <- ldetect_chr$stop[-n_block0]
  block_num0 <- rowSums(outer(map$physical.pos, all_ind0, '>=')) + 1L
  stopifnot(length(block_num0) == ncol(corr))
  (cost0 <- bigsnpr:::compute_cost(block_num0, Matrix::tril(corr), THR_R2))
} else {
  n_block0 <- ncol(corr) / 600
  cost0 <- NA_real_
}



res <- runonce::save_run(
  bigsnpr::snp_ldsplit(corr, thr_r2 = THR_R2, min_size = 500, max_size = 10e3,
                       max_K = n_block0),
  file = paste0("tmp-data/split_", pop, "_chr", CHR, ".rds")
)
 # 11 sec for chr 22 (24,624 variants)
print(res, n = 20)

qplot(n_block, cost, data = res) +
  theme_bw(16) +
  scale_x_continuous(breaks = 0:20 * 10, minor_breaks = 0:100 * 2) +
  scale_y_log10() +
  geom_point(data = data.frame(cost = cost0, n_block = n_block0),
             color = "red", size = 3) +
  labs(y = "Sum of squared correlations outside blocks (log-scale)", x = "Number of blocks")
ggsave(paste0("figures/cost_split_", pop, "_chr", CHR, ".pdf"), width = 10, height = 7)

best_res <- subset(res, n_block == 20)
all_ind <- head(best_res$all_last[[1]], -1)
map$physical.pos[all_ind] / 1e6
block_num <- best_res$block_num[[1]]

ind_group <- split(seq_along(block_num), block_num)
ind_two <- lapply(seq(0, length(ind_group) - 2L), function(add) 1:2 + add)

STEP <- 200
all_cost <- list()

all_plots <- lapply(ind_two, function(ind2) {

  lims <- all_ind[ind2[1]] + 0.5

  all_cost[[ind2[1]]] <<- corr[ind_group[[ind2[1]]], ind_group[[ind2[2]]]] %>%
    { .@x^2 } %>%
    { sum(.[. >= THR_R2]) }

  ind <- unlist(ind_group[ind2], use.names = FALSE) %>%
    { .[between(., lims - STEP, lims + STEP)] }

  df2 <- corr[ind, ind] %>%
    as("dgTMatrix") %>%
    { data.frame(i = .@i + ind[1], j = .@j + ind[1], r2 = .@x^2) } %>%
    filter(r2 >= THR_R2)

  ggplot(df2, aes(x = i, y = j)) +
    geom_raster(aes(fill = r2)) +
    scale_fill_viridis_b(direction = -1, breaks = c(0.02, 0.1, 0.3, 0.8)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    geom_vline(xintercept = lims, linetype = 2) +
    geom_hline(yintercept = lims, linetype = 2) +
    scale_y_reverse(breaks = NULL) +
    coord_equal() +
    theme(axis.text.y = element_blank(), strip.text.x = element_blank())
})
(all_cost2 <- unlist(all_cost[seq_along(all_plots)]))
best_res$all_size[[1]]

cowplot::plot_grid(
  plotlist = c(lapply(all_plots, function(p) p + theme(legend.position = "none")),
               list(cowplot::get_legend(all_plots[[1]]))),
  nrow = 4
)

