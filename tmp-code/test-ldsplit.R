library(dplyr)
library(ggplot2)

# chr 22
corr <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/24928052", dir = "../bigsnpr/tmp-data"))
CHR <- 22

# chr 12
corr <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/24927977", dir = "../bigsnpr/tmp-data"))
CHR <- 12

# chr 6
corr <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/24927920", dir = "../bigsnpr/tmp-data"))
CHR <- 6

# chr 2
corr <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/24927572", dir = "../bigsnpr/tmp-data"))
CHR <- 2


map <- filter(readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/25503788", dir = "../bigsnpr/tmp-data")),
  chr == CHR)
# POS2 <- with(map, bigsnpr::snp_asGeneticPos(chr, pos, dir = "../bigsnpr/tmp-data"))
# stopifnot(length(POS2) == ncol(corr))

THR_R2 <- 0.01

# ldetect
ldetect_chr <- bigreadr::fread2(paste0(
  "https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7ff3e90be31f278a7b6a61daaba0dc/EUR/fourier_ls-chr", CHR, ".bed"))
all_ind0 <- head(ldetect_chr$stop, -1)
n_block0 <- length(all_ind0) + 1L
block_num0 <- rowSums(outer(map$pos, all_ind0, '>=')) + 1L
(cost0 <- bigsnpr:::compute_cost(block_num0, Matrix::tril(corr), THR_R2))


system.time(
  res <- bigsnpr::snp_ldsplit(
    corr,
    thr_r2 = THR_R2,
    min_size = 500,
    max_size = 10e3,
    max_K = n_block0
  )
) # 121 sec for chr2 (87,521 variants)
print(res[1:5], n = 20)

qplot(n_block, cost, data = res) +
  theme_bw(16) +
  scale_x_continuous(breaks = 0:20 * 10, minor_breaks = 0:100 * 2) +
  scale_y_log10() +
  geom_point(data = data.frame(cost = cost0, n_block = n_block0),
             color = "red", size = 3) +
  labs(y = "Sum of squared correlations outside blocks (log-scale)", x = "Number of blocks")
# ggsave("cost_chr2.pdf", width = 10, height = 7)

best_res <- subset(res, n_block == 50)
all_ind <- head(best_res$all_last[[1]], -1)
map$pos[all_ind] / 1e6
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
  nrow = 7
)
