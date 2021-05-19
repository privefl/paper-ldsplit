sumstats <- bigreadr::fread2(runonce::download_file(
  "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/all_sumstats/PASS_ADHD_Demontis2018.sumstats",
  dir = "tmp-data/sumstats"))

map_ldref <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/25503788",
  dir = "tmp-data", fname = "map_ldref.rds"))

res <- tibble::tibble(chr = 1:22)
res$split <- lapply(res$chr, function(chr) {
  print(chr)
  corr <- readRDS(paste0("../../Downloads/LD_REF/ld-ref/LD_chr", chr, ".rds"))
  res <- bigsnpr::snp_ldsplit(corr, thr_r2 = 0.05, min_size = 3000, max_size = 10e3, max_K = 50)
})

library(ggplot2)
qplot(data = tidyr::unnest(res, "split"), n_block, cost) +
  theme_bw(12) +
  scale_x_continuous(breaks = 0:20 * 5, minor_breaks = 0:100) +
  scale_y_log10() +
  facet_wrap(~ chr, nrow = 7, scales = "free_x") +
  theme(legend.position = "top") +
  labs(y = "Sum of squared correlations outside blocks (log-scale)",
       x = "Number of blocks")

blocks_chosen <- data.frame(
  chr = 1:22,
  n_block = c(20, 20, 16, 15, 16, 14, 15, 12, 10, 11,
              12, 11, 9, 8, 8, 8, 7, 9, 6, 7, 4, 4)
)
sum(blocks_chosen$n_block) # 242


res_all <- tidyr::unnest(res, "split")
qplot(data = res_all, n_block, cost) +
  theme_bw(14) +
  scale_x_continuous(breaks = 0:20 * 5, minor_breaks = 0:100) +
  scale_y_log10() +
  facet_wrap(~ chr, nrow = 7, scales = "free_x") +
  theme(legend.position = "top") +
  labs(y = "Sum of squared correlations outside blocks (log-scale)",
       x = "Number of blocks") +
  geom_vline(aes(xintercept = n_block), data = blocks_chosen,
             color = "red", linetype = 3)
# ggsave("figures/cost_split_ldref.pdf", width = 13, height = 10)

blocks <- dplyr::semi_join(res_all, blocks_chosen)
blocks$cost

offset <- 0
block_num <- unlist(lapply(blocks$block_num, function(num) {
  res <- num + offset
  offset <<- offset + tail(num, 1)
  res
}))
tail(block_num)  # 242
rle(block_num)

map_ldref$block <- block_num
# saveRDS(map_ldref, "map_blocks.rds")
