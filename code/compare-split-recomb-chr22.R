library(bigsnpr)
library(dplyr)
library(ggplot2)

obj.bed <- bed(bedfile <- download_1000G(dir = "../bigsnpr/tmp-data/"))
fam2 <- bigreadr::fread2(sub_bed(bedfile, ".fam2"))
THR_R2 <- 0.05

CHR <- 22


# ldetect
url <- paste0(
  "https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7ff3e90be31f278a7b6a61daaba0dc/",
  "EUR/fourier_ls-chr", CHR, ".bed")
ldetect_chr <- bigreadr::fread2(url)
n_block0 <- nrow(ldetect_chr)
pos_bounds0 <- ldetect_chr$stop[-n_block0]


# prepare for snp_splitld
ind.row <-  which(fam2$`Super Population` == "EUR")
ind.col <- which(obj.bed$map$chromosome == CHR)
obj.bigsnp <- snp_attach(snp_readBed2(bedfile, tempfile(), ind.row, ind.col, ncores = 4))
G <- obj.bigsnp$genotypes
ind.col2 <- which(snp_MAF(G, ncores = 4) > 0.05)
map <- obj.bigsnp$map[ind.col2, ]
POS <- map$physical.pos
POS2 <- snp_asGeneticPos(map$chromosome, POS, dir = "../bigsnpr/tmp-data/")
corr <- snp_cor(G, ind.col = ind.col2, infos.pos = POS2, size = 3 / 1000,
                thr_r2 = THR_R2, ncores = 4)
dim(corr)

# snp_splitld
res <- bigsnpr::snp_ldsplit(
  corr, thr_r2 = THR_R2, min_size = 200, max_size = 3000, max_K = n_block0)

qplot(data = res, n_block, cost) +
  theme_bw(12) +
  scale_x_continuous(breaks = 0:100 * 10, minor_breaks = 0:500 * 2) +
  scale_y_log10() +
  theme(legend.position = "top") +
  labs(y = "Sum of squared correlations outside blocks (log-scale)", x = "Number of blocks")

id <- head(filter(res, n_block == n_block0)$all_last[[1]], -1)
pos_bounds <- (POS[id] + POS[id + 1]) / 2


# https://drive.google.com/drive/folders/1Tgt_7GsDO0-o02vcYSfwqHFd3JNF6R06
# untar("tmp-data/hg19_maps.tar.gz", exdir = "tmp-data")
all_maps <- paste0(
  "tmp-data/hg19/CEU/CEU_recombination_map_hapmap_format_hg19_chr_", CHR, ".txt")
readLines(all_maps, n = 5)
recomb_map <- bigreadr::fread2(all_maps, col.names = c("chr", "pos", "rate", "gen_pos"))

levels(cut(recomb_map$pos, 6))
breaks <- c(0, 2.2, 2.8, 3.4, 4, 4.5, Inf) * 1e7

qplot(pos, bigutilsr::rollmean(rate, 20), geom = "line", data = recomb_map) +
  theme_bw(15) +
  geom_vline(aes(xintercept = pos), data = data.frame(pos = pos_bounds),
                 linetype = 2, color = "#00BFC4") +
  geom_vline(aes(xintercept = pos), data = data.frame(pos = pos_bounds0),
             linetype = 2, color = "#F8766D") +
  labs(x = "Physical position (bp)", y = "Recombination rate (smoothed)") +
  facet_wrap(~ cut(pos, breaks), scales = "free_x", ncol = 1) +
  theme(strip.background = element_blank(), strip.text.x = element_blank())
# ggsave("figures/compare_recomb_chr22.pdf", width = 10, height = 10)
