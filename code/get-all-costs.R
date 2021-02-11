library(dplyr)
library(ggplot2)
library(bigsnpr)

obj.bed <- bed(bedfile <- download_1000G(dir = "../bigsnpr/tmp-data/"))
fam2 <- bigreadr::fread2(sub_bed(bedfile, ".fam2"))
THR_R2 <- 0.05

# computing PCA for Africans, to remove admixed/outlier individuals
obj.svd <- runonce::save_run(
  bed_autoSVD(obj.bed, ind.row = which(fam2$`Super Population` == "AFR"),
              k = 10, ncores = 4),
  file = "tmp-data/svdAFR.rds"
)
plot(obj.svd, type = "scores", scores = 1:8)
PC <- predict(obj.svd)
center <- bigutilsr::geometric_median(PC)
dist <- sqrt(bigutilsr:::rowSumsSq(sweep(PC, 2, center, '-')))
hist(dist, "FD")
mean(dist < 150) # 95%
ind.afr2 <- which(fam2$`Super Population` == "AFR")[dist < 150]

for (pop in c(unique(fam2$`Super Population`), "AFR2")) {

  print(pop)

  all_res <- lapply(1:22, function(CHR) {

    print(CHR)

    # prepare for snp_splitld
    ind.row <- `if`(pop == "AFR2", ind.afr2, which(fam2$`Super Population` == pop))
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
    ldetect_chr <- try(bigreadr::fread2(url), silent = TRUE)
    if (is.data.frame(ldetect_chr)) {
      n_block0 <- nrow(ldetect_chr)
      all_ind0 <- ldetect_chr$stop[-n_block0]
      block_num0 <- rowSums(outer(map$physical.pos, all_ind0, '>=')) + 1L
      stopifnot(length(block_num0) == ncol(corr))
      (cost0 <- bigsnpr:::compute_cost(block_num0, Matrix::tril(corr), THR_R2))
    } else {
      n_block0 <- ncol(corr) / 500
      cost0 <- NA_real_
    }

    # snp_splitld
    res <- runonce::save_run(
      bigsnpr::snp_ldsplit(corr, thr_r2 = THR_R2, min_size = 500, max_size = 10e3,
                           max_K = n_block0),
      file = paste0("tmp-data/split_", pop, "_chr", CHR, ".rds")
    )

    one_res <- slice_sample(res, n = 1)
    one_cost <- bigsnpr:::compute_cost(one_res$block_num[[1]],
                                       Matrix::tril(corr), THR_R2)
    stopifnot(all.equal(one_cost, one_res$cost, tolerance = 1e-5))

    bind_cols(
      bind_rows(
        bind_cols(res[1:2], method = "snp_ldsplit"),
        tibble(cost = cost0, n_block = n_block0, method = "ldetect")
      ),
      chr = CHR
    )
  })

  qplot(data = bind_rows(all_res), n_block, cost, color = method, size = method) +
    scale_size_manual(values = c(snp_ldsplit = 1, ldetect = 3)) +
    theme_bw(12) +
    scale_x_continuous(breaks = 0:20 * 10, minor_breaks = 0:100 * 2) +
    scale_y_log10() +
    facet_wrap(~ chr, nrow = 8, scales = "free_x") +
    theme(legend.position = "top") +
    labs(y = "Sum of squared correlations outside blocks (log-scale)", x = "Number of blocks")
  ggsave(paste0("figures/cost_split_", pop, ".pdf"), width = 13, height = 13)
}
