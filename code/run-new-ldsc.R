library(bigsnpr)
library(dplyr)

map_rds <- runonce::download_file(
  "https://www.dropbox.com/s/hdui60p9ohyhvv5/map_blocks.rds?dl=1",
  dir = "tmp-data")

pheno_files <- c(list.files("data/ukbb-quant-pheno", full.names = TRUE),
                 list.files("data/ukbb-binary-pheno", full.names = TRUE))

files <- tibble(
  pheno_file = pheno_files,
  gwas_file = file.path("GWAS", basename(pheno_file)),
  res_file = file.path("ldsc_blocks", basename(pheno_file))
) %>%
  filter(file.exists(gwas_file), !file.exists(res_file),
         !startsWith(basename(pheno_file), "LTFH_")) %>%
  print(n = Inf)

bigassertr::assert_dir("ldsc_blocks")

NCORES <- 16
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "20g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files, function(pheno_file, gwas_file, res_file) {

  # Get effective sample size from GWAS
  fam <- snp_attach("data/UKBB_HM3.rds")$fam
  ind_csv <- match(fam$eid, readRDS("data/csv_eid.rds"))
  y <- readRDS(print(pheno_file))[ind_csv]
  ind.train <- which(!is.na(y) & (fam$set == "train"))
  N <- length(ind.train)

  ## Information for the variants provided in the LD reference
  map_ldref <- readRDS(map_rds)
  map_ukbb <- snp_attach("data/UKBB_HM3.rds")$map %>%
    transmute(chr = as.integer(chromosome), pos = physical.pos,
              a0 = allele1, a1 = allele2)

  sumstats <- readRDS(gwas_file) %>%
    transmute(beta = estim, beta_se = std.err, n_eff = N) %>%
    bind_cols(map_ukbb)

  df_beta <- snp_match(sumstats, map_ldref) %>%
    as_tibble() %>%
    select(-starts_with("pos_")) %>%
    print()

  # LD score regression
  (ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                  chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff,
                                  blocks = 242,
                                  ncores = NCORES)))

  (ldsc2 <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                   chi2 = (beta / beta_se)^2,
                                   sample_size = n_eff,
                                   blocks = block,
                                   ncores = NCORES)))

  bad_groups <- sort(df_beta$block + sample(0:1, nrow(df_beta), replace = TRUE))
  (ldsc3 <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                   chi2 = (beta / beta_se)^2,
                                   sample_size = n_eff,
                                   blocks = bad_groups,
                                   ncores = NCORES)))

  saveRDS(rbind(normal = ldsc, good = ldsc2, bad = ldsc3),
          res_file)
})


library(bigsnpr)
library(dplyr)
all_res2 <-
  tibble(res_file = list.files("ldsc_blocks", full.names = TRUE),
         name = sub("\\.rds$", "", basename(res_file))) %>%
  mutate(res = lapply(res_file, function(file) reshape2::melt(readRDS(file)))) %>%
  select(-res_file) %>%
  tidyr::unnest(res) %>%
  tidyr::pivot_wider(names_from = c(Var1, Var2))

all_res2 %>%
  filter(good_h2 != bad_h2)

all_res2 %>%
  arrange(desc(good_h2)) %>%
  select(name, starts_with("good_")) %>%
  print()

all_res2 %>%
  arrange(desc(abs(log(good_h2_se / normal_h2_se)))) %>%
  select(name, ends_with("_h2_se")) %>%
  print()

library(ggplot2)
theme_set(theme_bw(10))
cowplot::plot_grid(

  qplot(bad_int, good_int, data = all_res2) +
    geom_abline(color = "red") + coord_equal() +
    labs(x = "LDSc regression intercept using non-independent blocks",
         y = "LDSc regression intercept using nearly-independent blocks"),

  qplot(bad_int_se, good_int_se, data = all_res2) +
    geom_abline(color = "red") + coord_equal() +
    scale_x_log10() + scale_y_log10() +
    labs(x = "SE of LDSc regression intercept using non-independent blocks",
         y = "SE of LDSc regression intercept using nearly-independent blocks"),

  qplot(bad_h2, good_h2, data = all_res2) +
    geom_abline(color = "red") + coord_equal() +
    labs(x = "LDSc regression heritability using non-independent blocks",
         y = "LDSc regression heritability using nearly-independent blocks"),

  qplot(bad_h2_se, good_h2_se, data = all_res2) +
    geom_abline(color = "red") + coord_equal() +
    scale_x_log10() + scale_y_log10() +
    labs(x = "SE of LDSc regression heritability using non-independent blocks",
         y = "SE of LDSc regression heritability using nearly-independent blocks"),

  scale = 0.95
)
# ggsave("figures/compare_new_ldsc.pdf", width = 10, height = 10)
