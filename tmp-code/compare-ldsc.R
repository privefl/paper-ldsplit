names <- c(
  "PASS_ADHD_Demontis2018.sumstats", "PASS_AdultOnsetAsthma_Ferreira2019.sumstats",
  "PASS_AgeFirstBirth.sumstats", "PASS_AgeOfInitiation_Liu2019.sumstats",
  "PASS_Alzheimer.sumstats", "PASS_AlzheimersMaternal_Marioni2018.sumstats",
  "PASS_AlzheimersPaternal_Marioni2018.sumstats", "PASS_AlzheimersProxy_Marioni2018.sumstats",
  "PASS_Alzheimers_Jansen2019.sumstats", "PASS_Anorexia.sumstats",
  "PASS_AnorexiaNervosa_Duncan_2017.sumstats", "PASS_AtrialFibrillation_Nielsen2018.sumstats",
  "PASS_Autism.sumstats", "PASS_Autism_Grove2019.sumstats", "PASS_BDSCZ_Ruderfer2018.sumstats",
  "PASS_BIP_Stahl2019.sumstats", "PASS_BMI1.sumstats", "PASS_BipolarDisorder_Ruderfer2018.sumstats",
  "PASS_Bipolar_Disorder.sumstats", "PASS_Bone_Mineral_Density_Kemp_2017.sumstats",
  "PASS_BreastCancer.sumstats", "PASS_CD_deLange2017.sumstats",
  "PASS_CardioembolicStroke_Malik2018.sumstats", "PASS_Celiac.sumstats",
  "PASS_ChildOnsetAsthma_Ferreira2019.sumstats", "PASS_CigarettesPerDay_Liu2019.sumstats",
  "PASS_Coronary_Artery_Disease.sumstats", "PASS_Coronary_Artery_Disease_Howson_2017.sumstats",
  "PASS_Crohns_Disease.sumstats", "PASS_DS.sumstats", "PASS_DepressedAffect_Nagel2018.sumstats",
  "PASS_Depression_Nagel2018.sumstats", "PASS_DrinksPerWeek_Liu2019.sumstats",
  "PASS_ENIGMA2_MeanPutamen.sumstats", "PASS_Epilepsy_Anney_2014.sumstats",
  "PASS_Ever_Smoked.sumstats", "PASS_FastingGlucose_Manning.sumstats",
  "PASS_FetalBirthWeight_Warrington2019.sumstats", "PASS_GeneralRiskTolerance_KarlssonLinner2019.sumstats",
  "PASS_HDL.sumstats", "PASS_HbA1C.sumstats", "PASS_Height1.sumstats",
  "PASS_IBD.sumstats", "PASS_IBD_deLange2017.sumstats", "PASS_Insomnia_Jansen2019.sumstats",
  "PASS_Intelligence_SavageJansen2018.sumstats", "PASS_IschemicStroke_Malik2018.sumstats",
  "PASS_LDL.sumstats", "PASS_LargeArteryStroke_Malik2018.sumstats",
  "PASS_Lean_Body_Mass_Zillikens_2017.sumstats", "PASS_LongSleepDuration_Dashti2019.sumstats",
  "PASS_Lupus.sumstats", "PASS_MDD_Howard2019.sumstats", "PASS_MDD_Wray2018.sumstats",
  "PASS_MaternalBirthWeight_Warrington2019.sumstats", "PASS_MedicationUse_Wu2019.sumstats",
  "PASS_Menarche2017.sumstats", "PASS_Multiple_sclerosis.sumstats",
  "PASS_Neuroticism.sumstats", "PASS_Neuroticism_Nagel2018.sumstats",
  "PASS_NumberChildrenEverBorn.sumstats", "PASS_OvarianCancer.sumstats",
  "PASS_Parkinsons_23andMe_Corces2020.sumstats", "PASS_Primary_biliary_cirrhosis.sumstats",
  "PASS_ProstateCancer.sumstats", "PASS_ReactionTime_Davies2018.sumstats",
  "PASS_Rheumatoid_Arthritis.sumstats", "PASS_SCZvsBD_Ruderfer2018.sumstats",
  "PASS_SWB.sumstats", "PASS_Schizophrenia.sumstats", "PASS_Schizophrenia_Pardinas2018.sumstats",
  "PASS_Schizophrenia_Ruderfer2018.sumstats", "PASS_ShortSleepDuration_Dashti2019.sumstats",
  "PASS_SleepDuration_Dashti2019.sumstats", "PASS_SmallVesselStroke_Malik2018.sumstats",
  "PASS_SmokingCessation_Liu2019.sumstats", "PASS_SmokingInitiation_Liu2019.sumstats",
  "PASS_Stroke_Malik2018.sumstats", "PASS_Triglycerides.sumstats",
  "PASS_Type_1_Diabetes.sumstats", "PASS_Type_2_Diabetes.sumstats",
  "PASS_UC_deLange2017.sumstats", "PASS_Ulcerative_Colitis.sumstats",
  "PASS_VerbalNumericReasoning_Davies2018.sumstats", "PASS_Worry_Nagel2018.sumstats",
  "PASS_Years_of_Education1.sumstats", "PASS_Years_of_Education2.sumstats"
)

NCORES <- 5
URL <- "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/all_sumstats/"
map_ldref <- readRDS("map_blocks.rds")
M <- nrow(map_ldref)

bigassertr::assert_dir("results_ldsc")

all_res <- lapply(names, function(.) {

  res_file <- paste0("results_ldsc/", print(sub("PASS_(.+)\\.sumstats", "\\1", .)), ".rds")

  print(runonce::save_run(file = res_file, {

    file <- runonce::download_file(paste0(URL, .), dir = "tmp-data/sumstats")
    sumstats <- na.omit(bigreadr::fread2(file, select = c("SNP", "N", "Z")))

    ind <- match(map_ldref$rsid, sumstats$SNP)
    nona <- which(!is.na(ind))

    ld   <- map_ldref$ld[nona]
    chi2 <- sumstats$Z[ind[nona]]^2
    N    <- sumstats$N[ind[nona]]

    ldsc <- bigsnpr::snp_ldsc(ld, M, chi2, N, ncores = NCORES,
                              blocks = 242)

    ldsc2 <- bigsnpr::snp_ldsc(ld, M, chi2, N, ncores = NCORES,
                               blocks = map_ldref$block[nona])

    bad_groups <- sort(map_ldref$block + sample(0:1, nrow(map_ldref), replace = TRUE))
    ldsc3 <- bigsnpr::snp_ldsc(ld, M, chi2, N, ncores = NCORES,
                               blocks = bad_groups[nona])

    rbind(normal = ldsc, good = ldsc2, bad = ldsc3)
  }))
})

library(magrittr)
all_res2 <-
  tibble::tibble(gwas = sub("PASS_(.+)\\.sumstats", "\\1", names),
                 res = lapply(all_res, reshape2::melt)) %>%
  tidyr::unnest(res) %>%
  tidyr::pivot_wider(names_from = c(Var1, Var2))

library(ggplot2)
theme_set(theme_bw(16))
cowplot::plot_grid(

  qplot(normal_int, good_int, data = all_res2) +
    geom_abline(color = "red"),

  qplot(normal_int_se, good_int_se, data = all_res2) +
    geom_abline(color = "red") +
    scale_x_log10() + scale_y_log10(),

  qplot(normal_h2, good_h2, data = all_res2) +
    geom_abline(color = "red"),

  qplot(normal_h2_se, good_h2_se, data = all_res2) +
    geom_abline(color = "red") +
    scale_x_log10() + scale_y_log10(),

  scale = 0.95
)

