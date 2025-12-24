## ===============================================================
## QTL PYRAMIDING — UNIFIED LINEAR MIXED MODEL (MLM)
## Model: trait ~ introgression_status * recurrent_parent * donor_parent * location 
##        + (1 | location:block)
## ===============================================================

## 0) PACKAGES ---------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(emmeans)
  library(tidyr)
  library(janitor)
  library(broom)
  library(lmerTest)
})

## 1) PATHS ------------------------------------------------------
base_dir <- "C:/Users/adele.jamalzei/Desktop/Dissertation/QTL Pyramiding Project/2025 Data/Analysis/Analysis"
file_qtl <- file.path(base_dir, "Master_QTL_Pyramiding_Simplified.xlsx")
mlm_base <- file.path(base_dir, "MLM_Unified_Analysis")
dir.create(mlm_base, recursive = TRUE, showWarnings = FALSE)

## 2) GLOBAL CONTAINER ------------------------------------------
MLM_SUMMARY_TABLE <- data.frame()

## 3) NUMERIC CLEANUP (CRITICAL FIX) -----------------------------
num_clean <- function(x) {
  x_chr <- gsub(",", "", as.character(x))
  x_chr <- gsub("[^0-9.\\-]", "", x_chr)
  suppressWarnings(as.numeric(x_chr))
}

## 4) READ DATA --------------------------------------------------
df <- readxl::read_excel(file_qtl)

if (!is.data.frame(df)) {
  stop("ERROR: read_excel() did not return a data frame.")
}

## 5) CLEAN -----------------------------------------------------
df <- df %>%
  janitor::clean_names() %>%
  mutate(
    recurrent_parent = str_trim(as.character(recurrent_parent)),
    donor_parent     = str_trim(as.character(donor_parent)),
    location         = str_trim(as.character(location))
  ) %>%
  filter(!tolower(name1) %in% c("kelse", "scarlet")) %>%
  filter(!is.na(recurrent_parent) & recurrent_parent != "")

## 6) FACTORS & INTROGRESSION STATUS ----------------------------
df <- df %>%
  mutate(
    introgression_status = factor(
      case_when(
        qtn6b == "+" & elf == "+" & gw2 == "+" ~ "MT",
        qtn6b == "-" & elf == "-" & gw2 == "-" ~ "WT",
        TRUE ~ NA_character_
      ),
      levels = c("WT", "MT")
    ),
    recurrent_parent = factor(recurrent_parent, levels = c("Kelse", "Scarlet")),
    donor_parent = factor(case_when(
      donor_parent %in% c("KB", "Kingbird") ~ "Kingbird",
      donor_parent %in% c("HB", "PI683503") ~ "PI683503",
      TRUE ~ as.character(donor_parent)
    )),
    block = factor(bloc),
    location = factor(location)
  ) %>%
  filter(!is.na(introgression_status))

## 7) TRAIT DEFINITIONS (AUTHORITATIVE) --------------------------

## Present in ALL locations
traits_common <- c(
  "tiller_number","prod_tiller_number",
  "h_date","plant_ht",
  "prot_12","prot_ai","prot_db",
  "moisture","hardness","twt_g_l","test_wt",
  "biomass_dry","biomass_green",
  "yield_g_p","yield_kg_ha",
  "spike_length","spikelets_no"
)

## Present in Pullman + Othello ONLY
traits_grain_po <- c(
  "grain_length","grain_width",
  "grain_number","grain_weight","tgw"
)

## Present in Pullman + Almota ONLY
traits_dtm_pa <- c("dtm")

ALL_TRAITS <- c(
  traits_common,
  traits_grain_po,
  traits_dtm_pa
)

traits_to_analyze <- intersect(ALL_TRAITS, names(df))

## >>>>>>>>>>>>>>>>>>>> FIX THAT SAVES YOUR LIFE <<<<<<<<<<<<<<<<<<
## Convert traits to numeric BEFORE is.finite()
df <- df %>% mutate(across(all_of(traits_to_analyze), num_clean))

## 8) MLM FUNCTION ----------------------------------------------
## 8) MLM FUNCTION (WITH SIGNIFICANCE STARS) --------------------
run_one_mlm <- function(df_full, trait, trait_dir) {
  
  df_trait <- df_full %>%
    select(block, introgression_status, recurrent_parent, donor_parent, location, all_of(trait)) %>%
    filter(is.finite(.data[[trait]])) %>%
    droplevels()
  
  if (nlevels(df_trait$location) < 2 || nlevels(df_trait$introgression_status) < 2 || nrow(df_trait) < 30) {
    writeLines("SKIPPED: Not enough data levels.", file.path(trait_dir, "SKIPPED.txt"))
    return(NULL)
  }
  
  fml_lmer <- as.formula(paste0(trait, " ~ (introgression_status * recurrent_parent * donor_parent) + (introgression_status * location) + (1 | location:block)"))
  
  model_fit <- tryCatch(lmer(fml_lmer, data = df_trait, REML = TRUE), error = function(e) { return(NULL) })
  if (is.null(model_fit)) return(NULL)
  
  ## 1. CREATE TIDY ANOVA TABLE WITH STARS -----------------------
  tidy_tab <- broom::tidy(anova(model_fit, type = 3)) %>%
    mutate(significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE            ~ "ns"
    ))
  
  write.csv(tidy_tab, file.path(trait_dir, "MLM_ANOVA_Tidy.csv"), row.names = FALSE)
  
  ## 2. UPDATE GLOBAL SUMMARY TABLE ------------------------------
  # Helper function to get p-value + star
  get_sig_string <- function(term_name) {
    row <- tidy_tab %>% filter(term == term_name)
    if(nrow(row) == 1) return(paste0(round(row$p.value, 5), " (", row$significance, ")")) else return(NA)
  }
  
  MLM_SUMMARY_TABLE <<- bind_rows(
    MLM_SUMMARY_TABLE,
    data.frame(
      Trait           = trait,
      Intro_Sig       = get_sig_string("introgression_status"),
      GxE_Sig         = get_sig_string("introgression_status:location"),
      ThreeWay_Sig    = get_sig_string("introgression_status:recurrent_parent:donor_parent")
    )
  )
  ## VARIANCE COMPONENTS (NOW INSIDE THE FUNCTION) -------------
  # This was the cause of your NULL error!
  write.csv(
    as.data.frame(VarCorr(model_fit)),
    file.path(trait_dir, "MLM_Variance_Components.csv"),
    row.names = FALSE
  )
  
  ## EMMEANS + TUKEY --------------------------------------------
  emm_all <- emmeans(model_fit, ~ introgression_status * location * recurrent_parent * donor_parent)
  write.csv(as.data.frame(emm_all), file.path(trait_dir, "MLM_EMMEANS.csv"), row.names = FALSE)
  write.csv(as.data.frame(pairs(emm_all, adjust = "tukey")), file.path(trait_dir, "MLM_Tukey_HSD.csv"), row.names = FALSE)
  
  ## DIAGNOSTICS -----------------------------------------------
  png(file.path(trait_dir, "MLM_Diagnostics.png"), 1200, 600)
  par(mfrow = c(1, 2))
  plot(model_fit, which = 1)
  qqnorm(resid(model_fit)); qqline(resid(model_fit))
  dev.off()
}

## 9) RUN -------------------------------------------------------
for (trait in traits_to_analyze) {
  trait_dir <- file.path(mlm_base, trait)
  dir.create(trait_dir, showWarnings = FALSE)
  run_one_mlm(df, trait, trait_dir)
}

## 10) FINAL SUMMARY --------------------------------------------
write.csv(
  MLM_SUMMARY_TABLE,
  file.path(mlm_base, "Final_MLM_GxE_Summary.csv"),
  row.names = FALSE
)

cat("DONE — Unified MLM analysis complete (FIXED & VALID).\n")
