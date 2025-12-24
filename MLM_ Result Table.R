## ===============================================================
## QTL PYRAMIDING — FINAL SUMMARY TABLES (FORCED SIGNS + UNITS)
## ===============================================================

## 0) PACKAGES ---------------------------------------------------
library(openxlsx)
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(emmeans)
  library(tidyr)
  library(janitor)
  library(lmerTest)
})

## 1) PATHS ------------------------------------------------------
base_dir <- "C:/Users/adele.jamalzei/Desktop/Dissertation/QTL Pyramiding Project/2025 Data/Analysis/Analysis"
file_qtl  <- file.path(base_dir, "Master_QTL_Pyramiding_Simplified.xlsx")
out_base  <- file.path(base_dir, "MLM_Summary_Tables_Results")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

## 2) DATA LOADING & CLEANING ------------------------------------
df <- read_excel(file_qtl) %>%
  janitor::clean_names() %>%
  mutate(
    recurrent_parent = str_trim(as.character(recurrent_parent)),
    donor_parent = case_when(
      donor_parent %in% c("KB", "Kingbird") ~ "Kingbird",
      donor_parent %in% c("HB", "PI683503", "PI683501") ~ "PI683503",
      TRUE ~ str_trim(as.character(donor_parent))
    ),
    location = str_trim(as.character(location)),
    block    = factor(bloc)
  ) %>%
  mutate(
    introgression_status = case_when(
      qtn6b == "+" & elf == "+" & gw2 == "+" ~ "MT",
      qtn6b == "-" & elf == "-" & gw2 == "-" ~ "WT",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(introgression_status), !is.na(location)) %>%
  mutate(
    introgression_status = factor(introgression_status, levels = c("WT", "MT")),
    recurrent_parent     = factor(recurrent_parent, levels = c("Kelse", "Scarlet")),
    location             = factor(location)
  )

## 3) TRAIT DEFINITIONS & UNITS (## NEW ##) -----------------------
# This dictionary maps your variable names to professional names with units
trait_units <- c(
  "yield_kg_ha"        = "Yield (kg/ha)",
  "yield_g_p"          = "Yield (g/plant)",
  "tgw"                = "TGW (g)",
  "grain_number"       = "Grain Number",
  "grain_weight"       = "Grain Weight (g)",
  "tiller_number"      = "Tiller Number",
  "prod_tiller_number" = "Productive Tillers",
  "plant_ht"           = "Plant Height (cm)",
  "spike_length"       = "Spike Length (cm)",
  "spikelets_no"       = "Spikelets per Spike",
  "biomass_dry"        = "Dry Biomass (g)",
  "biomass_green"      = "Green Biomass (g)",
  "prot_12"            = "Protein 12% (mb)",
  "prot_ai"            = "Protein (as is)",
  "prot_db"            = "Protein (dry basis)",
  "hardness"           = "Grain Hardness",
  "twt_g_l"            = "Test Weight (g/L)",
  "test_wt"            = "Test Weight (lb/bu)",
  "moisture"           = "Moisture (%)",
  "h_date"             = "Heading Date (Julian)",
  "dtm"                = "Days to Maturity",
  "grain_length"       = "Grain Length (mm)",
  "grain_width"        = "Grain Width (mm)"
)

ALL_TRAITS <- names(trait_units)
traits_to_analyze <- intersect(ALL_TRAITS, names(df))

df <- df %>%
  mutate(across(all_of(traits_to_analyze), ~ as.numeric(gsub("[^0-9.\\-]", "", .))))

## 4) ANALYSIS LOOP ---------------------------------------------
ALL_LOCATION_RESULTS <- data.frame()

for (tr in traits_to_analyze) {
  df_tr <- df %>%
    select(location, recurrent_parent, donor_parent, introgression_status, block, all_of(tr)) %>%
    filter(is.finite(.data[[tr]])) %>%
    droplevels()
  
  if (nrow(df_tr) < 20) next
  
  fml <- as.formula(paste0(tr, " ~ (introgression_status * recurrent_parent * donor_parent) + (introgression_status * location) + (1 | location:block)"))
  model <- try(lmer(fml, data = df_tr, REML = TRUE), silent = TRUE)
  if (inherits(model, "try-error")) next
  
  emm_obj <- emmeans(model, ~ introgression_status * location * recurrent_parent * donor_parent)
  con_df <- as.data.frame(contrast(emm_obj, method = "pairwise", by = c("recurrent_parent", "donor_parent", "location"), adjust = "tukey")) %>%
    rename(p_val = p.value)
  
  trait_perf <- as.data.frame(emm_obj) %>%
    select(location, recurrent_parent, donor_parent, introgression_status, emmean) %>%
    pivot_wider(names_from = introgression_status, values_from = emmean) %>%
    mutate(
      Trait = tr,
      Percent_Diff = ((MT - WT) / abs(WT)) * 100
    ) %>%
    left_join(con_df, by = c("location", "recurrent_parent", "donor_parent"))
  
  ALL_LOCATION_RESULTS <- bind_rows(ALL_LOCATION_RESULTS, trait_perf)
}

## 5) EXPORT FINAL TABLES + RAW DATA VALIDATION (## UPDATED ##) ---
for (loc in levels(df$location)) {
  loc_results <- ALL_LOCATION_RESULTS %>% filter(location == loc)
  if (nrow(loc_results) == 0) next
  
  raw_means <- df %>%
    filter(location == loc) %>%
    group_by(recurrent_parent, donor_parent, introgression_status) %>%
    summarise(across(all_of(traits_to_analyze), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(cols = all_of(traits_to_analyze), names_to = "Trait", values_to = "Raw_Mean") %>%
    pivot_wider(names_from = introgression_status, values_from = Raw_Mean, names_prefix = "Raw_")
  
  validation_data <- loc_results %>%
    select(location, recurrent_parent, donor_parent, Trait, WT_MLM = WT, MT_MLM = MT, Percent_Diff, p_val) %>%
    left_join(raw_means, by = c("recurrent_parent", "donor_parent", "Trait")) %>%
    mutate(
      Trait_Name = trait_units[Trait], ## NEW ##
      Sig = case_when(p_val < 0.01 ~ "**", p_val < 0.05 ~ "*", p_val < 0.10 ~ "+", TRUE ~ ""),
      Value = paste0(if_else(Percent_Diff >= 0, "+", "-"), sprintf("%.2f%%", abs(Percent_Diff)), Sig)
    )
  
  final_summary <- validation_data %>%
    mutate(Parent_Combo = paste0(recurrent_parent, " | ", donor_parent)) %>%
    select(Parent_Combo, Trait_Name, Value) %>%
    pivot_wider(names_from = Trait_Name, values_from = Value)
  
  # Ensure the traits appear in the logical order from your trait_units list
  ordered_traits <- intersect(unname(trait_units), colnames(final_summary))
  write.xlsx(final_summary %>% select(Parent_Combo, all_of(ordered_traits)), 
             file.path(out_base, paste0("Final_Summary_", loc, ".xlsx")))
  
  # Save validation data with units as well
  write.xlsx(validation_data %>% select(-Trait) %>% rename(Trait = Trait_Name), 
             file.path(out_base, paste0("Validation_Means_", loc, ".xlsx")))
}

cat("\nDONE! Your tables now include full trait names and units.\n")