## ===============================================================
## QTL PYRAMIDING — UNIFIED LINEAR MIXED MODEL (MLM) — FIXED
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
  library(ggplot2)
  library(lmerTest)
})

## 1) PATHS ------------------------------------------------------
base_dir <- "C:/Users/adele.jamalzei/Desktop/Dissertation/QTL Pyramiding Project/2025 Data/Analysis/Analysis"
file_qtl  <- file.path(base_dir, "Master_QTL_Pyramiding_Simplified.xlsx")
mlm_base  <- file.path(base_dir, "MLM_Unified_Analysis")
dir.create(mlm_base, recursive = TRUE, showWarnings = FALSE)

## 2) CONSTANTS (COLORS + LABELS) -------------------------------
COL_WT <- "#00BFC4"   # Wild type
COL_MT <- "#F8766D"   # Mutant type


intro_labels <- c(
  WT = "Wild type (without pyramided QTL)",
  MT = "Mutant type (with pyramided QTL)"
)

## 3) NUMERIC CLEANUP -------------------------------------------
num_clean <- function(x) {
  x_chr <- gsub(",", "", as.character(x))
  x_chr <- gsub("[^0-9.\\-]", "", x_chr)
  out <- suppressWarnings(as.numeric(x_chr))
  out[is.infinite(out)] <- NA
  out
}

## 4) READ + CLEAN ----------------------------------------------
df <- read_excel(file_qtl) %>%
  janitor::clean_names() %>%
  mutate(
    recurrent_parent = str_trim(as.character(recurrent_parent)),
    donor_parent     = str_trim(as.character(donor_parent)),
    location         = str_trim(as.character(location)),
    qtn6b            = str_trim(as.character(qtn6b)),
    elf              = str_trim(as.character(elf)),
    gw2              = str_trim(as.character(gw2))
  ) %>%
  filter(!tolower(name1) %in% c("kelse","scarlet")) %>%
  filter(!is.na(recurrent_parent) & recurrent_parent != "") %>%
  filter(!is.na(donor_parent) & donor_parent != "") %>%
  filter(!is.na(location) & location != "")

## 5) FACTORS & INTROGRESSION STATUS ----------------------------
df <- df %>%
  mutate(
    introgression_status = case_when(
      qtn6b == "+" & elf == "+" & gw2 == "+" ~ "MT",
      qtn6b == "-" & elf == "-" & gw2 == "-" ~ "WT",
      TRUE ~ NA_character_
    ),
    introgression_status = factor(introgression_status, levels = c("WT","MT")),
    recurrent_parent = factor(recurrent_parent, levels = c("Kelse","Scarlet")),
    donor_parent = factor(case_when(
      donor_parent %in% c("KB","Kingbird")  ~ "Kingbird",
      donor_parent %in% c("HB","PI683503")  ~ "PI683503",
      TRUE ~ as.character(donor_parent)
    )),
    block    = factor(bloc),
    location = factor(location)
  ) %>%
  filter(!is.na(introgression_status))

## 6) TRAITS -----------------------------------------------------
## Traits present in ALL THREE locations (Pullman + Othello + Almota)
traits_common <- c(
  "tiller_number","prod_tiller_number",
  "h_date","plant_ht",
  "prot_12","prot_ai","prot_db",
  "moisture","hardness","twt_g_l","test_wt",
  "biomass_dry","biomass_green",
  "yield_g_p","yield_kg_ha",
  "spike_length","spikelets_no"
)

## Traits present in Pullman + Othello (NOT Almota)
traits_grain_po <- c(
  "grain_length","grain_width",
  "grain_number","grain_weight","tgw"
)

## Traits present in Pullman + Almota (NOT Othello)
traits_dtm_pa <- c("dtm")

ALL_TRAITS <- c(traits_common, traits_grain_po, traits_dtm_pa)

traits_to_analyze <- intersect(ALL_TRAITS, names(df))

## Convert only those traits that exist
df <- df %>% mutate(across(all_of(traits_to_analyze), num_clean))

## 7) GLOBAL CONTAINER ------------------------------------------
MLM_SUMMARY_TABLE <- data.frame()

## 8) MLM FUNCTION ----------------------------------------------
run_one_mlm <- function(df_full, trait, trait_dir) {
  
  df_trait <- df_full %>%
    select(block, introgression_status, recurrent_parent, donor_parent, location, all_of(trait)) %>%
    filter(is.finite(.data[[trait]])) %>%
    droplevels()
  
  n_rows  <- nrow(df_trait)
  n_loc   <- nlevels(df_trait$location)
  n_intro <- nlevels(df_trait$introgression_status)
  
  if (n_loc < 2 || n_intro < 2 || n_rows < 20) {
    writeLines(
      paste0(
        "SKIPPED (insufficient data)\n",
        "Trait: ", trait, "\n",
        "Rows (finite): ", n_rows, " (need >= 20)\n",
        "Locations: ", n_loc, " (need >= 2)\n",
        "Introgression levels: ", n_intro, " (need >= 2)\n"
      ),
      file.path(trait_dir, "SKIPPED.txt")
    )
    return(invisible(NULL))
  }
  
  fml_lmer <- as.formula(
    paste0(
      trait, 
      " ~ (introgression_status * recurrent_parent * donor_parent) + 
          (introgression_status * location) + 
          (1 | location:block)"
    )
  )
  
  model_fit <- tryCatch(
    lmer(fml_lmer, data = df_trait, REML = TRUE),
    error = function(e) {
      writeLines(conditionMessage(e), file.path(trait_dir, "MLM_ERROR.txt"))
      return(NULL)
    }
  )
  if (is.null(model_fit)) return(invisible(NULL))
  
  ## Type III ANOVA
  mlm_anova <- anova(model_fit, type = 3)
  
  tidy_tab <- broom::tidy(mlm_anova) %>%
    mutate(Sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ "ns"
    ))
  
  write.csv(tidy_tab, file.path(trait_dir, "MLM_ANOVA_Tidy.csv"), row.names = FALSE)
  
  ## Fix p-value extractor (your old one was broken)
  extract_p <- function(term_name) {
    r <- tidy_tab %>% filter(.data$term == term_name)
    if (nrow(r) >= 1) return(r$p.value[1])
    NA_real_
  }
  
  MLM_SUMMARY_TABLE <<- bind_rows(
    MLM_SUMMARY_TABLE,
    data.frame(
      Trait   = trait,
      P_Intro = extract_p("introgression_status"),
      P_GxE   = extract_p("introgression_status:location"),
      P_4Way  = extract_p("introgression_status:recurrent_parent:donor_parent:location")
    )
  )
  
  ## EMMEANS for plotting
  emm <- emmeans(model_fit, ~ introgression_status * donor_parent * location | recurrent_parent)
  plot_df <- as.data.frame(emm)
  
  ## Make sure introgression_status is WT/MT (keep levels)
  plot_df$introgression_status <- factor(
    as.character(plot_df$introgression_status),
    levels = c("WT","MT")
  )
  
  ## One plot per recurrent parent
  for (rp in levels(plot_df$recurrent_parent)) {
    
    df_rp <- plot_df %>% filter(recurrent_parent == rp)
    
    p <- ggplot(
      df_rp,
      aes(
        x = location,
        y = emmean,
        color = introgression_status,
        group = introgression_status
      )
    ) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15) +
      facet_grid(recurrent_parent ~ donor_parent) +
      scale_color_manual(
        values = c(WT = COL_WT, MT = COL_MT),
        labels = intro_labels,
        name   = "Introgression Status"
      ) +
      labs(
        title = paste("Introgression Effect Across Locations —", trait),
        x = "Environment (Location)",
        y = paste0("Estimated Marginal Mean (", trait, ")")
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        legend.box = "horizontal",
        strip.background.y = element_rect(fill = "grey85", color = "black"),
        strip.text.y = element_text(face = "bold", size = 13, angle = -90),
        strip.background.x = element_rect(fill = "grey90", color = "black"),
        strip.text.x = element_text(face = "bold", size = 12)
      )
    
    ggsave(
      file.path(trait_dir, paste0("MLM_GxE_", rp, ".png")),
      p, width = 11, height = 7, dpi = 320
    )
  }
  
  ## Diagnostics
  png(file.path(trait_dir, "MLM_Diagnostics.png"), 1200, 600)
  par(mfrow = c(1, 2))
  plot(model_fit, which = 1)
  qqnorm(resid(model_fit)); qqline(resid(model_fit))
  dev.off()
  
  invisible(TRUE)
}

## 9) RUN ALL TRAITS --------------------------------------------
cat("Starting Unified MLM...\n")

for (trait in traits_to_analyze) {
  message("Running trait: ", trait)
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

cat("DONE — Unified MLM analysis complete.\n")



#------------------------
#joined
#-----------------------



## ===============================================================
## 6.5) ENFORCE JOINT-ANALYSIS TRAIT VALIDITY (CRITICAL FIX)
## ===============================================================
## ===============================================================
## QTL PYRAMIDING — UNIFIED LINEAR MIXED MODEL (MLM)
## Model: trait ~ introgression_status * recurrent_parent * donor_parent * location 
##         + (1 | location:block)
## ===============================================================

# ---------------- PACKAGES -------------------------------------
library(ggplot2)
library(ggtext)
library(lmerTest)
library(emmeans)
library(dplyr)
library(broom)

# ---------------- MAIN FUNCTION (REVISED FOR JOINED PLOT) ----------------
run_one_mlm_joined <- function(df_full, trait, trait_dir) {
  
  df_trait <- df_full %>%
    select(block, introgression_status, recurrent_parent, donor_parent, location, all_of(trait)) %>%
    filter(is.finite(.data[[trait]])) %>%
    droplevels()
  
  # Basic data check
  if (nlevels(df_trait$location) < 2 || nlevels(df_trait$introgression_status) < 2 || nrow(df_trait) < 30) {
    return(NULL)
  }
  
  # MODEL: Stabilized formula
  fml_lmer <- as.formula(paste0(trait, " ~ (introgression_status * recurrent_parent * donor_parent) + (introgression_status * location) + (1 | location:block)"))
  
  model_fit <- try(lmer(fml_lmer, data = df_trait, REML = TRUE), silent = TRUE)
  if (inherits(model_fit, "try-error")) return(NULL)
  
  # EMMs for Plotting
  emm_all <- emmeans(model_fit, ~ introgression_status * location * recurrent_parent * donor_parent)
  plot_df <- as.data.frame(emm_all) %>%
    mutate(introgression_status = recode(introgression_status, 
                                         WT = "Wild type", 
                                         MT = "Mutant type"))
  
  # THE PLOT
  p <- ggplot(plot_df, aes(x = location, y = emmean, group = introgression_status, color = introgression_status)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
    facet_grid(recurrent_parent ~ donor_parent) +
    scale_color_manual(values = c("Wild type" = "#00BFC4", "Mutant type" = "#F8766D")) +
    labs(
      title = paste("Effect of QTL Introgression Across Environments —", trait),
      # THE UNIFIED LEGEND LINE: Bold, clean, and pipe-separated
      subtitle = paste0(
        "**Parent Status:** ",
        "<span style='color:", col_recurrent, ";'>**Recurrent Parent**</span> | ",
        "<span style='color:", col_donor, ";'>**Donor Parent**</span> | ",
        "**QTL Status:** ",
        "<span style='color:#00BFC4;'>**Wild Type**</span> | ",
        "<span style='color:#F8766D;'>**Mutant Type**</span>"
      ),
      x = "Environment (Location)", 
      y = paste0("Estimated Mean (", trait, ")")
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_markdown(hjust = 0.5, size = 11), 
      legend.position = "none", # Legend is now handled by the subtitle
      
      # Recurrent Parent (Rows) - BLUE STRIPS
      strip.background.y = element_rect(fill = col_recurrent),
      strip.text.y = element_text(face = "bold", color = "white"),
      
      # Donor Parent (Cols) - GREEN STRIPS
      strip.background.x = element_rect(fill = col_donor),
      strip.text.x = element_text(face = "bold", color = "white")
    )
  
  ggsave(file.path(trait_dir, paste0("MLM_GxE_", trait, "_Combined.png")), p, width = 13, height = 8, dpi = 300)
}
# ---------------------------------------------------------------
# EXECUTION: This part actually runs the function above
# ---------------------------------------------------------------
cat("Starting Joint Figure Plotting...\n")

for (trait in traits_to_analyze) {
  # Define the directory where the plots should go
  trait_dir <- file.path(mlm_base, trait)
  
  # Ensure the directory exists
  if (!dir.exists(trait_dir)) dir.create(trait_dir, recursive = TRUE)
  
  message("Generating Combined Plot for: ", trait)
  
  # CALL THE FUNCTION
  run_one_mlm_joined(df, trait, trait_dir)
}

cat("DONE — Check your trait folders for 'MLM_Joint_Plot' files.\n")