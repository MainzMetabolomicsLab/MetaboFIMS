
{
library("ggplot2") 
library("dplyr")
library("titanic")
library("cowplot")
library("missForest")
library("matrixStats")
library("tidyverse")
library("KODAMA")
library("vegan")
library("MASS")
library("reshape2")
library("knitr")
library("ggfortify")
library("ggdendroplot")
library("pheatmap")
library("S4Vectors")
library("SummarizedExperiment")
library("pmp")
library("reshape2")
library("gridExtra")
library("mdatools")
library("ropls")
library("mixOmics")
library("ggrepel")
library("pheatmap")
library("RColorBrewer")
library("scales")
library("colorRamp2")
library("circlize")
library("grid")
library("readr")  
library("stringr")  
}

library(dplyr)
library(tools)

library(dplyr)
library(stringr)

# ============================================================
# --- Define File Paths Configuration
# ============================================================

if (.Platform$OS.type == "unix") {
  base_ms_path   <- "/Volumes/T7/Arbeit/FIA/MSKohorte_NEU/10ppm annotation/"
  metadata_path  <- "/Users/fabianschmitt/Seafile/Meine Bibliothek/FIA/EmDia_MS_Kohorte_FIMS_2025 batch and metadata/Metadata_combined_analysis.csv"
} else {
  base_ms_path   <- "D:/Arbeit/FIA/MSKohorte_NEU/10ppm annotation/"
  metadata_path  <- "C:/Users/fabia/Seafile/Meine Bibliothek/FIA/EmDia_MS_Kohorte_FIMS_2025 batch and metadata/Metadata_combined_analysis.csv"
}

modes <- list(
  Positive = file.path(base_ms_path, "CSV/feature_matrix_transposed_MSKohorte_with_metadata.csv"),
  Negative = file.path(base_ms_path, "CSV_NEG/feature_matrix_transposed_MSKohorte_with_metadata_neg.csv")
)

# ============================================================
# --- Loop through both Positive and Negative Modes
# ============================================================

for (mode_name in names(modes)) {
  
  MS_data_path <- modes[[mode_name]]
  message("\n--------------------------------------------------")
  message("🚀 Processing ", mode_name, " mode dataset...")
  message("--------------------------------------------------")
  
  # --- 1. Load MS data and clean metadata rows ---
  if (!file.exists(MS_data_path)) {
    warning("❌ File not found, skipping: ", MS_data_path)
    next
  }
  
  MS_data_raw <- read.csv2(MS_data_path, header = TRUE, sep = ",", dec = ".")
  
  existing_meta <- MS_data_raw[1:3, ]       
  MS_data       <- MS_data_raw[-c(1:3), ]    
  
  sample_cols <- colnames(MS_data)[-1]       
  
  MetabolomicsID <- str_extract(sample_cols, "(?<=_(Mainz|Bochum)_)(\\d+)$")
  
  sample_info <- data.frame(
    column_name = sample_cols,
    MetabolomicsID = MetabolomicsID,
    stringsAsFactors = FALSE
  )
  
  if (!file.exists(metadata_path)) {
    stop("❌ External metadata file missing at: ", metadata_path)
  }
  MS_metadata <- read.csv2(metadata_path, header = TRUE, sep = ",", dec = ".")
  
  MS_metadata$MetabolomicsID <- as.character(MS_metadata$MetabolomicsID)
  sample_info$MetabolomicsID <- as.character(sample_info$MetabolomicsID)
  
  sample_info <- sample_info %>%
    left_join(MS_metadata, by = "MetabolomicsID")
  
  sample_info <- sample_info[match(sample_cols, sample_info$column_name), ]
  
  sex_row      <- data.frame(SampleID = "Sex",      t(sample_info$Sex),      stringsAsFactors = FALSE)
  group_row    <- data.frame(SampleID = "Group",    t(sample_info$Group),    stringsAsFactors = FALSE)
  visit_row    <- data.frame(SampleID = "Visit",    t(sample_info$Visit),    stringsAsFactors = FALSE)
  PseudoID_row <- data.frame(SampleID = "PseudoID", t(sample_info$PseudoID), stringsAsFactors = FALSE)
  age_row      <- data.frame(SampleID = "Age",      t(sample_info$Age),      stringsAsFactors = FALSE)
  
  colnames(sex_row)       <- colnames(MS_data)
  colnames(group_row)     <- colnames(MS_data)
  colnames(visit_row)     <- colnames(MS_data)
  colnames(PseudoID_row)  <- colnames(MS_data)
  colnames(age_row)       <- colnames(MS_data)
  
  meta_keep <- c("class", "index", "batch")
  
  MS_data <- MS_data[!MS_data$SampleID %in% c(meta_keep, "Sex", "Group", "Visit", "PseudoID", "Age"), ]
  
  MS_data_new <- rbind(
    existing_meta,   
    sex_row,
    group_row,
    visit_row,
    PseudoID_row,
    age_row,         
    MS_data          
  )
  
  MS_data_new[is.na(MS_data_new)] <- ""
  
  unmatched <- sample_info %>% filter(!is.na(MetabolomicsID) & (is.na(Visit) | is.na(Group)))
  
  if (nrow(unmatched) > 0) {
    message("⚠️ Warning: Some Metabolomics IDs did not match metadata in ", mode_name, ":")
    print(unmatched$column_name)
  }
  
  row_headers <- MS_data_new$SampleID
  MS_data_transposed <- as.data.frame(t(MS_data_new[, -1]))
  colnames(MS_data_transposed) <- row_headers
  
  MS_data_transposed <- cbind(SampleID = rownames(MS_data_transposed), MS_data_transposed)
  rownames(MS_data_transposed) <- NULL
  
  MS_data_dir          <- dirname(MS_data_path)
  MS_data_filename     <- basename(MS_data_path)
  MS_data_new_filename <- sub("\\.csv$", "_with_age.csv", MS_data_filename)
  MS_data_new_path     <- file.path(MS_data_dir, MS_data_new_filename)
  
  write.csv(MS_data_transposed, file = MS_data_new_path, row.names = FALSE, quote = TRUE)
  
  message("✅ Success! Transposed file saved to: ", MS_data_new_path)
}






# ── File paths ────────────────────────────────────────────────────────────────
is_unix <- .Platform$OS.type == "unix"

file_paths <- list(
  pos = if (is_unix) {
    "/Volumes/T7/Arbeit/FIA/MSKohorte_NEU/10ppm annotation/CSV/feature_matrix_transposed_MSKohorte_with_metadata_with_age.csv"
  } else {
    "D:/Arbeit/FIA/MSKohorte_NEU/10ppm annotation/CSV/feature_matrix_transposed_MSKohorte_with_metadata_with_age.csv"
  },
  neg = if (is_unix) {
    "/Volumes/T7/Arbeit/FIA/MSKohorte_NEU/10ppm annotation/CSV_NEG/feature_matrix_transposed_MSKohorte_with_metadata_neg_with_age.csv"
  } else {
    "D:/Arbeit/FIA/MSKohorte_NEU/10ppm annotation/CSV_NEG/feature_matrix_transposed_MSKohorte_with_metadata_neg_with_age.csv"
  }
)

meta_cols <- c("class", "index", "batch", "Sex", "Group", "Visit", "PseudoID", "Age")

# ── Per-mode processing ───────────────────────────────────────────────────────
for (mode in names(file_paths)) {
  
  input_file <- file_paths[[mode]]
  message("\n", strrep("─", 60))
  message("▶  Processing ionisation mode: ", toupper(mode))
  message("   File: ", input_file)
  
  if (!file.exists(input_file)) {
    warning("File not found, skipping: ", input_file)
    next
  }
  
  # 1. Load ──────────────────────────────────────────────────────────────────
  dat <- read.csv(
    input_file,
    header      = TRUE,
    sep         = ",",
    dec         = ".",
    check.names = FALSE,   
    colClasses  = "character"
  )
  
  rownames(dat) <- dat[, 1]
  dat           <- dat[, -1, drop = FALSE]
  
  colnames(dat) <- trimws(colnames(dat))
  
  feature_cols <- setdiff(colnames(dat), meta_cols)
  message("   Features before filtering : ", length(feature_cols))
  
  # 2. QC CV filter (< 25 %) ─────────────────────────────────────────────────
  qc_data <- dat %>% filter(class == "QC")
  
  if (nrow(qc_data) == 0) {
    warning("No 'QC' samples found. Skipping mode: ", toupper(mode))
    next
  }
  
  qc_cvs <- sapply(feature_cols, function(f) {
    v <- as.numeric(qc_data[[f]])
    m <- mean(v, na.rm = TRUE)
    s <- sd(v,   na.rm = TRUE)
    if (is.na(m) || m == 0) return(NA_real_)
    (s / m) * 100
  })
  
  features_pass_cv <- names(qc_cvs)[!is.na(qc_cvs) & qc_cvs < 25]
  message("   Features passing QC CV < 25 %  : ", length(features_pass_cv))
  
  # 3. Blank filter (QC/Blank ≥ 4×) ─────────────────────────────────────────
  blank_data <- dat %>% filter(class == "Blank_processing")
  
  qc_means    <- sapply(feature_cols, function(f) mean(as.numeric(qc_data[[f]]),    na.rm = TRUE))
  blank_means <- sapply(feature_cols, function(f) mean(as.numeric(blank_data[[f]]), na.rm = TRUE))
  blank_means[blank_means == 0 | is.na(blank_means)] <- NA_real_
  
  ratio               <- qc_means / blank_means
  features_pass_blank <- names(ratio)[!is.na(ratio) & ratio >= 4]
  message("   Features passing QC/Blank ≥ 4× : ", length(features_pass_blank))
  
  # 4. Combine & subset ──────────────────────────────────────────────────────
  features_to_keep <- intersect(features_pass_cv, features_pass_blank)
  message("   Features retained (both filters): ", length(features_to_keep))
  
  valid_cols   <- intersect(c(meta_cols, features_to_keep), colnames(dat))
  dat_filtered <- dat[, valid_cols, drop = FALSE]
  
  # Replace NA / "NA" with empty string
  dat_filtered[is.na(dat_filtered) | dat_filtered == "NA"] <- ""
  
  message("   Final dimensions: ",
          nrow(dat_filtered), " samples × ",
          ncol(dat_filtered), " columns")
  
  # 5. Save ──────────────────────────────────────────────────────────────────
  out_file <- file.path(
    dirname(input_file),
    paste0(tools::file_path_sans_ext(basename(input_file)), "_CV_Blank_filtered.csv")
  )
  
  # write.table with quote = TRUE wraps all fields in double quotes,
  # so metabolite names containing commas (e.g. "1,2-Ethanediamine") are
  # preserved as a single field and not split by the CSV delimiter.
  write.table(
    cbind(SampleID = rownames(dat_filtered), dat_filtered),
    file      = out_file,
    sep       = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote     = TRUE
  )
  
  message("✅ Saved: ", out_file)
}

message("\n", strrep("─", 60))
message("Done — both modes processed.")



# ─────────────────────────────────────────────────────────────────────────────
# Load libraries
# ─────────────────────────────────────────────────────────────────────────────
library(tidyverse)
library(ggpubr)
library(patchwork)
library(lme4)
library(lmerTest)
library(pheatmap)
library(RColorBrewer)

# ─────────────────────────────────────────────────────────────────────────────
# File paths for positive and negative ion mode
# ─────────────────────────────────────────────────────────────────────────────
ion_modes <- list(
  POS = if (is_unix) {
    "/Volumes/T7/Arbeit/FIA/MSKohorte_NEU/10ppm annotation/CSV/feature_matrix_transposed_MSKohorte_with_metadata_with_age_CV_Blank_filtered.csv"
  } else {
    "D:/Arbeit/FIA/MSKohorte_NEU/10ppm annotation/CSV/feature_matrix_transposed_MSKohorte_with_metadata_with_age_CV_Blank_filtered.csv"
  },
  NEG = if (is_unix) {
    "/Volumes/T7/Arbeit/FIA/MSKohorte_NEU/10ppm annotation/CSV_NEG/feature_matrix_transposed_MSKohorte_with_metadata_neg_with_age_CV_Blank_filtered.csv"
  } else {
    "D:/Arbeit/FIA/MSKohorte_NEU/10ppm annotation/CSV_NEG/feature_matrix_transposed_MSKohorte_with_metadata_neg_with_age_CV_Blank_filtered.csv"
  }
)

# ─────────────────────────────────────────────────────────────────────────────
# Parameters (shared across both modes)
# ─────────────────────────────────────────────────────────────────────────────
meta_cols    <- c("class", "index", "batch", "Sex", "Group", "Visit", "PseudoID", "Age")
valid_visits <- c("T1", "T2", "T3", "T4")

# ─────────────────────────────────────────────────────────────────────────────
# Main loop
# ─────────────────────────────────────────────────────────────────────────────
for (mode in names(ion_modes)) {
  
  MSKohorte <- ion_modes[[mode]]
  out_dir   <- dirname(MSKohorte)
  
  cat("\n\n")
  cat("╔══════════════════════════════════════════════╗\n")
  cat(sprintf("  Processing: %s ion mode\n", mode))
  cat(sprintf("  File: %s\n", basename(MSKohorte)))
  cat("╚══════════════════════════════════════════════╝\n\n")
  
  # ── Load ──────────────────────────────────────────────────────────────────
  MSKohorte_data_transposed <- read.csv(
    MSKohorte,
    header      = TRUE,
    sep         = ",",
    dec         = ".",
    check.names = FALSE,
    colClasses  = "character"
  )
  
  rownames(MSKohorte_data_transposed) <- MSKohorte_data_transposed[, 1]
  MSKohorte_data_transposed           <- MSKohorte_data_transposed[, -1, drop = FALSE]
  MSKohorte_data_transposed[["class"]] <- as.character(MSKohorte_data_transposed[["class"]])
  
  feature_cols <- setdiff(colnames(MSKohorte_data_transposed), meta_cols)
  message("Loaded: ", nrow(MSKohorte_data_transposed), " samples × ", length(feature_cols), " features")
  
  # ── PCA of full dataset ───────────────────────────────────────────────────
  pca_features <- MSKohorte_data_transposed[-c(1:8)]
  pca_features <- apply(pca_features, 2, as.numeric)
  
  pca_res <- prcomp(pca_features, center = TRUE, scale. = TRUE)
  
  pca_df <- data.frame(
    MSKohorte_data_transposed[, intersect(meta_cols, colnames(MSKohorte_data_transposed)), drop = FALSE],
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2]
  )
  
  var_expl <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
  pc1_lab  <- paste0("PC1 (", var_expl[1], "%)")
  pc2_lab  <- paste0("PC2 (", var_expl[2], "%)")
  
  n_groups   <- length(unique(pca_df$class))
  brewer_cols <- colorRampPalette(brewer.pal(min(n_groups, 8), "Dark2"))(n_groups)
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = class)) +
    geom_point(shape = 21, size = 5, color = "black", alpha = 0.85, stroke = 0.3) +
    scale_fill_manual(values = brewer_cols,
                      guide = guide_legend(override.aes = list(shape = 22, size = 5))) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 1) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 1) +
    labs(title = paste0("PCA MSKohorte [ESI", mode, "]"), x = pc1_lab, y = pc2_lab) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid   = element_blank(),
      axis.line    = element_line(color = "black"),
      axis.ticks   = element_line(color = "black"),
      axis.title   = element_text(size = 30),
      axis.text    = element_text(size = 30),
      plot.title   = element_text(size = 30, face = "bold"),
      legend.title = element_blank(),
      legend.text  = element_text(size = 30)
    )
  
  print(p_pca)
  ggsave(file.path(out_dir, paste0("PCA_MSKohorte_full_data_", mode, ".png")),
         p_pca, width = 9, height = 7, dpi = 400)
  
  # ── CV plot ───────────────────────────────────────────────────────────────
  MSKohorte_data_transposed_CV <- MSKohorte_data_transposed %>%
    filter(xor(class == "QC", class == "Sample")) %>%
    dplyr::select(-c(index, Sex, Group, Visit, PseudoID, Age)) %>%
    rownames_to_column(var = "SampleID")
  
  data_long_cv <- pivot_longer(MSKohorte_data_transposed_CV,
                               cols = -c(SampleID, batch, class),
                               names_to = "Metabolite", values_to = "Value") %>%
    mutate(Value = as.numeric(Value))
  
  cv_data <- data_long_cv %>%
    group_by(batch, class, Metabolite) %>%
    summarise(
      Mean_Value = mean(Value, na.rm = TRUE),
      SD_Value   = sd(Value,   na.rm = TRUE),
      CV         = (SD_Value / Mean_Value) * 100,
      .groups    = "drop"
    ) %>%
    filter(Mean_Value != 0)
  
  p_cv <- ggplot(cv_data, aes(x = class, y = CV, fill = batch)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 position = position_dodge(width = 0.9), color = "black") +
    geom_hline(yintercept = 25, linetype = "dotted", color = "red", linewidth = 1) +
    ggplot2::annotate("text", x = 1.5, y = 10, label = "25% threshold",
                      color = "red", size = 5, hjust = 0) +
    coord_cartesian(ylim = c(0, 200)) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title    = paste0("CV of Metabolites by Batch and Class (", mode, ")"),
      subtitle = "Violin and boxplot (CV capped at 200%)",
      x = "Class", y = "Coefficient of Variation (%)", fill = "Batch"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x  = element_text(size = 16, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 16),
      axis.title   = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
      plot.title   = element_text(size = 20, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 16, hjust = 0.5)
    )
  
  print(p_cv)
  ggsave(file.path(out_dir, paste0("CVs_MSKohorte_", mode, ".png")),
         p_cv, width = 9, height = 7, dpi = 400)
  
  # ── Prepare long-format data for LMM ─────────────────────────────────────
  # ── Prepare long-format data for LMM ─────────────────────────────────────
  MS <- as.data.frame(MSKohorte_data_transposed, stringsAsFactors = FALSE) %>%
    rownames_to_column("Sample") %>%
    as_tibble() %>%
    filter(Visit %in% valid_visits & !is.na(Visit))
  
  lmm_meta_cols <- c("class", "index", "batch", "Sex", "Group", "Visit", "PseudoID")
  metadata_available <- intersect(c("Sample", lmm_meta_cols), colnames(MS))
  metadata        <- MS[, metadata_available]
  metabolite_cols <- setdiff(colnames(MS), colnames(metadata))
  
  data_long <- MS[, metabolite_cols] %>%
    mutate(Sample = MS$Sample) %>%
    pivot_longer(cols = -Sample, names_to = "Metabolite", values_to = "Intensity") %>%
    left_join(metadata, by = "Sample") %>%
    dplyr::filter(class == "Sample") %>%
    dplyr::select(PseudoID, Metabolite, Visit, Intensity) %>%
    mutate(
      # Crucial Fix: Force the main Intensity column to be numeric
      # so your original downstream plotting loops do not crash
      Intensity     = as.numeric(Intensity),
      # Create the log-transformed column strictly for the LMM analysis
      Intensity_Log = log2(Intensity + 1), 
      Visit         = factor(Visit, levels = valid_visits, ordered = TRUE)
    ) %>%
    filter(!is.na(PseudoID) & PseudoID != "" & !is.na(Intensity))
  
  # ── Linear Mixed Models (Run on transformed data) ────────────────────────
  cat("Running Linear Mixed Models on log2-transformed data...\n")
  
  lmm_results <- data_long %>%
    group_by(Metabolite) %>%
    summarise(
      n_obs      = n(),
      n_patients = n_distinct(PseudoID),
      n_visits   = n_distinct(Visit),
      p_value    = {
        if (n_obs >= 20 && n_patients >= 5 && n_visits >= 2) {
          tryCatch({
            # Model fits on the stabilized log2 scale
            fit <- lmer(Intensity_Log ~ Visit + (1 | PseudoID), data = pick(everything()))
            anova(fit)["Visit", "Pr(>F)"]
          }, error = function(e) NA_real_)
        } else NA_real_
      },
      .groups = "drop"
    ) %>%
    filter(!is.na(p_value)) %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
  
  cat("Total metabolites tested:", nrow(lmm_results), "\n")
  cat("Significant at FDR < 0.05:", sum(lmm_results$FDR < 0.05, na.rm = TRUE), "\n")
  cat("Significant at FDR < 0.10:", sum(lmm_results$FDR < 0.10, na.rm = TRUE), "\n")
  print(head(lmm_results, 20))
  
  write_csv(lmm_results, file.path(out_dir, paste0("LMM_Visit_Results_", mode, ".csv")))
  
  sig_mets <- lmm_results %>% filter(FDR < 0.05) %>% pull(Metabolite)
  
  # ── Plots ─────────────────────────────────────────────────────────────────
  # ── Helper: fixed-format scientific notation for axis labels ────────────────
# Scientific notation formatted with 0 decimal places (e.g., 4E+06) — matches dual-axis version
format_scientific <- function(x) {
  ifelse(is.na(x), "", ifelse(x == 0, "0", sprintf("%.0E", x)))
}

  
  if (length(sig_mets) > 0) {
    
    sig_mets_top15 <- lmm_results %>%
      filter(FDR < 0.05) %>% arrange(p_value) %>% slice_head(n = 15) %>% pull(Metabolite)
    
    cat("\nGenerating plots for top", length(sig_mets_top15), "significant metabolites...\n")
    
    summary_stats <- data_long %>%
      filter(Metabolite %in% sig_mets) %>%
      group_by(Metabolite, Visit) %>%
      summarise(
        mean_int = mean(Intensity, na.rm = TRUE),
        se       = sd(Intensity, na.rm = TRUE) / sqrt(n()),
        .groups  = "drop"
      )
    
    n_sig       <- length(sig_mets_top15)
    ncol_plot   <- min(3, n_sig)
    height_plot <- max(3, min(9, ceiling(n_sig / ncol_plot) * 1.5))
    
    # Boxplots
    p_box <- ggplot(data_long %>% filter(Metabolite %in% sig_mets_top15),
                    aes(x = Visit, y = Intensity)) +
      geom_boxplot(aes(fill = Visit), outlier.size = 1, outlier.alpha = 0.5, linewidth = 0.5) +
      facet_wrap(~Metabolite, scales = "free_y", ncol = ncol_plot) +
      scale_y_continuous(n.breaks = 4, labels = format_scientific)  +
      labs(y = "Relative intensity", x = "Visit") +
      theme_bw(base_size = 11) +
      theme(
        strip.text       = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "grey60"),
        panel.grid.major = element_line(linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position  = "none",
        axis.text        = element_text(size = 8),
        axis.title       = element_text(size = 8, face = "bold"),
        plot.margin      = margin(8, 8, 8, 8)
      )
    
    print(p_box)
    ggsave(file.path(out_dir, paste0("LMM_Significant_Metabolites_Boxplots_", mode, ".png")),
           p_box, width = 7, height = height_plot, units = "in", dpi = 600)
    ggsave(file.path(out_dir, paste0("LMM_Significant_Metabolites_Boxplots_", mode, ".pdf")),
           p_box, width = 7, height = height_plot, units = "in")
    
    
    # ── Individual spaghetti plots (all significant metabolites) ────────────────
    spaghetti_dir <- file.path(out_dir, paste0("spaghetti_", mode))
    if (!dir.exists(spaghetti_dir)) dir.create(spaghetti_dir)
    
    wrap_text <- function(text, width = 30) {
      strwrap(text, width = width, simplify = FALSE)[[1]] %>% paste(collapse = "\n")
    }
    
    cat("\nExporting", length(sig_mets), "individual spaghetti plots (matching dual-axis dimensions)...\n")
    
    for (met in sig_mets) {
      met_data    <- data_long     %>% filter(Metabolite == met)
      met_summary <- summary_stats %>% filter(Metabolite == met)
      
      p_individual <- ggplot(met_data, aes(x = Visit, y = Intensity)) +
        geom_line(aes(group = PseudoID), alpha = 0.15, linewidth = 0.35, color = "#1B9E77") +
        geom_ribbon(data = met_summary,
                    aes(x = Visit, y = mean_int,
                        ymin = mean_int - se, ymax = mean_int + se, group = 1),
                    fill = "#1B9E77", alpha = 0.3, inherit.aes = FALSE) +
        geom_line(data = met_summary,
                  aes(x = Visit, y = mean_int, group = 1),
                  color = "#1B9E77", linewidth = 1.2, inherit.aes = FALSE) +
        geom_point(data = met_summary,
                   aes(x = Visit, y = mean_int, group = 1),
                   color = "#1B9E77", size = 3, shape = 21, fill = "white",
                   stroke = 1.3, inherit.aes = FALSE) +
        scale_y_continuous(
          name     = "FI-MS Intensity",
          labels   = format_scientific,
          sec.axis = sec_axis(~ ., name = "LC-MS Intensity", labels = format_scientific),
          n.breaks = 4
        ) +
        labs(title = wrap_text(met), x = "Visit") +
        theme_bw(base_size = 9) +
        theme(
          axis.title.x       = element_text(size = 9, face = "bold", margin = margin(t = 2)),
          axis.title.y.left  = element_text(color = "#1B9E77", face = "bold", size = 9, margin = margin(r = 2)),
          axis.text.y.left   = element_text(color = "#1B9E77", size = 8),
          axis.title.y.right = element_text(color = "white", face = "bold", size = 9, margin = margin(l = 2)),
          axis.text.y.right  = element_text(color = "white", size = 8),
          axis.text.x        = element_text(size = 8, color = "grey20"),
          axis.ticks         = element_line(linewidth = 0.25),
          plot.title         = element_text(hjust = 0.5, face = "bold", size = 8, lineheight = 0.85, margin = margin(b = 3)),
          panel.grid.major   = element_line(linewidth = 0.2, color = "grey88"),
          panel.grid.minor   = element_blank(),
          panel.border       = element_rect(color = "grey60", linewidth = 0.4),
          plot.margin        = margin(3, 3, 3, 3)
        )
      
      safe_name <- gsub("[/\\\\:*?\"<>|]", "_", met)
      ggsave(file.path(spaghetti_dir, paste0(safe_name, "_spaghetti.png")),
             p_individual, width = 2.3, height = 2, units = "in", dpi = 600)
    }
    
    cat("Exported", length(sig_mets), "individual plots to:", spaghetti_dir, "\n")
    
    # ── ALL Significant Metabolites Faceted (Mean) ───────────────────────────
    cat("\nGenerating combined high-DPI plot for all", length(sig_mets), "metabolites (5 columns)...\n")
    
    n_all_sig  <- length(sig_mets)
    all_rows   <- ceiling(n_all_sig / 5)
    all_height <- max(4, all_rows * 2)
    
    p_spaghetti_all <- ggplot(data_long %>% filter(Metabolite %in% sig_mets), 
                              aes(x = Visit, y = Intensity)) +
      geom_line(aes(group = PseudoID), alpha = 0.15, linewidth = 0.35, color = "#1B9E77") +
      geom_ribbon(data = summary_stats %>% filter(Metabolite %in% sig_mets),
                  aes(x = Visit, y = mean_int,
                      ymin = mean_int - se, ymax = mean_int + se, group = 1),
                  fill = "#1B9E77", alpha = 0.3, inherit.aes = FALSE) +
      geom_line(data = summary_stats %>% filter(Metabolite %in% sig_mets),
                aes(x = Visit, y = mean_int, group = 1),
                color = "#1B9E77", linewidth = 1.2, inherit.aes = FALSE) +
      geom_point(data = summary_stats %>% filter(Metabolite %in% sig_mets),
                 aes(x = Visit, y = mean_int, group = 1),
                 color = "#1B9E77", size = 2.5, shape = 21, fill = "white",
                 stroke = 1.2, inherit.aes = FALSE) +
      facet_wrap(~Metabolite, scales = "free_y", ncol = 5) + 
      scale_y_continuous(
        name     = "FI-MS Intensity",
        labels   = format_scientific,
        sec.axis = sec_axis(~ ., name = "LC-MS Intensity", labels = format_scientific),
        n.breaks = 4
      ) +
      labs(title = paste0("All Significant Metabolites [", mode, "]"), x = "Visit") +
      theme_bw(base_size = 9) +
      theme(
        strip.text         = element_text(size = 8, face = "bold", color = "grey20"),
        strip.background   = element_rect(fill = "grey92", color = "grey60", linewidth = 0.5),
        panel.grid.major   = element_line(linewidth = 0.25, color = "grey85"),
        panel.grid.minor   = element_blank(),
        panel.border       = element_rect(color = "grey60", linewidth = 0.5),
        legend.position    = "none",
        axis.text.x        = element_text(size = 8, color = "grey20"),
        axis.text.y.left   = element_text(size = 8, color = "#1B9E77"),
        axis.title.y.left  = element_text(size = 9, face = "bold", color = "#1B9E77", margin = margin(r = 2)),
        axis.title.y.right = element_text(size = 9, face = "bold", color = "white", margin = margin(l = 2)),
        axis.text.y.right  = element_text(size = 8, color = "white"),
        axis.title.x       = element_text(size = 9, face = "bold"),
        axis.ticks         = element_line(linewidth = 0.4),
        plot.margin        = margin(8, 8, 8, 8)
      )
    
    ggsave(file.path(spaghetti_dir, paste0("ALL_Significant_Spaghetti_Mean_", mode, ".png")),
           p_spaghetti_all, width = 12, height = all_height, units = "in", dpi = 600)
    
    ggsave(file.path(spaghetti_dir, paste0("ALL_Significant_Spaghetti_Mean_", mode, ".pdf")),
           p_spaghetti_all, width = 12, height = all_height, units = "in")
    
    # ── Individual spaghetti plots, MEDIAN version (extra output) ───────────
    summary_stats_median <- data_long %>%
      filter(Metabolite %in% sig_mets) %>%
      group_by(Metabolite, Visit) %>%
      summarise(
        median_int = median(Intensity, na.rm = TRUE),
        q25        = quantile(Intensity, 0.25, na.rm = TRUE),
        q75        = quantile(Intensity, 0.75, na.rm = TRUE),
        .groups    = "drop"
      )
    
    spaghetti_median_dir <- file.path(out_dir, paste0("spaghetti_", mode, "_median"))
    if (!dir.exists(spaghetti_median_dir)) dir.create(spaghetti_median_dir)
    
    cat("\nExporting", length(sig_mets), "individual spaghetti plots (median, matching dual layout)...\n")
    
    for (met in sig_mets) {
      met_data    <- data_long           %>% filter(Metabolite == met)
      met_summary <- summary_stats_median   %>% filter(Metabolite == met)
      
      p_individual_median <- ggplot(met_data, aes(x = Visit, y = Intensity)) +
        geom_line(aes(group = PseudoID), alpha = 0.2, linewidth = 0.4, color = "grey60") +
        geom_ribbon(data = met_summary,
                    aes(x = Visit, y = median_int,
                        ymin = q25, ymax = q75, group = 1),
                    fill = "#1B9E77", alpha = 0.3, inherit.aes = FALSE) +
        geom_line(data = met_summary,
                  aes(x = Visit, y = median_int, group = 1),
                  color = "#1B9E77", linewidth = 1.2, inherit.aes = FALSE) +
        geom_point(data = met_summary,
                   aes(x = Visit, y = median_int, group = 1),
                   color = "#1B9E77", size = 3, shape = 21, fill = "white",
                   stroke = 1.3, inherit.aes = FALSE) +
        scale_y_continuous(
          name     = "FI-MS Intensity",
          labels   = format_scientific,
          sec.axis = sec_axis(~ ., name = "LC-MS Intensity", labels = format_scientific),
          n.breaks = 4
        ) +
        labs(title = wrap_text(met), x = "Visit") +
        theme_bw(base_size = 9) +
        theme(
          axis.title.x       = element_text(size = 9, face = "bold", margin = margin(t = 2)),
          axis.title.y.left  = element_text(color = "#1B9E77", face = "bold", size = 9, margin = margin(r = 2)),
          axis.text.y.left   = element_text(color = "#1B9E77", size = 8),
          axis.title.y.right = element_text(color = "white", face = "bold", size = 9, margin = margin(l = 2)),
          axis.text.y.right  = element_text(color = "white", size = 8),
          axis.text.x        = element_text(size = 8, color = "grey20"),
          axis.ticks         = element_line(linewidth = 0.25),
          plot.title         = element_text(hjust = 0.5, face = "bold", size = 8, lineheight = 0.85, margin = margin(b = 3)),
          panel.grid.major   = element_line(linewidth = 0.2, color = "grey88"),
          panel.grid.minor   = element_blank(),
          panel.border       = element_rect(color = "grey60", linewidth = 0.4),
          plot.margin        = margin(3, 3, 3, 3)
        )
      
      safe_name <- gsub("[/\\\\:*?\"<>|]", "_", met)
      ggsave(file.path(spaghetti_median_dir, paste0(safe_name, "_spaghetti_median.png")),
             p_individual_median, width = 2.3, height = 2, units = "in", dpi = 600)
    }
    
    cat("Exported", length(sig_mets), "individual median plots to:", spaghetti_median_dir, "\n")
    
    # ── Trajectory clustering heatmap ────────────────────────────────────────
    trajectory_data <- data_long %>%
      filter(Metabolite %in% sig_mets_top15) %>%
      group_by(Metabolite, Visit) %>%
      summarise(mean_intensity = mean(Intensity, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Visit, values_from = mean_intensity) %>%
      column_to_rownames("Metabolite")
    
    trajectory_scaled <- t(scale(t(trajectory_data)))
    hc <- hclust(dist(trajectory_scaled), method = "ward.D2")
    
    p_heatmap <- pheatmap(
      trajectory_scaled,
      cluster_rows     = TRUE,
      cluster_cols     = FALSE,
      clustering_method = "ward.D2",
      color            = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      breaks           = seq(-2, 2, length.out = 101),
      main             = paste0("Metabolite Trajectory Patterns [", mode, "]"),
      fontsize         = 8, fontsize_row = 8, fontsize_col = 8,
      angle_col        = 0, border_color = "grey80",
      cellwidth = 15, cellheight = 15,
      legend_breaks  = c(-2, -1, 0, 1, 2),
      legend_labels  = c("-2 SD", "-1 SD", "Mean", "+1 SD", "+2 SD"),
      silent = TRUE
    )
    
    ggsave(file.path(out_dir, paste0("LMM_Trajectory_Heatmap_", mode, ".png")),
           p_heatmap, width = 5, height = max(6, nrow(trajectory_data) * 0.3),
           units = "in", dpi = 600)
    ggsave(file.path(out_dir, paste0("LMM_Trajectory_Heatmap_", mode, ".pdf")),
           p_heatmap, width = 5, height = max(6, nrow(trajectory_data) * 0.3),
           units = "in")
    
    # Trajectory cluster line plot
    n_clusters  <- min(3, max(2, ceiling(n_sig / 5)))
    clusters    <- cutree(hc, k = n_clusters)
    
    trajectory_clustered <- trajectory_data %>%
      rownames_to_column("Metabolite") %>%
      mutate(Cluster = paste0("Pattern ", clusters[Metabolite])) %>%
      pivot_longer(cols = any_of(valid_visits), names_to = "Visit", values_to = "Mean_Intensity") %>%
      mutate(Visit = factor(Visit, levels = valid_visits, ordered = TRUE))
    
    p_clusters <- ggplot(trajectory_clustered,
                         aes(x = Visit, y = Mean_Intensity, group = Metabolite, color = Cluster)) +
      geom_line(linewidth = 1, alpha = 0.8) +
      geom_point(size = 2, alpha = 0.8) +
      facet_wrap(~Cluster, scales = "free_y", ncol = n_clusters) +
      scale_color_brewer(palette = "Set1") +
      labs(y = "Mean relative intensity", x = "Visit",
           title = paste0("Distinct Temporal Patterns [", mode, "]")) +
      theme_bw(base_size = 11) +
      theme(
        strip.text       = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "grey60"),
        panel.grid.major = element_line(linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position  = "none",
        axis.text        = element_text(size = 9),
        axis.title       = element_text(size = 10, face = "bold"),
        plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.margin      = margin(8, 8, 8, 8)
      )
    
    print(p_clusters)
    ggsave(file.path(out_dir, paste0("LMM_Trajectory_Clusters_", mode, ".png")),
           p_clusters, width = 7, height = 3, units = "in", dpi = 600)
    ggsave(file.path(out_dir, paste0("LMM_Trajectory_Clusters_", mode, ".pdf")),
           p_clusters, width = 7, height = 3, units = "in")
    
    cluster_summary <- trajectory_clustered %>%
      dplyr::select(Metabolite, Cluster) %>% distinct() %>% arrange(Cluster, Metabolite)
    write_csv(cluster_summary, file.path(out_dir, paste0("LMM_Trajectory_Clusters_", mode, ".csv")))
    
    # ── PCA with patient trajectories ────────────────────────────────────────
    MS_pca <- as.data.frame(MSKohorte_data_transposed, stringsAsFactors = FALSE) %>%
      rownames_to_column("Sample") %>% as_tibble() %>%
      filter(Visit %in% valid_visits & !is.na(Visit),
             class == "Sample",
             !is.na(PseudoID) & PseudoID != "")
    
    pca_meta_cols  <- c("Sample", "class", "index", "batch", "Sex", "Group", "Visit", "PseudoID", "Age")
    pca_met_cols   <- setdiff(colnames(MS_pca), pca_meta_cols)
    pca_matrix     <- MS_pca %>% dplyr::select(all_of(pca_met_cols)) %>%
      mutate(across(everything(), as.numeric)) %>% as.matrix()
    pca_meta_sub   <- MS_pca %>% dplyr::select(PseudoID, Visit)
    
    complete_rows  <- complete.cases(pca_matrix)
    pca_matrix     <- pca_matrix[complete_rows, ]
    pca_meta_sub   <- pca_meta_sub[complete_rows, ]
    valid_cols     <- apply(pca_matrix, 2, function(x) !is.na(sd(x, na.rm=TRUE)) && sd(x, na.rm=TRUE) > 0)
    pca_matrix     <- pca_matrix[, valid_cols]
    
    pca_result     <- prcomp(pca_matrix, scale. = TRUE, center = TRUE)
    var_explained  <- summary(pca_result)$importance[2, 1:2] * 100
    pc1_lab_lmm    <- sprintf("PC1 (%.1f%%)", var_explained[1])
    pc2_lab_lmm    <- sprintf("PC2 (%.1f%%)", var_explained[2])
    
    pca_traj_df <- as.data.frame(pca_result$x[, 1:2]) %>%
      mutate(PseudoID = pca_meta_sub$PseudoID,
             Visit    = factor(pca_meta_sub$Visit, levels = valid_visits))
    
    n_v          <- length(valid_visits)
    visit_cols   <- colorRampPalette(brewer.pal(min(n_v, 8), "Dark2"))(n_v)
    
    pca_mean <- pca_traj_df %>%
      group_by(Visit) %>%
      summarise(PC1_mean = mean(PC1), PC2_mean = mean(PC2),
                PC1_se   = sd(PC1) / sqrt(n()),
                PC2_se   = sd(PC2) / sqrt(n()), .groups = "drop")
    
    p_pca_traj <- ggplot(pca_traj_df, aes(x = PC1, y = PC2, fill = Visit)) +
      geom_path(aes(group = PseudoID), color = "grey60", linewidth = 0.5, alpha = 0.4,
                arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
      geom_point(shape = 21, size = 5, color = "black", alpha = 0.85, stroke = 0.3) +
      scale_fill_manual(values = visit_cols,
                        guide = guide_legend(override.aes = list(shape = 22, size = 5))) +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 1) +
      geom_vline(xintercept = 0, color = "grey50", linewidth = 1) +
      labs(title = paste0("PCA: Patient Trajectories [", mode, "]"),
           x = pc1_lab_lmm, y = pc2_lab_lmm) +
      theme_minimal(base_size = 16) +
      theme(panel.grid = element_blank(), axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"), axis.title = element_text(size = 30),
            axis.text = element_text(size = 30), plot.title = element_text(size = 30, face = "bold"),
            legend.title = element_blank(), legend.text = element_text(size = 30))
    
    print(p_pca_traj)
    ggsave(file.path(out_dir, paste0("LMM_PCA_Trajectories_", mode, ".png")),
           p_pca_traj, width = 12, height = 9, units = "in", dpi = 600)
    
    p_pca_mean <- ggplot(pca_traj_df, aes(x = PC1, y = PC2, fill = Visit)) +
      geom_point(shape = 21, size = 3, color = "black", alpha = 0.3, stroke = 0.2) +
      geom_path(data = pca_mean, aes(x = PC1_mean, y = PC2_mean, group = 1),
                color = "black", linewidth = 1.5,
                arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
      geom_point(data = pca_mean, aes(x = PC1_mean, y = PC2_mean, fill = Visit),
                 shape = 21, size = 8, color = "black", stroke = 0.5) +
      geom_errorbar(data = pca_mean,
                    aes(x = PC1_mean, y = PC2_mean,
                        ymin = PC2_mean - PC2_se, ymax = PC2_mean + PC2_se),
                    width = 0, linewidth = 0.8, color = "black") +
      geom_errorbarh(data = pca_mean,
                     aes(x = PC1_mean, y = PC2_mean,
                         xmin = PC1_mean - PC1_se, xmax = PC1_mean + PC1_se),
                     height = 0, linewidth = 0.8, color = "black") +
      scale_fill_manual(values = visit_cols,
                        guide = guide_legend(override.aes = list(shape = 22, size = 5))) +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 1) +
      geom_vline(xintercept = 0, color = "grey50", linewidth = 1) +
      labs(title = paste0("PCA: Mean Trajectory [", mode, "]"),
           x = pc1_lab_lmm, y = pc2_lab_lmm) +
      theme_minimal(base_size = 16) +
      theme(panel.grid = element_blank(), axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"), axis.title = element_text(size = 30),
            axis.text = element_text(size = 30), plot.title = element_text(size = 30, face = "bold"),
            legend.title = element_blank(), legend.text = element_text(size = 30))
    
    print(p_pca_mean)
    ggsave(file.path(out_dir, paste0("LMM_PCA_Mean_Trajectory_", mode, ".png")),
           p_pca_mean, width = 12, height = 9, units = "in", dpi = 600)
    ggsave(file.path(out_dir, paste0("LMM_PCA_Mean_Trajectory_", mode, ".pdf")),
           p_pca_mean, width = 12, height = 9, units = "in")
    
  } else {
    cat("\nNo significant metabolites found at FDR < 0.05 for", mode, "\n")
  }
  
  cat("\n── Done:", mode, "──────────────────────────────────────────────\n")
}

cat("\n╔══════════════════════════════════════════════╗\n")
cat("  All ion modes complete!\n")
cat("╚══════════════════════════════════════════════╝\n")































####Paper Review

####LC-MS analysis:


#####
##### Raw data import, metadata rearrangement, and QC/Blank filtering
##### (LC-MS, Positive and Negative ion modes)
#####

library(ggplot2) 
library(dplyr)
library(titanic)
library(cowplot)
library(missForest)
library(matrixStats)
library(tidyverse)
library(KODAMA)
library(vegan)
library(MASS)
library(reshape2)
library(knitr)
library(ggfortify)
library(ggdendroplot)
library(pheatmap)
library(S4Vectors)
library(SummarizedExperiment)
library(pmp)
library(gridExtra)
library(mdatools)
library(ropls)
library(mixOmics)
library(ggrepel)
library(RColorBrewer)
library(scales)
library(colorRamp2)
library(circlize)
library(grid)
library(readr)
library(stringr)
library(tools)

# ============================================================
# --- Define File Paths Configuration
# ============================================================
if (.Platform$OS.type == "unix") {
  # macOS paths
  LCMS_dir      <- "/Volumes/T7/Arbeit/FIA/MSKohorte_NEU/MSKohorte LCMS"
  LC_dir_neg    <- "/Volumes/T7/Arbeit/FIA/MSKohorte_NEU/MSKohorte LCMS Neg"
  metadata_path <- "/Users/fabianschmitt/Seafile/Meine Bibliothek/FIA/EmDia_MS_Kohorte_FIMS_2025 batch and metadata/Metadata_combined_analysis.csv"
} else {
  # Windows paths
  LCMS_dir      <- "D:/Arbeit/FIA/MSKohorte_NEU/MSKohorte LCMS"
  LC_dir_neg    <- "D:/Arbeit/FIA/MSKohorte_NEU/MSKohorte LCMS Neg"
  metadata_path <- "C:/Users/fabia/Seafile/Meine Bibliothek/FIA/EmDia_MS_Kohorte_FIMS_2025 batch and metadata/Metadata_combined_analysis.csv"
}

MS_data_path_pos <- file.path(LCMS_dir, "MSKohorte_LCMS_data_short.csv")
MS_data_path_neg <- file.path(LC_dir_neg, "MSKohorte_LCMS_data_neg_short.csv")
# ============================================================
# --- Core Processing Function
# ============================================================
# This function handles the entire pipeline for a single mode.
# We will run it once for Positive and once for Negative.
process_lcms_mode <- function(MS_data_path, mode_name, metadata_path) {
  
  message("\n", strrep("=", 60))
  message("🚀 Processing ", mode_name, " (LC-MS) mode dataset...")
  message(strrep("=", 60))
  
  # --- 1. Load MS data and clean metadata rows ---
  if (!file.exists(MS_data_path)) stop("❌ File not found: ", MS_data_path)
  
  MS_data_raw <- read.csv2(MS_data_path, header = TRUE, sep = ",", dec = ".")
  
  # First 3 rows are metadata: class, index, batch
  existing_meta <- MS_data_raw[1:3, ]        # keep these
  MS_data       <- MS_data_raw[-c(1:3), ]    # remove metadata rows
  
  # --- 2. Extract sample column names ---
  sample_cols <- colnames(MS_data)[-1]       # skip SampleID column
  
  # --- 3. Extract MetabolomicsID from column names ---
  # Pattern example: XXX_Mainz_001234 or XXX_Bochum_987654
  MetabolomicsID <- str_extract(sample_cols, "(?<=_(Mainz|Bochum)_)(\\d+)$")
  
  sample_info <- data.frame(
    column_name = sample_cols,
    MetabolomicsID = MetabolomicsID,
    stringsAsFactors = FALSE
  )
  
  # --- 4. Load external metadata ---
  if (!file.exists(metadata_path)) stop("❌ External metadata file missing at: ", metadata_path)
  MS_metadata <- read.csv2(metadata_path, header = TRUE, sep = ",", dec = ".")
  
  MS_metadata$MetabolomicsID <- as.character(MS_metadata$MetabolomicsID)
  sample_info$MetabolomicsID <- as.character(sample_info$MetabolomicsID)
  
  # --- 5. Join metadata ---
  sample_info <- sample_info %>%
    left_join(MS_metadata, by = "MetabolomicsID")
  
  # Ensure correct ordering matches the original column layout
  sample_info <- sample_info[match(sample_cols, sample_info$column_name), ]
  
  # --- 6. Build metadata rows for Sex, Group, Visit, PseudoID, Age ---
  sex_row      <- data.frame(SampleID = "Sex",      t(sample_info$Sex),      stringsAsFactors = FALSE)
  group_row    <- data.frame(SampleID = "Group",    t(sample_info$Group),    stringsAsFactors = FALSE)
  visit_row    <- data.frame(SampleID = "Visit",    t(sample_info$Visit),    stringsAsFactors = FALSE)
  PseudoID_row <- data.frame(SampleID = "PseudoID", t(sample_info$PseudoID), stringsAsFactors = FALSE)
  age_row      <- data.frame(SampleID = "Age",      t(sample_info$Age),      stringsAsFactors = FALSE)
  
  # Fix column names so they match MS_data
  colnames(sex_row)      <- colnames(MS_data)
  colnames(group_row)    <- colnames(MS_data)
  colnames(visit_row)    <- colnames(MS_data)
  colnames(PseudoID_row) <- colnames(MS_data)
  colnames(age_row)      <- colnames(MS_data)
  
  # --- 7. Clean old metadata rows from feature space ---
  meta_keep <- c("class", "index", "batch")
  MS_data <- MS_data[!MS_data$SampleID %in% c(meta_keep, "Sex", "Group", "Visit", "PseudoID", "Age"), ]
  
  # --- 8. Combine everything together ---
  MS_data_new <- rbind(
    existing_meta,   # original 3 metadata rows
    sex_row,
    group_row,
    visit_row,
    PseudoID_row,
    age_row,         # added new metadata row
    MS_data          # clean profile data matrix
  )
  
  # Fill NA elements safely with empty characters
  MS_data_new[is.na(MS_data_new)] <- ""
  
  # --- 9. Detect unmatched sample IDs (warning block) ---
  unmatched <- sample_info %>% filter(!is.na(MetabolomicsID) & (is.na(Visit) | is.na(Group)))
  if (nrow(unmatched) > 0) {
    message("⚠️ Warning: Some Metabolomics IDs did not match metadata in ", mode_name, ":")
    print(unmatched$column_name)
  }
  
  # --- 10. Transpose the Output Data Frame ---
  row_headers <- MS_data_new$SampleID
  MS_data_transposed <- as.data.frame(t(MS_data_new[, -1]))
  colnames(MS_data_transposed) <- row_headers
  
  MS_data_transposed <- cbind(SampleID = rownames(MS_data_transposed), MS_data_transposed)
  rownames(MS_data_transposed) <- NULL
  
  # --- 11. Save output file ---
  MS_data_dir          <- dirname(MS_data_path)
  MS_data_filename     <- basename(MS_data_path)
  MS_data_new_filename <- sub("\\.csv$", "_with_age.csv", MS_data_filename)
  MS_data_new_path     <- file.path(MS_data_dir, MS_data_new_filename)
  
  write.csv(MS_data_transposed, file = MS_data_new_path, row.names = FALSE, quote = TRUE)
  message("✅ Transposed file saved to: ", MS_data_new_path)
  
  
  # ============================================================
  # --- Filtering block (CV < 25%, Blank > 4x)
  # ============================================================
  message("\n▶  Filtering mode: ", mode_name)
  
  meta_cols <- c("class", "index", "batch", "Sex", "Group", "Visit", "PseudoID", "Age")
  
  dat <- read.csv(
    MS_data_new_path,
    header      = TRUE,
    sep         = ",",
    dec         = ".",
    check.names = FALSE,
    colClasses  = "character"
  )
  
  rownames(dat) <- dat[, 1]
  dat           <- dat[, -1, drop = FALSE]
  colnames(dat) <- trimws(colnames(dat))
  
  feature_cols <- setdiff(colnames(dat), meta_cols)
  message("   Features before filtering : ", length(feature_cols))
  
  # --- QC CV filter (< 25 %) ---
  qc_data <- dat %>% filter(class == "QC")
  if (nrow(qc_data) == 0) stop("No 'QC' samples found in: ", MS_data_new_path)
  
  qc_cvs <- sapply(feature_cols, function(f) {
    v <- as.numeric(qc_data[[f]])
    m <- mean(v, na.rm = TRUE)
    s <- sd(v,   na.rm = TRUE)
    if (is.na(m) || m == 0) return(NA_real_)
    (s / m) * 100
  })
  
  features_pass_cv <- names(qc_cvs)[!is.na(qc_cvs) & qc_cvs < 25]
  message("   Features passing QC CV < 25 %  : ", length(features_pass_cv))
  
  # --- Blank filter (QC/Blank >= 4x) ---
  blank_data <- dat %>% filter(class == "Blank_processing")
  
  qc_means    <- sapply(feature_cols, function(f) mean(as.numeric(qc_data[[f]]),    na.rm = TRUE))
  blank_means <- sapply(feature_cols, function(f) mean(as.numeric(blank_data[[f]]), na.rm = TRUE))
  blank_means[blank_means == 0 | is.na(blank_means)] <- NA_real_
  
  ratio               <- qc_means / blank_means
  features_pass_blank <- names(ratio)[!is.na(ratio) & ratio >= 4]
  message("   Features passing QC/Blank >= 4x : ", length(features_pass_blank))
  
  # --- Combine & subset ---
  features_to_keep <- intersect(features_pass_cv, features_pass_blank)
  message("   Features retained (both filters): ", length(features_to_keep))
  
  valid_cols   <- intersect(c(meta_cols, features_to_keep), colnames(dat))
  dat_filtered <- dat[, valid_cols, drop = FALSE]
  
  # Replace NA / "NA" with empty string
  dat_filtered[is.na(dat_filtered) | dat_filtered == "NA"] <- ""
  
  message("   Final dimensions: ", nrow(dat_filtered), " samples x ", ncol(dat_filtered), " columns")
  
  # --- Save final filtered file ---
  out_file <- file.path(
    dirname(MS_data_new_path),
    paste0(tools::file_path_sans_ext(basename(MS_data_new_path)), "_CV_Blank_filtered.csv")
  )
  
  write.table(
    cbind(SampleID = rownames(dat_filtered), dat_filtered),
    file      = out_file,
    sep       = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote     = TRUE
  )
  
  message("✅ Saved final filtered matrix: ", out_file)
  message(strrep("─", 60))
  message("Done — ", mode_name, " mode processed.\n")
}

# ============================================================
# --- Execute Pipeline for Both Modes
# ============================================================

# 1. Process Positive Mode
process_lcms_mode(
  MS_data_path  = MS_data_path_pos, 
  mode_name     = "Positive", 
  metadata_path = metadata_path
)

# 2. Process Negative Mode
process_lcms_mode(
  MS_data_path  = MS_data_path_neg, 
  mode_name     = "Negative", 
  metadata_path = metadata_path
)






# ─────────────────────────────────────────────────────────────────────────────
# Load libraries
# ─────────────────────────────────────────────────────────────────────────────
library(tidyverse)
library(ggpubr)
library(patchwork)
library(lme4)
library(lmerTest)
library(pheatmap)
library(RColorBrewer)

# ─────────────────────────────────────────────────────────────────────────────
# File paths (same source data as the import/filtering script) — redefined
# here so this script can also be run on its own, not only chained after it
# in the same session.
# ─────────────────────────────────────────────────────────────────────────────
if (.Platform$OS.type == "unix") {
  # macOS paths
  LCMS_dir   <- "/Volumes/T7/Arbeit/FIA/MSKohorte_NEU/MSKohorte LCMS"
  LC_dir_neg <- "/Volumes/T7/Arbeit/FIA/MSKohorte_NEU/MSKohorte LCMS Neg"
} else {
  # Windows paths
  LCMS_dir   <- "D:/Arbeit/FIA/MSKohorte_NEU/MSKohorte LCMS"
  LC_dir_neg <- "D:/Arbeit/FIA/MSKohorte_NEU/MSKohorte LCMS Neg"
}

MS_data_path_pos <- file.path(LCMS_dir,   "MSKohorte_LCMS_data_short.csv")
MS_data_path_neg <- file.path(LC_dir_neg, "MSKohorte_LCMS_data_neg_short.csv")

# The import/filtering script derives "..._with_age_CV_Blank_filtered.csv"
# from the original file name — reproduce that here so both scripts always
# agree on the filename without hardcoding it twice.
filtered_path <- function(MS_data_path) {
  file.path(
    dirname(MS_data_path),
    sub("\\.csv$", "_with_age_CV_Blank_filtered.csv", basename(MS_data_path))
  )
}

MSKohorte_path_pos <- filtered_path(MS_data_path_pos)
MSKohorte_path_neg <- filtered_path(MS_data_path_neg)
# ─────────────────────────────────────────────────────────────────────────────
# Parameters (shared across both modes)
# ─────────────────────────────────────────────────────────────────────────────
meta_cols    <- c("class", "index", "batch", "Sex", "Group", "Visit", "PseudoID", "Age")
valid_visits <- c("T1", "T2", "T3", "T4")

# ─────────────────────────────────────────────────────────────────────────────
# Core Analysis Function
# ─────────────────────────────────────────────────────────────────────────────
# This function handles the entire PCA / CV-plot / LMM / trajectory analysis
# for a single mode. We will run it once for Positive and once for Negative.
process_lcms_analysis <- function(MSKohorte, mode) {
  
  out_dir <- dirname(MSKohorte)
  
  cat("\n\n")
  cat("╔══════════════════════════════════════════════╗\n")
  cat(sprintf("  Processing: %s ion mode (LC-MS)\n", mode))
  cat(sprintf("  File: %s\n", basename(MSKohorte)))
  cat("╚══════════════════════════════════════════════╝\n\n")
  
  if (!file.exists(MSKohorte)) {
    warning("File not found, skipping: ", MSKohorte)
    return(invisible(NULL))
  }
  
  # ── Load ──────────────────────────────────────────────────────────────────
  MSKohorte_data_transposed <- read.csv(
    MSKohorte,
    header      = TRUE,
    sep         = ",",
    dec         = ".",
    check.names = FALSE,
    colClasses  = "character"
  )
  
  rownames(MSKohorte_data_transposed) <- MSKohorte_data_transposed[, 1]
  MSKohorte_data_transposed           <- MSKohorte_data_transposed[, -1, drop = FALSE]
  MSKohorte_data_transposed[["class"]] <- as.character(MSKohorte_data_transposed[["class"]])
  
  feature_cols <- setdiff(colnames(MSKohorte_data_transposed), meta_cols)
  message("Loaded: ", nrow(MSKohorte_data_transposed), " samples × ", length(feature_cols), " features")
  
  # ── PCA of full dataset ───────────────────────────────────────────────────
  pca_features <- MSKohorte_data_transposed[-c(1:8)]
  pca_features <- apply(pca_features, 2, as.numeric)
  
  pca_res <- prcomp(pca_features, center = TRUE, scale. = TRUE)
  
  pca_df <- data.frame(
    MSKohorte_data_transposed[, intersect(meta_cols, colnames(MSKohorte_data_transposed)), drop = FALSE],
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2]
  )
  
  var_expl <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
  pc1_lab  <- paste0("PC1 (", var_expl[1], "%)")
  pc2_lab  <- paste0("PC2 (", var_expl[2], "%)")
  
  n_groups   <- length(unique(pca_df$class))
  brewer_cols <- colorRampPalette(brewer.pal(min(n_groups, 8), "Dark2"))(n_groups)
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = class)) +
    geom_point(shape = 21, size = 5, color = "black", alpha = 0.85, stroke = 0.3) +
    scale_fill_manual(values = brewer_cols,
                      guide = guide_legend(override.aes = list(shape = 22, size = 5))) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 1) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 1) +
    labs(title = paste0("PCA MSKohorte LC-MS [ESI", mode, "]"), x = pc1_lab, y = pc2_lab) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid   = element_blank(),
      axis.line    = element_line(color = "black"),
      axis.ticks   = element_line(color = "black"),
      axis.title   = element_text(size = 30),
      axis.text    = element_text(size = 30),
      plot.title   = element_text(size = 30, face = "bold"),
      legend.title = element_blank(),
      legend.text  = element_text(size = 30)
    )
  
  print(p_pca)
  ggsave(file.path(out_dir, paste0("PCA_MSKohorte_LCMS_full_data_", mode, ".png")),
         p_pca, width = 9, height = 7, dpi = 400)
  
  # ── CV plot ───────────────────────────────────────────────────────────────
  MSKohorte_data_transposed_CV <- MSKohorte_data_transposed %>%
    filter(xor(class == "QC", class == "Sample")) %>%
    dplyr::select(-c(index, Sex, Group, Visit, PseudoID, Age)) %>%
    rownames_to_column(var = "SampleID")
  
  data_long_cv <- pivot_longer(MSKohorte_data_transposed_CV,
                               cols = -c(SampleID, batch, class),
                               names_to = "Metabolite", values_to = "Value") %>%
    mutate(Value = as.numeric(Value))
  
  cv_data <- data_long_cv %>%
    group_by(batch, class, Metabolite) %>%
    summarise(
      Mean_Value = mean(Value, na.rm = TRUE),
      SD_Value   = sd(Value,   na.rm = TRUE),
      CV         = (SD_Value / Mean_Value) * 100,
      .groups    = "drop"
    ) %>%
    filter(Mean_Value != 0)
  
  p_cv <- ggplot(cv_data, aes(x = class, y = CV, fill = batch)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 position = position_dodge(width = 0.9), color = "black") +
    geom_hline(yintercept = 25, linetype = "dotted", color = "red", linewidth = 1) +
    ggplot2::annotate("text", x = 1.5, y = 10, label = "25% threshold",
                      color = "red", size = 5, hjust = 0) +
    coord_cartesian(ylim = c(0, 200)) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title    = paste0("CV of Metabolites by Batch and Class (LC-MS, ", mode, ")"),
      subtitle = "Violin and boxplot (CV capped at 200%)",
      x = "Class", y = "Coefficient of Variation (%)", fill = "Batch"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x  = element_text(size = 16, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 16),
      axis.title   = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
      plot.title   = element_text(size = 20, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 16, hjust = 0.5)
    )
  
  print(p_cv)
  ggsave(file.path(out_dir, paste0("CVs_MSKohorte_LCMS_", mode, ".png")),
         p_cv, width = 9, height = 7, dpi = 400)
  
  # ── Prepare long-format data for LMM ─────────────────────────────────────
  MS <- as.data.frame(MSKohorte_data_transposed, stringsAsFactors = FALSE) %>%
    rownames_to_column("Sample") %>%
    as_tibble() %>%
    filter(Visit %in% valid_visits & !is.na(Visit))
  
  lmm_meta_cols <- c("class", "index", "batch", "Sex", "Group", "Visit", "PseudoID")
  metadata_available <- intersect(c("Sample", lmm_meta_cols), colnames(MS))
  metadata        <- MS[, metadata_available]
  metabolite_cols <- setdiff(colnames(MS), colnames(metadata))
  
  data_long <- MS[, metabolite_cols] %>%
    mutate(Sample = MS$Sample) %>%
    pivot_longer(cols = -Sample, names_to = "Metabolite", values_to = "Intensity") %>%
    left_join(metadata, by = "Sample") %>%
    dplyr::filter(class == "Sample") %>%
    dplyr::select(PseudoID, Metabolite, Visit, Intensity) %>%
    mutate(
      # Crucial Fix: Force the main Intensity column to be numeric
      # so your original downstream plotting loops do not crash
      Intensity     = as.numeric(Intensity),
      # Create the log-transformed column strictly for the LMM analysis
      Intensity_Log = log2(Intensity + 1), 
      Visit         = factor(Visit, levels = valid_visits, ordered = TRUE)
    ) %>%
    filter(!is.na(PseudoID) & PseudoID != "" & !is.na(Intensity))
  
  # ── Linear Mixed Models (Run on transformed data) ────────────────────────
  cat("Running Linear Mixed Models on log2-transformed data...\n")
  
  lmm_results <- data_long %>%
    group_by(Metabolite) %>%
    summarise(
      n_obs      = n(),
      n_patients = n_distinct(PseudoID),
      n_visits   = n_distinct(Visit),
      p_value    = {
        if (n_obs >= 20 && n_patients >= 5 && n_visits >= 2) {
          tryCatch({
            # Model fits on the stabilized log2 scale
            fit <- lmer(Intensity_Log ~ Visit + (1 | PseudoID), data = pick(everything()))
            anova(fit)["Visit", "Pr(>F)"]
          }, error = function(e) NA_real_)
        } else NA_real_
      },
      .groups = "drop"
    ) %>%
    filter(!is.na(p_value)) %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
  
  cat("Total metabolites tested:", nrow(lmm_results), "\n")
  cat("Significant at FDR < 0.05:", sum(lmm_results$FDR < 0.05, na.rm = TRUE), "\n")
  cat("Significant at FDR < 0.10:", sum(lmm_results$FDR < 0.10, na.rm = TRUE), "\n")
  print(head(lmm_results, 20))
  
  write_csv(lmm_results, file.path(out_dir, paste0("LMM_Visit_Results_LCMS_", mode, ".csv")))
  
  sig_mets <- lmm_results %>% filter(FDR < 0.05) %>% pull(Metabolite)
  
  # ── Plots ─────────────────────────────────────────────────────────────────
  
  # ── Helper: fixed-format scientific notation for axis labels ────────────────
  format_scientific <- function(x) {
    ifelse(x == 0, "0", sprintf("%.1E", x))
  }
  
  if (length(sig_mets) > 0) {
    
    sig_mets_top15 <- lmm_results %>%
      filter(FDR < 0.05) %>% arrange(p_value) %>% slice_head(n = 15) %>% pull(Metabolite)
    
    cat("\nGenerating plots for top", length(sig_mets_top15), "significant metabolites...\n")
    
    summary_stats <- data_long %>%
      filter(Metabolite %in% sig_mets) %>%
      group_by(Metabolite, Visit) %>%
      summarise(
        mean_int = mean(Intensity, na.rm = TRUE),
        se       = sd(Intensity, na.rm = TRUE) / sqrt(n()),
        .groups  = "drop"
      )
    
    n_sig       <- length(sig_mets_top15)
    ncol_plot   <- min(3, n_sig)
    height_plot <- max(3, min(9, ceiling(n_sig / ncol_plot) * 1.5))
    
    # Boxplots
    p_box <- ggplot(data_long %>% filter(Metabolite %in% sig_mets_top15),
                    aes(x = Visit, y = Intensity)) +
      geom_boxplot(aes(fill = Visit), outlier.size = 1, outlier.alpha = 0.5, linewidth = 0.5) +
      facet_wrap(~Metabolite, scales = "free_y", ncol = ncol_plot) +
      scale_y_continuous(n.breaks = 4, labels = format_scientific) +
      labs(y = "Relative intensity", x = "Visit") +
      theme_bw(base_size = 11) +
      theme(
        strip.text       = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "grey60"),
        panel.grid.major = element_line(linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position  = "none",
        axis.text        = element_text(size = 8),
        axis.title       = element_text(size = 8, face = "bold"),
        plot.margin      = margin(8, 8, 8, 8)
      )
    
    print(p_box)
    ggsave(file.path(out_dir, paste0("LMM_Significant_Metabolites_Boxplots_LCMS_", mode, ".png")),
           p_box, width = 7, height = height_plot, units = "in", dpi = 600)
    ggsave(file.path(out_dir, paste0("LMM_Significant_Metabolites_Boxplots_LCMS_", mode, ".pdf")),
           p_box, width = 7, height = height_plot, units = "in")
    
    # Combined spaghetti (top 15)
    summary_stats_top15 <- summary_stats %>% filter(Metabolite %in% sig_mets_top15)
    
    p_spaghetti <- ggplot(data_long %>% filter(Metabolite %in% sig_mets_top15),
                          aes(x = Visit, y = Intensity)) +
      geom_line(aes(group = PseudoID), alpha = 0.15, linewidth = 0.35, color = "grey60") +
      geom_ribbon(data = summary_stats_top15,
                  aes(x = Visit, y = mean_int,
                      ymin = mean_int - se, ymax = mean_int + se, group = 1),
                  fill = "#1B9E77", alpha = 0.25, inherit.aes = FALSE) +
      geom_line(data = summary_stats_top15,
                aes(x = Visit, y = mean_int, group = 1),
                color = "#1B9E77", linewidth = 1, inherit.aes = FALSE) +
      geom_point(data = summary_stats_top15,
                 aes(x = Visit, y = mean_int, group = 1),
                 color = "#1B9E77", size = 2.5, shape = 21, fill = "white",
                 stroke = 1.2, inherit.aes = FALSE) +
      facet_wrap(~Metabolite, scales = "free_y", ncol = ncol_plot) +
      scale_y_continuous(n.breaks = 4, labels = format_scientific) +
      labs(y = "Relative intensity", x = "Visit") +
      theme_bw(base_size = 11) +
      theme(
        strip.text       = element_text(size = 8, face = "bold", color = "grey20"),
        strip.background = element_rect(fill = "grey92", color = "grey60", linewidth = 0.5),
        panel.grid.major = element_line(linewidth = 0.25, color = "grey85"),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(color = "grey60", linewidth = 0.5),
        legend.position  = "none",
        axis.text        = element_text(size = 8, color = "grey20"),
        axis.title       = element_text(size = 8, face = "bold"),
        axis.ticks       = element_line(linewidth = 0.4),
        plot.margin      = margin(8, 8, 8, 8)
      )
    
    print(p_spaghetti)
    ggsave(file.path(out_dir, paste0("LMM_Significant_Metabolites_Spaghetti_LCMS_", mode, ".png")),
           p_spaghetti, width = 7, height = height_plot, units = "in", dpi = 600)
    ggsave(file.path(out_dir, paste0("LMM_Significant_Metabolites_Spaghetti_LCMS_", mode, ".pdf")),
           p_spaghetti, width = 7, height = height_plot, units = "in")
    
    # Individual spaghetti plots (all significant metabolites)
    spaghetti_dir <- file.path(out_dir, paste0("spaghetti_LCMS_", mode))
    if (!dir.exists(spaghetti_dir)) dir.create(spaghetti_dir)
    
    wrap_text <- function(text, width = 30) {
      strwrap(text, width = width, simplify = FALSE)[[1]] %>% paste(collapse = "\n")
    }
    
    cat("\nExporting", length(sig_mets), "individual spaghetti plots...\n")
    
    for (met in sig_mets) {
      met_data    <- data_long     %>% filter(Metabolite == met)
      met_summary <- summary_stats %>% filter(Metabolite == met)
      
      p_individual <- ggplot(met_data, aes(x = Visit, y = Intensity)) +
        geom_line(aes(group = PseudoID), alpha = 0.2, linewidth = 0.4, color = "grey60") +
        geom_ribbon(data = met_summary,
                    aes(x = Visit, y = mean_int,
                        ymin = mean_int - se, ymax = mean_int + se, group = 1),
                    fill = "#1B9E77", alpha = 0.3, inherit.aes = FALSE) +
        geom_line(data = met_summary,
                  aes(x = Visit, y = mean_int, group = 1),
                  color = "#1B9E77", linewidth = 1.2, inherit.aes = FALSE) +
        geom_point(data = met_summary,
                   aes(x = Visit, y = mean_int, group = 1),
                   color = "#1B9E77", size = 3, shape = 21, fill = "white",
                   stroke = 1.3, inherit.aes = FALSE) +
        scale_y_continuous(n.breaks = 4, labels = format_scientific) +
        labs(title = wrap_text(met), y = "Relative intensity", x = "Visit") +
        theme_bw(base_size = 9) +
        theme(
          plot.title       = element_text(size = 8, face = "bold", hjust = 0.5, lineheight = 0.75),
          panel.grid.major = element_line(linewidth = 0.25, color = "grey85"),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(color = "grey60", linewidth = 0.6),
          axis.text        = element_text(size = 8, color = "grey20"),
          axis.title       = element_text(size = 9, face = "bold"),
          axis.ticks       = element_line(linewidth = 0.4),
          plot.margin      = margin(15, 3, 3, 3)
        )
      
      safe_name <- gsub("[/\\\\:*?\"<>|]", "_", met)
      ggsave(file.path(spaghetti_dir, paste0(safe_name, "_spaghetti.png")),
             p_individual, width = 2.3, height = 2, units = "in", dpi = 600)
    }
    
    cat("Exported", length(sig_mets), "individual plots to:", spaghetti_dir, "\n")
    
    # ── Individual spaghetti plots, MEDIAN version (extra output) ───────────
    summary_stats_median <- data_long %>%
      filter(Metabolite %in% sig_mets) %>%
      group_by(Metabolite, Visit) %>%
      summarise(
        median_int = median(Intensity, na.rm = TRUE),
        q25        = quantile(Intensity, 0.25, na.rm = TRUE),
        q75        = quantile(Intensity, 0.75, na.rm = TRUE),
        .groups    = "drop"
      )
    
    spaghetti_median_dir <- file.path(out_dir, paste0("spaghetti_LCMS_", mode, "_median"))
    if (!dir.exists(spaghetti_median_dir)) dir.create(spaghetti_median_dir)
    
    cat("\nExporting", length(sig_mets), "individual spaghetti plots (median)...\n")
    
    for (met in sig_mets) {
      met_data    <- data_long              %>% filter(Metabolite == met)
      met_summary <- summary_stats_median   %>% filter(Metabolite == met)
      
      p_individual_median <- ggplot(met_data, aes(x = Visit, y = Intensity)) +
        geom_line(aes(group = PseudoID), alpha = 0.2, linewidth = 0.4, color = "grey60") +
        geom_ribbon(data = met_summary,
                    aes(x = Visit, y = median_int,
                        ymin = q25, ymax = q75, group = 1),
                    fill = "#1B9E77", alpha = 0.3, inherit.aes = FALSE) +
        geom_line(data = met_summary,
                  aes(x = Visit, y = median_int, group = 1),
                  color = "#1B9E77", linewidth = 1.2, inherit.aes = FALSE) +
        geom_point(data = met_summary,
                   aes(x = Visit, y = median_int, group = 1),
                   color = "#1B9E77", size = 3, shape = 21, fill = "white",
                   stroke = 1.3, inherit.aes = FALSE) +
        scale_y_continuous(n.breaks = 4, labels = format_scientific) +
        labs(title = wrap_text(met), y = "Relative intensity", x = "Visit") +
        theme_bw(base_size = 9) +
        theme(
          plot.title       = element_text(size = 8, face = "bold", hjust = 0.5, lineheight = 0.75),
          panel.grid.major = element_line(linewidth = 0.25, color = "grey85"),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(color = "grey60", linewidth = 0.6),
          axis.text        = element_text(size = 8, color = "grey20"),
          axis.title       = element_text(size = 9, face = "bold"),
          axis.ticks       = element_line(linewidth = 0.4),
          plot.margin      = margin(15, 3, 3, 3)
        )
      
      safe_name <- gsub("[/\\\\:*?\"<>|]", "_", met)
      ggsave(file.path(spaghetti_median_dir, paste0(safe_name, "_spaghetti_median.png")),
             p_individual_median, width = 2.3, height = 2, units = "in", dpi = 600)
    }
    
    cat("Exported", length(sig_mets), "individual median plots to:", spaghetti_median_dir, "\n")
    
    # ── Trajectory clustering heatmap ────────────────────────────────────────
    trajectory_data <- data_long %>%
      filter(Metabolite %in% sig_mets_top15) %>%
      group_by(Metabolite, Visit) %>%
      summarise(mean_intensity = mean(Intensity, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Visit, values_from = mean_intensity) %>%
      column_to_rownames("Metabolite")
    
    trajectory_scaled <- t(scale(t(trajectory_data)))
    hc <- hclust(dist(trajectory_scaled), method = "ward.D2")
    
    p_heatmap <- pheatmap(
      trajectory_scaled,
      cluster_rows     = TRUE,
      cluster_cols     = FALSE,
      clustering_method = "ward.D2",
      color            = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      breaks           = seq(-2, 2, length.out = 101),
      main             = paste0("Metabolite Trajectory Patterns [LC-MS ", mode, "]"),
      fontsize         = 8, fontsize_row = 8, fontsize_col = 8,
      angle_col        = 0, border_color = "grey80",
      cellwidth = 15, cellheight = 15,
      legend_breaks  = c(-2, -1, 0, 1, 2),
      legend_labels  = c("-2 SD", "-1 SD", "Mean", "+1 SD", "+2 SD"),
      silent = TRUE
    )
    
    ggsave(file.path(out_dir, paste0("LMM_Trajectory_Heatmap_LCMS_", mode, ".png")),
           p_heatmap, width = 5, height = max(6, nrow(trajectory_data) * 0.3),
           units = "in", dpi = 600)
    ggsave(file.path(out_dir, paste0("LMM_Trajectory_Heatmap_LCMS_", mode, ".pdf")),
           p_heatmap, width = 5, height = max(6, nrow(trajectory_data) * 0.3),
           units = "in")
    
    # Trajectory cluster line plot
    n_clusters  <- min(3, max(2, ceiling(n_sig / 5)))
    clusters    <- cutree(hc, k = n_clusters)
    
    trajectory_clustered <- trajectory_data %>%
      rownames_to_column("Metabolite") %>%
      mutate(Cluster = paste0("Pattern ", clusters[Metabolite])) %>%
      pivot_longer(cols = any_of(valid_visits), names_to = "Visit", values_to = "Mean_Intensity") %>%
      mutate(Visit = factor(Visit, levels = valid_visits, ordered = TRUE))
    
    p_clusters <- ggplot(trajectory_clustered,
                         aes(x = Visit, y = Mean_Intensity, group = Metabolite, color = Cluster)) +
      geom_line(linewidth = 1, alpha = 0.8) +
      geom_point(size = 2, alpha = 0.8) +
      facet_wrap(~Cluster, scales = "free_y", ncol = n_clusters) +
      scale_color_brewer(palette = "Set1") +
      labs(y = "Mean relative intensity", x = "Visit",
           title = paste0("Distinct Temporal Patterns [LC-MS ", mode, "]")) +
      theme_bw(base_size = 11) +
      theme(
        strip.text       = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "grey60"),
        panel.grid.major = element_line(linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position  = "none",
        axis.text        = element_text(size = 9),
        axis.title       = element_text(size = 10, face = "bold"),
        plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.margin      = margin(8, 8, 8, 8)
      )
    
    print(p_clusters)
    ggsave(file.path(out_dir, paste0("LMM_Trajectory_Clusters_LCMS_", mode, ".png")),
           p_clusters, width = 7, height = 3, units = "in", dpi = 600)
    ggsave(file.path(out_dir, paste0("LMM_Trajectory_Clusters_LCMS_", mode, ".pdf")),
           p_clusters, width = 7, height = 3, units = "in")
    
    cluster_summary <- trajectory_clustered %>%
      dplyr::select(Metabolite, Cluster) %>% distinct() %>% arrange(Cluster, Metabolite)
    write_csv(cluster_summary, file.path(out_dir, paste0("LMM_Trajectory_Clusters_LCMS_", mode, ".csv")))
    
    # ── PCA with patient trajectories ────────────────────────────────────────
    MS_pca <- as.data.frame(MSKohorte_data_transposed, stringsAsFactors = FALSE) %>%
      rownames_to_column("Sample") %>% as_tibble() %>%
      filter(Visit %in% valid_visits & !is.na(Visit),
             class == "Sample",
             !is.na(PseudoID) & PseudoID != "")
    
    pca_meta_cols  <- c("Sample", "class", "index", "batch", "Sex", "Group", "Visit", "PseudoID", "Age")
    pca_met_cols   <- setdiff(colnames(MS_pca), pca_meta_cols)
    pca_matrix     <- MS_pca %>% dplyr::select(all_of(pca_met_cols)) %>%
      mutate(across(everything(), as.numeric)) %>% as.matrix()
    pca_meta_sub   <- MS_pca %>% dplyr::select(PseudoID, Visit)
    
    complete_rows  <- complete.cases(pca_matrix)
    pca_matrix     <- pca_matrix[complete_rows, ]
    pca_meta_sub   <- pca_meta_sub[complete_rows, ]
    valid_cols     <- apply(pca_matrix, 2, function(x) !is.na(sd(x, na.rm=TRUE)) && sd(x, na.rm=TRUE) > 0)
    pca_matrix     <- pca_matrix[, valid_cols]
    
    pca_result     <- prcomp(pca_matrix, scale. = TRUE, center = TRUE)
    var_explained  <- summary(pca_result)$importance[2, 1:2] * 100
    pc1_lab_lmm    <- sprintf("PC1 (%.1f%%)", var_explained[1])
    pc2_lab_lmm    <- sprintf("PC2 (%.1f%%)", var_explained[2])
    
    pca_traj_df <- as.data.frame(pca_result$x[, 1:2]) %>%
      mutate(PseudoID = pca_meta_sub$PseudoID,
             Visit    = factor(pca_meta_sub$Visit, levels = valid_visits))
    
    n_v          <- length(valid_visits)
    visit_cols   <- colorRampPalette(brewer.pal(min(n_v, 8), "Dark2"))(n_v)
    
    pca_mean <- pca_traj_df %>%
      group_by(Visit) %>%
      summarise(PC1_mean = mean(PC1), PC2_mean = mean(PC2),
                PC1_se   = sd(PC1) / sqrt(n()),
                PC2_se   = sd(PC2) / sqrt(n()), .groups = "drop")
    
    p_pca_traj <- ggplot(pca_traj_df, aes(x = PC1, y = PC2, fill = Visit)) +
      geom_path(aes(group = PseudoID), color = "grey60", linewidth = 0.5, alpha = 0.4,
                arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
      geom_point(shape = 21, size = 5, color = "black", alpha = 0.85, stroke = 0.3) +
      scale_fill_manual(values = visit_cols,
                        guide = guide_legend(override.aes = list(shape = 22, size = 5))) +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 1) +
      geom_vline(xintercept = 0, color = "grey50", linewidth = 1) +
      labs(title = paste0("PCA: Patient Trajectories [LC-MS ", mode, "]"),
           x = pc1_lab_lmm, y = pc2_lab_lmm) +
      theme_minimal(base_size = 16) +
      theme(panel.grid = element_blank(), axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"), axis.title = element_text(size = 30),
            axis.text = element_text(size = 30), plot.title = element_text(size = 30, face = "bold"),
            legend.title = element_blank(), legend.text = element_text(size = 30))
    
    print(p_pca_traj)
    ggsave(file.path(out_dir, paste0("LMM_PCA_Trajectories_LCMS_", mode, ".png")),
           p_pca_traj, width = 12, height = 9, units = "in", dpi = 600)
    
    p_pca_mean <- ggplot(pca_traj_df, aes(x = PC1, y = PC2, fill = Visit)) +
      geom_point(shape = 21, size = 3, color = "black", alpha = 0.3, stroke = 0.2) +
      geom_path(data = pca_mean, aes(x = PC1_mean, y = PC2_mean, group = 1),
                color = "black", linewidth = 1.5,
                arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
      geom_point(data = pca_mean, aes(x = PC1_mean, y = PC2_mean, fill = Visit),
                 shape = 21, size = 8, color = "black", stroke = 0.5) +
      geom_errorbar(data = pca_mean,
                    aes(x = PC1_mean, y = PC2_mean,
                        ymin = PC2_mean - PC2_se, ymax = PC2_mean + PC2_se),
                    width = 0, linewidth = 0.8, color = "black") +
      geom_errorbarh(data = pca_mean,
                     aes(x = PC1_mean, y = PC2_mean,
                         xmin = PC1_mean - PC1_se, xmax = PC1_mean + PC1_se),
                     height = 0, linewidth = 0.8, color = "black") +
      scale_fill_manual(values = visit_cols,
                        guide = guide_legend(override.aes = list(shape = 22, size = 5))) +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 1) +
      geom_vline(xintercept = 0, color = "grey50", linewidth = 1) +
      labs(title = paste0("PCA: Mean Trajectory [LC-MS ", mode, "]"),
           x = pc1_lab_lmm, y = pc2_lab_lmm) +
      theme_minimal(base_size = 16) +
      theme(panel.grid = element_blank(), axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"), axis.title = element_text(size = 30),
            axis.text = element_text(size = 30), plot.title = element_text(size = 30, face = "bold"),
            legend.title = element_blank(), legend.text = element_text(size = 30))
    
    print(p_pca_mean)
    ggsave(file.path(out_dir, paste0("LMM_PCA_Mean_Trajectory_LCMS_", mode, ".png")),
           p_pca_mean, width = 12, height = 9, units = "in", dpi = 600)
    ggsave(file.path(out_dir, paste0("LMM_PCA_Mean_Trajectory_LCMS_", mode, ".pdf")),
           p_pca_mean, width = 12, height = 9, units = "in")
    
  } else {
    cat("\nNo significant metabolites found at FDR < 0.05 for", mode, "\n")
  }
  
  cat("\n── Done:", mode, "──────────────────────────────────────────────\n")
  invisible(NULL)
}

# ─────────────────────────────────────────────────────────────────────────────
# Execute analysis for both modes
# ─────────────────────────────────────────────────────────────────────────────

# 1. Process Positive Mode
process_lcms_analysis(MSKohorte = MSKohorte_path_pos, mode = "POS")

# 2. Process Negative Mode
process_lcms_analysis(MSKohorte = MSKohorte_path_neg, mode = "NEG")

cat("\n╔══════════════════════════════════════════════╗\n")
cat("  LC-MS positive and negative mode analysis complete!\n")
cat("╚══════════════════════════════════════════════╝\n")








## ===========================================================================
## MSKohorte: Load FIMS + LC-MS data (POS & NEG), match samples by shared
## MetabolomicsID, harmonize feature names, and perform method comparison
## ===========================================================================
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(broom)
library(mcr)             # For Deming regression
library(DescTools)       # For Lin's CCC (DescTools::CCC)
library(ggVennDiagram)   # For the visually pleasing feature-overlap Venn diagrams

# Defensive reassignment - avoid masking by other packages
select <- dplyr::select
filter <- dplyr::filter
count  <- dplyr::count
sym    <- rlang::sym
# ── OS-specific base paths ──────────────────────────────────────────────────
is_unix <- .Platform$OS.type == "unix"

if (is_unix) {
  # macOS path (with leading slash for mounted volume)
  base_ms_path <- "/Volumes/T7/Arbeit/FIA/MSKohorte_NEU/"
} else {
  # Windows path
  base_ms_path <- "D:/Arbeit/FIA/MSKohorte_NEU/"
}

# ── File Paths (Positive and Negative Modes) ────────────────────────────────
fims_pos_path <- file.path(base_ms_path, "10ppm annotation/CSV/feature_matrix_transposed_MSKohorte_with_metadata_with_age_CV_Blank_filtered.csv")
fims_neg_path <- file.path(base_ms_path, "10ppm annotation/CSV_NEG/feature_matrix_transposed_MSKohorte_with_metadata_neg_with_age_CV_Blank_filtered.csv")

lcms_pos_path <- file.path(base_ms_path, "MSKohorte LCMS/MSKohorte_LCMS_data_short_with_age_CV_Blank_filtered.csv")
lcms_neg_path <- file.path(base_ms_path, "MSKohorte LCMS Neg/MSKohorte_LCMS_data_neg_short_with_age_CV_Blank_filtered.csv")

# Harmonization file location
lcms_folder       <- file.path(base_ms_path, "MSKohorte LCMS")
harmonization_csv <- file.path(lcms_folder, "harmonized_metabolite_matches.csv")

meta_cols <- c("class", "index", "batch", "Sex", "Group", "Visit", "PseudoID", "Age")
# ── Directory setup ─────────────────────────────────────────────────────────
output_dir <- file.path(base_ms_path, "cross_platform_correlation_results_MSKohorte")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
message("📁 Output directory: ", output_dir)

ba_individual_dir <- file.path(output_dir, "individual_bland_altman")
scatter_dir       <- file.path(output_dir, "individual_scatter")
venn_dir          <- file.path(output_dir, "venn_diagrams")
dir.create(ba_individual_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(scatter_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(venn_dir, showWarnings = FALSE, recursive = TRUE)

# ── Helper: extract numeric MetabolomicsID ──────────────────────────────────
extract_match_id <- function(sample_id) {
  str_extract(sample_id, "(?<=_(Mainz|Bochum)_)(\\d+)$")
}

# ── ROBUST FILE READER (Fixes multibyte string <a0> errors) ──────────────────
read_csv_robust <- function(file_path, is_matrix = FALSE) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  
  df <- tryCatch({
    read.csv(file_path, header = TRUE, sep = ",", dec = ".",
             stringsAsFactors = FALSE, check.names = FALSE, 
             colClasses = "character", na.strings = c("", "NA"))
  }, error = function(e) {
    read.csv(file_path, header = TRUE, sep = ",", dec = ".",
             stringsAsFactors = FALSE, check.names = FALSE, 
             colClasses = "character", na.strings = c("", "NA"),
             fileEncoding = "latin1")
  })
  
  # Sanitize text: replace non-breaking spaces safely across encodings
  df[] <- lapply(df, function(col) {
    if (is.character(col)) {
      col <- gsub("\u00A0", " ", col, fixed = TRUE)
      col <- gsub("\xa0", " ", col, fixed = TRUE, useBytes = TRUE)
      col <- str_trim(col)
      col[col == ""] <- NA
    }
    col
  })
  
  if (is_matrix) {
    rownames(df) <- df[, 1]
    return(df[, -1, drop = FALSE])
  } else {
    return(df)
  }
}

# ── Sample Matching Function ────────────────────────────────────────────────
match_lcms_to_fims <- function(lcms_data, fims_data, mode_name = "pos") {
  lcms_ids <- tibble(
    tubecode_V2   = rownames(lcms_data),
    match_id      = extract_match_id(rownames(lcms_data)),
    PseudoID      = lcms_data$PseudoID,
    Visit         = lcms_data$Visit,
    Group         = lcms_data$Group,
    Sex           = lcms_data$Sex,
    Age           = lcms_data$Age
  )
  
  fims_ids <- tibble(
    tubeCode = rownames(fims_data),
    match_id = extract_match_id(rownames(fims_data))
  )
  
  dup_check <- lcms_ids %>%
    filter(!is.na(match_id)) %>%
    count(match_id) %>%
    filter(n > 1)
  
  if (nrow(dup_check) > 0) {
    warning("  [", toupper(mode_name), "] ", nrow(dup_check),
            " match_id value(s) appear more than once in LC-MS data.")
  }
  
  lookup <- lcms_ids %>%
    left_join(fims_ids, by = "match_id", relationship = "many-to-many") %>%
    mutate(mode = mode_name)
  
  n_before <- nrow(lookup)
  lookup <- lookup %>% distinct(tubecode_V2, .keep_all = TRUE)
  n_after <- nrow(lookup)
  
  if (n_before != n_after) {
    warning("  [", toupper(mode_name), "] dropped ", n_before - n_after,
            " duplicate row(s) after join.")
  }
  
  unmatched <- sum(is.na(lookup$tubeCode))
  if (unmatched > 0) {
    message("  [", toupper(mode_name), "] ", unmatched, " of ", nrow(lookup),
            " LC-MS samples could not be matched to a FIMS sample.")
  } else {
    message("  [", toupper(mode_name), "] all ", nrow(lookup), " LC-MS samples matched.")
  }
  
  lookup
}

# ── Helper: Safe Log2 Transformation ───────────────────────────────────────
safe_log2 <- function(x) {
  min_val <- suppressWarnings(min(x, na.rm = TRUE))
  shift <- if (is.finite(min_val) && min_val <= 0) abs(min_val) + 1 else 1
  log2(x + shift)
}

# ── Load Harmonization File ────────────────────────────────────────────────
if (!file.exists(harmonization_csv)) stop("Harmonization file not found: ", harmonization_csv)
df_harm <- read_csv_robust(harmonization_csv, is_matrix = FALSE)

# ── Single Mode Processing Function ─────────────────────────────────────────
# Returns a list(matched_data, venn) — matched_data is the paired-sample
# intensity table (or NULL if no samples/features matched), venn is the
# feature-space identity vectors used for the overlap Venn diagram/stats.
# venn is computed independently of sample matching, since it's a question
# about which metabolites were identified on each platform, not which
# patients paired up.
process_single_mode <- function(mode_name, fims_path, lcms_path, df_harm) {
  if (!file.exists(fims_path)) {
    warning("FIMS file for mode '", mode_name, "' not found at: ", fims_path)
    return(NULL)
  }
  if (!file.exists(lcms_path)) {
    warning("LC-MS file for mode '", mode_name, "' not found at: ", lcms_path)
    return(NULL)
  }
  
  message("\n--------------------------------------------------")
  message("⚡ Processing Ion Mode: ", toupper(mode_name))
  message("--------------------------------------------------")
  
  fims_data <- read_csv_robust(fims_path, is_matrix = TRUE)
  lcms_data <- read_csv_robust(lcms_path, is_matrix = TRUE)
  
  # ── Feature-level overlap (independent of sample matching) ────────────────
  fims_feature_cols <- setdiff(colnames(fims_data), meta_cols)
  lcms_feature_cols <- setdiff(colnames(lcms_data), meta_cols)
  
  harm_pairs_mode <- df_harm %>%
    filter(!is.na(fims_feature_matched) & fims_feature_matched != "" &
             !is.na(lcms_feature_full) & lcms_feature_full != "") %>%
    mutate(
      mode                   = trimws(tolower(ion_mode)),
      fims_feature           = trimws(fims_feature_matched),
      lcms_feature            = trimws(lcms_feature_full),
      harmonized_annotation  = trimws(harmonized_annotation)
    ) %>%
    filter(mode == tolower(mode_name)) %>%
    select(mode, fims_feature, lcms_feature, harmonized_annotation, molecular_formula, monoisotopic_mass)
  
  valid_fims_feats <- intersect(harm_pairs_mode$fims_feature, fims_feature_cols)
  valid_lcms_feats <- intersect(harm_pairs_mode$lcms_feature, lcms_feature_cols)
  
  features_existing <- harm_pairs_mode %>%
    filter(fims_feature %in% valid_fims_feats, lcms_feature %in% valid_lcms_feats)
  
  # Give every retained feature a shared "identity": the harmonized
  # annotation if it's matched to the other platform, otherwise a
  # platform-private label that can never collide across platforms — this
  # is what lets the Venn diagram treat "matched" as a genuine set
  # intersection rather than a raw name comparison (FIMS and LC-MS feature
  # names use completely different conventions, e.g. "123.4567" vs
  # "Glutamine", so they'd never overlap on their own).
  fims_annot_lookup <- features_existing %>% distinct(fims_feature, harmonized_annotation) %>% deframe()
  lcms_annot_lookup <- features_existing %>% distinct(lcms_feature, harmonized_annotation) %>% deframe()
  
  fims_identity <- unname(ifelse(fims_feature_cols %in% names(fims_annot_lookup),
                                 fims_annot_lookup[fims_feature_cols],
                                 paste0("FIMS_only::", fims_feature_cols)))
  lcms_identity <- unname(ifelse(lcms_feature_cols %in% names(lcms_annot_lookup),
                                 lcms_annot_lookup[lcms_feature_cols],
                                 paste0("LCMS_only::", lcms_feature_cols)))
  
  venn_result <- list(mode = mode_name, fims_set = fims_identity, lcms_set = lcms_identity)
  
  # ── Sample matching / intensity extraction ─────────────────────────────────
  match_table <- match_lcms_to_fims(lcms_data, fims_data, mode_name)
  
  n_before_drop <- nrow(match_table)
  match_table <- match_table %>%
    filter(!is.na(tubeCode), !is.na(PseudoID), PseudoID != "")
  n_after_drop <- nrow(match_table)
  
  message("  [", toupper(mode_name), "] Retained ", n_after_drop, " matched sample pairs.")
  
  if (nrow(match_table) == 0) {
    return(list(matched_data = NULL, venn = venn_result))
  }
  if (nrow(features_existing) == 0) {
    warning("No matching feature pairs found for mode: ", mode_name)
    return(list(matched_data = NULL, venn = venn_result))
  }
  
  fims_raw <- fims_data[match_table$tubeCode, features_existing$fims_feature, drop = FALSE]
  fims_long <- as.data.frame(fims_raw) %>%
    mutate(tubeCode = rownames(fims_raw)) %>%
    pivot_longer(cols = -tubeCode, names_to = "fims_feature", values_to = "fims_intensity") %>%
    mutate(fims_intensity = suppressWarnings(as.numeric(fims_intensity)))
  
  lcms_raw <- lcms_data[match_table$tubecode_V2, features_existing$lcms_feature, drop = FALSE]
  lcms_long <- as.data.frame(lcms_raw) %>%
    mutate(tubecode_V2 = rownames(lcms_raw)) %>%
    pivot_longer(cols = -tubecode_V2, names_to = "lcms_feature", values_to = "lcms_intensity") %>%
    mutate(lcms_intensity = suppressWarnings(as.numeric(lcms_intensity)))
  
  fims_long <- fims_long %>% mutate(fims_intensity_log2 = safe_log2(fims_intensity))
  lcms_long <- lcms_long %>% mutate(lcms_intensity_log2 = safe_log2(lcms_intensity))
  
  matched_data <- match_table %>%
    select(PseudoID, Visit, Group, Sex, Age, tubeCode, tubecode_V2) %>%
    left_join(fims_long, by = "tubeCode", relationship = "many-to-many") %>%
    left_join(lcms_long, by = "tubecode_V2", relationship = "many-to-many") %>%
    inner_join(features_existing, by = c("fims_feature", "lcms_feature")) %>%
    mutate(ion_mode = mode_name) %>%
    select(
      PseudoID, Visit, Group, Sex, Age, ion_mode,
      harmonized_annotation, fims_feature, lcms_feature,
      molecular_formula, monoisotopic_mass,
      fims_intensity, lcms_intensity,
      fims_intensity_log2, lcms_intensity_log2
    )
  
  return(list(matched_data = matched_data, venn = venn_result))
}

# ── Execute Import and Matching for both POS and NEG Modes ────────────────────
mode_configs <- list(
  pos = list(fims = fims_pos_path, lcms = lcms_pos_path),
  neg = list(fims = fims_neg_path, lcms = lcms_neg_path)
)

matched_list <- list()
venn_list    <- list()

for (m in names(mode_configs)) {
  res <- process_single_mode(
    mode_name = m,
    fims_path = mode_configs[[m]]$fims,
    lcms_path = mode_configs[[m]]$lcms,
    df_harm   = df_harm
  )
  if (!is.null(res)) {
    if (!is.null(res$matched_data)) matched_list[[m]] <- res$matched_data
    if (!is.null(res$venn))         venn_list[[m]]    <- res$venn
  }
}

if (length(matched_list) == 0) {
  stop("No datasets were successfully processed for either mode.")
}

matched_data <- bind_rows(matched_list)

# Export Raw and Log2 Datasets
matched_data_raw <- matched_data %>% select(-fims_intensity_log2, -lcms_intensity_log2)
write.csv(matched_data_raw, file.path(output_dir, "matched_raw_intensities_MSKohorte.csv"), row.names = FALSE)

matched_data_log2 <- matched_data %>% select(-fims_intensity, -lcms_intensity)
write.csv(matched_data_log2, file.path(output_dir, "matched_log2_intensities_MSKohorte.csv"), row.names = FALSE)
message("✅ Success! Raw and Log2 matched datasets exported.")

## ==============================================================================
## CROSS-PLATFORM METHOD COMPARISON (Log2, Deming, CCC, Bland-Altman, Diagnostics)
## ==============================================================================
message("\n🚀 Running Method Comparison Analysis...")

log2_transform_grouped <- function(x) {
  min_val <- suppressWarnings(min(x, na.rm = TRUE))
  shift   <- if (is.finite(min_val) && min_val <= 0) abs(min_val) else 0
  offset  <- max(shift, 1e-6) + shift * 0.1
  log2(x + offset)
}

concordance_data <- matched_data_raw %>%
  filter(!is.na(fims_intensity), !is.na(lcms_intensity)) %>%
  group_by(ion_mode, fims_feature) %>%
  mutate(fims_log2 = log2_transform_grouped(fims_intensity)) %>%
  ungroup() %>%
  group_by(ion_mode, lcms_feature) %>%
  mutate(lcms_log2 = log2_transform_grouped(lcms_intensity)) %>%
  ungroup() %>%
  filter(Visit %in% c("T1", "T2", "T3", "T4")) %>% 
  mutate(
    ba_mean = (lcms_log2 + fims_log2) / 2,
    ba_diff = lcms_log2 - fims_log2
  )

write.csv(concordance_data, file.path(output_dir, "concordance_data_transformed_T1_to_T4.csv"), row.names = FALSE)

# ── Per-Metabolite Method Comparison Statistics ──────────────────────────────
per_metabolite_stats <- concordance_data %>%
  group_by(ion_mode, harmonized_annotation, fims_feature, lcms_feature) %>%
  filter(n() >= 5, sd(fims_log2, na.rm = TRUE) > 0, sd(lcms_log2, na.rm = TRUE) > 0) %>%
  group_modify(~{
    r_val  <- suppressWarnings(cor(.x$fims_log2, .x$lcms_log2, method = "pearson"))
    r2_val <- r_val^2
    
    dm_fit <- tryCatch({
      mcr::mcreg(.x$fims_log2, .x$lcms_log2, method.reg = "Deming", error.ratio = 1)
    }, error = function(e) NULL)
    
    dm_slope     <- if(!is.null(dm_fit)) dm_fit@para["Slope", "EST"] else NA
    dm_intercept <- if(!is.null(dm_fit)) dm_fit@para["Intercept", "EST"] else NA
    
    ccc_res <- tryCatch({
      DescTools::CCC(.x$fims_log2, .x$lcms_log2, ci = "z-transform", conf.level = 0.95)
    }, error = function(e) NULL)
    
    ccc_val <- if(!is.null(ccc_res)) ccc_res$rho.c$est else NA
    
    mean_bias <- mean(.x$ba_diff, na.rm = TRUE)
    sd_diff   <- sd(.x$ba_diff, na.rm = TRUE)
    
    tibble(
      n                = nrow(.x),
      pearson_r        = r_val,
      r_squared        = r2_val,
      deming_slope     = dm_slope,
      deming_intercept = dm_intercept,
      lins_ccc         = ccc_val,
      mean_bias        = mean_bias,
      sd_diff          = sd_diff,
      lower_loa        = mean_bias - (1.96 * sd_diff),
      upper_loa        = mean_bias + (1.96 * sd_diff)
    )
  }) %>%
  ungroup() %>%
  arrange(desc(r_squared))

write.csv(per_metabolite_stats, file.path(output_dir, "fims_vs_lcms_per_metabolite_method_comparison_T1_to_T4.csv"), row.names = FALSE)

# ── Individual Plots Loop & Combined 5-Column Graphics ───────────────────────
# Load patchwork for assembling the combined plots

format_scientific <- function(x) {
  ifelse(x == 0, "0", sprintf("%.1E", x))
}


if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

unique_metabs <- per_metabolite_stats %>% distinct(ion_mode, harmonized_annotation, fims_feature, lcms_feature)

# Initialize lists to store scatter plots for the combined graphic
scatter_list_pos <- list()
scatter_list_neg <- list()

for(i in seq_len(nrow(unique_metabs))) {
  target <- unique_metabs[i, ]
  
  plot_df <- concordance_data %>% filter(ion_mode == target$ion_mode, fims_feature == target$fims_feature, lcms_feature == target$lcms_feature)
  stats_df <- per_metabolite_stats %>% filter(ion_mode == target$ion_mode, fims_feature == target$fims_feature, lcms_feature == target$lcms_feature)
  
  if(nrow(plot_df) == 0 || nrow(stats_df) == 0) next
  
  sanitized_name <- str_replace_all(stats_df$harmonized_annotation, "[^[:alnum:]]", "_")
  safe_ccc   <- if(is.na(stats_df$lins_ccc)) "NA" else round(stats_df$lins_ccc, 2)
  
  # Bland-Altman
  p_ba <- ggplot(plot_df, aes(x = ba_mean, y = ba_diff)) +
    geom_point(alpha = 0.6, color = "#1B9E77", size = 2) +
    geom_hline(yintercept = stats_df$mean_bias, color = "darkred", linewidth = 0.8) +
    geom_hline(yintercept = stats_df$upper_loa, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = stats_df$lower_loa, linetype = "dashed", color = "grey40") +
    geom_smooth(method = "lm", se = FALSE, color = "darkorange", linetype = "twodash", linewidth = 1) +
    labs(
      title = paste("Bland-Altman:", stats_df$harmonized_annotation),
      subtitle = paste0("Mode: ", toupper(stats_df$ion_mode), " | CCC: ", safe_ccc),
      x = "Mean Concentration [log2]", y = "Difference (LC-MS - FI-MS) [log2]"
    ) + theme_bw()
  ggsave(file.path(ba_individual_dir, paste0("BA_", stats_df$ion_mode, "_", sanitized_name, ".png")), p_ba, width = 6, height = 4.5, dpi = 200)
  
  # Scatter
  p_scatter <- ggplot(plot_df, aes(x = fims_log2, y = lcms_log2)) +
    geom_point(alpha = 0.6, color = "#1B9E77", size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.6) +
    labs(
      title = paste("", stats_df$harmonized_annotation),
      subtitle = paste0("Mode: ", toupper(stats_df$ion_mode), " | R\u00b2 = ", round(stats_df$r_squared, 3)),
      x = "FI-MS Intensity [log2]", y = "LC-MS Intensity [log2]"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title       = element_text(size = 7, face = "bold", hjust = 0.5, lineheight = 0.75),
      plot.subtitle    = element_text(size = 6, hjust = 0.5),
      axis.text        = element_text(size = 6, color = "grey20"),
      axis.title       = element_text(size = 6, face = "bold"),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey85"),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "grey60", linewidth = 0.6),
      axis.ticks       = element_line(linewidth = 0.4),
      plot.margin      = margin(15, 5, 3, 3)
    )
  
  # Save individual scatter
  ggsave(file.path(scatter_dir, paste0("Scatter_", stats_df$ion_mode, "_", sanitized_name, ".png")),
         p_scatter, width = 2.3, height = 2.0, units = "in", dpi = 600)
  
  # Append to the corresponding list for the combined 5-column export
  if (stats_df$ion_mode == "pos") {
    scatter_list_pos[[sanitized_name]] <- p_scatter
  } else if (stats_df$ion_mode == "neg") {
    scatter_list_neg[[sanitized_name]] <- p_scatter
  }
}

# ── Export Combined Scatter Plots (5 Columns per Mode) ──────────────────────
message("🔵 Exporting combined 5-column scatter plots per mode...")

# Define grid parameters based on your individual plot dimensions
plot_width  <- 2.3
plot_height <- 2.0
cols        <- 5

if (length(scatter_list_pos) > 0) {
  rows_pos <- ceiling(length(scatter_list_pos) / cols)
  p_combined_pos <- patchwork::wrap_plots(scatter_list_pos, ncol = cols)
  
  # limitsize = FALSE allows saving very tall files if you have >250 metabolites
  ggsave(file.path(scatter_dir, "Scatter_Combined_POS.png"), p_combined_pos, 
         width = cols * plot_width, height = rows_pos * plot_height, units = "in", 
         dpi = 600, limitsize = FALSE)
}

if (length(scatter_list_neg) > 0) {
  rows_neg <- ceiling(length(scatter_list_neg) / cols)
  p_combined_neg <- patchwork::wrap_plots(scatter_list_neg, ncol = cols)
  
  ggsave(file.path(scatter_dir, "Scatter_Combined_NEG.png"), p_combined_neg, 
         width = cols * plot_width, height = rows_neg * plot_height, units = "in", 
         dpi = 600, limitsize = FALSE)
}

# ── Global Plots & Diagnostics ───────────────────────────────────────────────
# Global Bland-Altman Per Mode
for (current_mode in unique(concordance_data$ion_mode)) {
  mode_data <- concordance_data %>% filter(ion_mode == current_mode)
  if(nrow(mode_data) == 0) next
  
  global_bias <- mean(mode_data$ba_diff, na.rm = TRUE)
  global_sd   <- sd(mode_data$ba_diff, na.rm = TRUE)
  
  p_global_ba <- ggplot(mode_data, aes(x = ba_mean, y = ba_diff)) +
    geom_point(alpha = 0.2, size = 1, aes(color = harmonized_annotation)) +
    geom_hline(yintercept = global_bias, color = "darkred", linewidth = 1) +
    geom_hline(yintercept = global_bias + (1.96 * global_sd), linetype = "dotdash", color = "darkred") +
    geom_hline(yintercept = global_bias - (1.96 * global_sd), linetype = "dotdash", color = "darkred") +
    theme_bw() + theme(legend.position = "none") + 
    labs(title = paste0("Global Bland-Altman — ", toupper(current_mode)), x = "Mean [log2]", y = "Difference [log2]")
  ggsave(file.path(output_dir, paste0("fims_vs_lcms_global_bland_altman_", current_mode, ".png")), p_global_ba, width = 7, height = 5.5, dpi = 300)
}

# Agreement Volcano Plot
p_volcano <- ggplot(per_metabolite_stats %>% filter(!is.na(lins_ccc)), 
                    aes(x = lins_ccc, y = mean_bias, fill = ion_mode, size = r_squared)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = c(0.7, 0.9), linetype = "dashed", color = "grey50") +
  geom_point(shape = 21, color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("pos" = "#7570B3", "neg" = "#1B9E77")) +
  theme_bw() + labs(title = "Global Metabolite Agreement Volcano", x = "Lin's CCC", y = "Mean Bias [log2]")
ggsave(file.path(output_dir, "global_agreement_volcano.png"), p_volcano, width = 8, height = 6, dpi = 300)

# Longitudinal Delta Analysis (T2 - T1)
delta_data <- concordance_data %>%
  filter(Visit %in% c("T1", "T2")) %>%
  select(PseudoID, ion_mode, harmonized_annotation, fims_feature, lcms_feature, Visit, fims_log2, lcms_log2) %>%
  pivot_wider(names_from = Visit, values_from = c(fims_log2, lcms_log2)) %>%
  filter(!is.na(fims_log2_T1), !is.na(fims_log2_T2), !is.na(lcms_log2_T1), !is.na(lcms_log2_T2)) %>%
  mutate(delta_fims = fims_log2_T2 - fims_log2_T1, delta_lcms = lcms_log2_T2 - lcms_log2_T1)

longitudinal_stats <- delta_data %>%
  group_by(ion_mode, harmonized_annotation, fims_feature, lcms_feature) %>%
  filter(n() >= 5) %>%
  summarize(patients_with_t1_t2 = n(), delta_pearson_r = tryCatch(cor(delta_fims, delta_lcms, method = "pearson"), error = function(e) NA), .groups = "drop") %>%
  arrange(desc(delta_pearson_r))
write.csv(longitudinal_stats, file.path(output_dir, "longitudinal_delta_concordance_T2_vs_T1.csv"), row.names = FALSE)

## ==============================================================================
## FEATURE-OVERLAP VENN DIAGRAMS (FIMS vs LC-MS), per mode + combined
## ==============================================================================
message("\n🔵 Building FIMS vs LC-MS feature-overlap Venn diagrams...")

# Counts how many retained metabolites are unique to FIMS, unique to LC-MS,
# or matched (harmonized) between both platforms, for a given pair of
# identity vectors (see fims_identity/lcms_identity construction above).
summarize_overlap <- function(fims_set, lcms_set, label) {
  fims_unique_set <- unique(fims_set)
  lcms_unique_set <- unique(lcms_set)
  matched_ids      <- intersect(fims_unique_set, lcms_unique_set)
  
  tibble(
    comparison = label,
    fims_total = length(fims_unique_set),
    lcms_total = length(lcms_unique_set),
    matched    = length(matched_ids),
    fims_only  = length(fims_unique_set) - length(matched_ids),
    lcms_only  = length(lcms_unique_set) - length(matched_ids)
  )
}

plot_overlap_venn <- function(fims_set, lcms_set, title, out_path) {
  venn_input <- list(`FI-MS` = unique(fims_set), `LC-MS` = unique(lcms_set))
  p <- ggVennDiagram(venn_input, label = "count", edge_size = 1) +
    scale_fill_gradient(low = "#EFF3FF", high = "#1B9E77", name = "Metabolites") +
    scale_x_continuous(expand = expansion(mult = 0.25)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    coord_cartesian(clip = "off") +
    labs(title = title) +
    theme(
      legend.position = "none",
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.margin     = margin(15, 25, 15, 25)
    )
  ggsave(out_path, p, width = 6.5, height = 6, dpi = 300)
}

overlap_stats <- list()

for (m in names(venn_list)) {
  v <- venn_list[[m]]
  overlap_stats[[m]] <- summarize_overlap(v$fims_set, v$lcms_set, toupper(m))
  plot_overlap_venn(
    v$fims_set, v$lcms_set,
    title    = paste0("FI-MS vs LC-MS Metabolite Overlap [", toupper(m), "]"),
    out_path = file.path(venn_dir, paste0("venn_fims_vs_lcms_", m, ".png"))
  )
}

# Combined overview across all processed modes
if (length(venn_list) > 0) {
  combined_fims <- unlist(lapply(venn_list, function(v) v$fims_set), use.names = FALSE)
  combined_lcms <- unlist(lapply(venn_list, function(v) v$lcms_set), use.names = FALSE)
  
  overlap_stats[["combined"]] <- summarize_overlap(combined_fims, combined_lcms, "COMBINED")
  plot_overlap_venn(
    combined_fims, combined_lcms,
    title    = "FI-MS vs LC-MS Metabolite Overlap [All Modes]",
    out_path = file.path(venn_dir, "venn_fims_vs_lcms_combined.png")
  )
}

overlap_stats_df <- bind_rows(overlap_stats)
write.csv(overlap_stats_df, file.path(output_dir, "fims_vs_lcms_feature_overlap_stats.csv"), row.names = FALSE)

message("✅ Venn diagrams saved to: ", normalizePath(venn_dir))
message("✅ Overlap statistics saved to: ", file.path(output_dir, "fims_vs_lcms_feature_overlap_stats.csv"))


## ==============================================================================
## DUAL-AXIS SPAGHETTI PLOTS (per matched metabolite): FI-MS (left) vs LC-MS (right)
## ==============================================================================
message("\n🍝 Building dual-axis FI-MS/LC-MS spaghetti plots per matched metabolite...")

# Scientific notation formatted with 0 decimal places (e.g., 4E+06)
format_scientific <- function(x) {
  ifelse(is.na(x), "", ifelse(x == 0, "0", sprintf("%.0E", x)))
}

# Helper to ensure Row 1 fits on a single line before forcing line break for mode
# Helper: create exactly a two-line title
wrap_annotation <- function(x, ion_mode, max_len = 22) {
  
  words <- strsplit(x, "\\s+")[[1]]
  
  # Fits on one line -> second line is only ion mode
  if (nchar(x) <= max_len) {
    return(paste0(x,
                  "\n(",
                  toupper(ion_mode),
                  ")"))
  }
  
  # Build first line using whole words
  line1 <- words[1]
  
  if (length(words) > 1) {
    for (w in words[-1]) {
      candidate <- paste(line1, w)
      if (nchar(candidate) <= max_len) {
        line1 <- candidate
      } else {
        break
      }
    }
  }
  
  # Everything else goes onto line 2
  line2 <- trimws(sub(paste0("^", line1), "", x))
  
  # Always exactly two rows
  paste0(
    line1,
    "\n",
    line2,
    " (",
    toupper(ion_mode),
    ")"
  )
}

dual_spaghetti_dir <- file.path(output_dir, "dual_axis_spaghetti")
dir.create(dual_spaghetti_dir, showWarnings = FALSE, recursive = TRUE)

rescale_to_range <- function(x, target_min, target_max) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  if (!is.finite(x_min) || !is.finite(x_max) || x_max == x_min) {
    return(rep((target_min + target_max) / 2, length(x)))
  }
  (x - x_min) / (x_max - x_min) * (target_max - target_min) + target_min
}

plot_dual_axis_spaghetti <- function(pair_data, fims_col, lcms_col,
                                     fims_axis_label, lcms_axis_label,
                                     title, out_path, center_fn = mean) {
  pair_data <- pair_data %>% filter(!is.na(.data[[fims_col]]), !is.na(.data[[lcms_col]]))
  if (nrow(pair_data) == 0) return(invisible(NULL))
  
  fims_vals  <- pair_data[[fims_col]]
  lcms_vals  <- pair_data[[lcms_col]]
  fims_range <- range(fims_vals, na.rm = TRUE)
  lcms_range <- range(lcms_vals, na.rm = TRUE)
  
  pair_data <- pair_data %>%
    mutate(
      fims_val    = fims_vals,
      lcms_scaled = rescale_to_range(lcms_vals, fims_range[1], fims_range[2])
    )
  
  fims_mean <- pair_data %>%
    group_by(Visit) %>%
    summarise(mean_val = center_fn(fims_val, na.rm = TRUE), .groups = "drop")
  lcms_mean <- pair_data %>%
    group_by(Visit) %>%
    summarise(mean_val = center_fn(lcms_scaled, na.rm = TRUE), .groups = "drop")
  
  inverse_rescale <- function(y) {
    (y - fims_range[1]) / (fims_range[2] - fims_range[1]) * diff(lcms_range) + lcms_range[1]
  }
  
  p <- ggplot(pair_data, aes(x = Visit)) +
    geom_line(aes(y = fims_val,    group = PseudoID), color = "#1B9E77", alpha = 0.2, linewidth = 0.4) +
    geom_line(aes(y = lcms_scaled, group = PseudoID), color = "#7570B3", alpha = 0.2, linewidth = 0.4) +
    geom_line(data = fims_mean, aes(y = mean_val, group = 1), color = "#1B9E77", linewidth = 1.2) +
    geom_point(data = fims_mean, aes(y = mean_val), color = "#1B9E77", size = 3, shape = 21, fill = "white", stroke = 1.3) +
    geom_line(data = lcms_mean, aes(y = mean_val, group = 1), color = "#7570B3", linewidth = 1.2, linetype = "dashed") +
    geom_point(data = lcms_mean, aes(y = mean_val), color = "#7570B3", size = 3, shape = 21, fill = "white", stroke = 1.3) +
    scale_y_continuous(
      name     = fims_axis_label,
      labels   = format_scientific,
      sec.axis = sec_axis(trans = ~ inverse_rescale(.), name = lcms_axis_label, labels = format_scientific),
      n.breaks = 4
    ) +
    labs(title = title, x = "Visit") +
    theme_bw(base_size = 9) +
    theme(
      axis.title.x       = element_text(size = 9, face = "bold", margin = margin(t = 2)),
      axis.title.y.left  = element_text(color = "#1B9E77", face = "bold", size = 9, margin = margin(r = 2)),
      axis.text.y.left   = element_text(color = "#1B9E77", size = 8),
      axis.title.y.right = element_text(color = "#7570B3", face = "bold", size = 9, margin = margin(l = 2)),
      axis.text.y.right  = element_text(color = "#7570B3", size = 8),
      axis.text.x        = element_text(size = 8, color = "grey20"),
      axis.ticks         = element_line(linewidth = 0.25),
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 8, lineheight = 0.85, margin = margin(b = 3)),
      panel.grid.major   = element_line(linewidth = 0.2, color = "grey88"),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "grey60", linewidth = 0.4),
      plot.margin        = margin(3, 3, 3, 3)
    )
  
  # Canvas footprint: 2.3 x 2.0 inches
  ggsave(out_path, p, width = 2.3, height = 2.0, units = "in", dpi = 600)
}
dual_axis_visits  <- c("T1", "T2", "T3", "T4")

# ── Raw-intensity version ────────────────────────────────────────────────────
dual_axis_targets <- matched_data_raw %>%
  filter(Visit %in% dual_axis_visits) %>%
  distinct(ion_mode, harmonized_annotation, fims_feature, lcms_feature)

for (i in seq_len(nrow(dual_axis_targets))) {
  target <- dual_axis_targets[i, ]
  
  pair_data <- matched_data_raw %>%
    filter(ion_mode == target$ion_mode,
           fims_feature == target$fims_feature,
           lcms_feature == target$lcms_feature,
           Visit %in% dual_axis_visits) %>%
    mutate(Visit = factor(Visit, levels = dual_axis_visits, ordered = TRUE))
  
  if (nrow(pair_data) == 0) next
  
  sanitized_name <- str_replace_all(target$harmonized_annotation, "[^[:alnum:]]", "_")
  out_filename   <- paste0("DualSpaghetti_", target$ion_mode, "_", sanitized_name, ".png")
  plot_title <- wrap_annotation(
    target$harmonized_annotation,
    target$ion_mode
  )
  
  plot_dual_axis_spaghetti(
    pair_data,
    fims_col        = "fims_intensity",
    lcms_col        = "lcms_intensity",
    fims_axis_label = "FI-MS Intensity",
    lcms_axis_label = "LC-MS Intensity",
    title           = plot_title,
    out_path        = file.path(dual_spaghetti_dir, out_filename),
    center_fn       = median
  )
}

message("✅ Dual-axis spaghetti plots (raw intensity, median) saved to: ", normalizePath(dual_spaghetti_dir))

# ── Raw-intensity version, mean center line ──────────────────────────────────
dual_spaghetti_mean_dir <- file.path(output_dir, "dual_axis_spaghetti_mean")
dir.create(dual_spaghetti_mean_dir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(nrow(dual_axis_targets))) {
  target <- dual_axis_targets[i, ]
  
  pair_data <- matched_data_raw %>%
    filter(ion_mode == target$ion_mode,
           fims_feature == target$fims_feature,
           lcms_feature == target$lcms_feature,
           Visit %in% dual_axis_visits) %>%
    mutate(Visit = factor(Visit, levels = dual_axis_visits, ordered = TRUE))
  
  if (nrow(pair_data) == 0) next
  
  sanitized_name <- str_replace_all(target$harmonized_annotation, "[^[:alnum:]]", "_")
  out_filename   <- paste0("DualSpaghetti_Mean_", target$ion_mode, "_", sanitized_name, ".png")
  plot_title <- wrap_annotation(
    target$harmonized_annotation,
    target$ion_mode
  )
  
  plot_dual_axis_spaghetti(
    pair_data,
    fims_col        = "fims_intensity",
    lcms_col        = "lcms_intensity",
    fims_axis_label = "FI-MS Intensity",
    lcms_axis_label = "LC-MS Intensity",
    title           = plot_title,
    out_path        = file.path(dual_spaghetti_mean_dir, out_filename),
    center_fn       = mean
  )
}

message("✅ Dual-axis spaghetti plots (raw intensity, mean) saved to: ", normalizePath(dual_spaghetti_mean_dir))

# ── Log2-transformed version ────────────────────────────────────────────────
dual_spaghetti_log2_dir <- file.path(output_dir, "dual_axis_spaghetti_log2")
dir.create(dual_spaghetti_log2_dir, showWarnings = FALSE, recursive = TRUE)

dual_axis_targets_log2 <- concordance_data %>%
  distinct(ion_mode, harmonized_annotation, fims_feature, lcms_feature)

for (i in seq_len(nrow(dual_axis_targets_log2))) {
  target <- dual_axis_targets_log2[i, ]
  
  pair_data <- concordance_data %>%
    filter(ion_mode == target$ion_mode,
           fims_feature == target$fims_feature,
           lcms_feature == target$lcms_feature) %>%
    mutate(Visit = factor(Visit, levels = dual_axis_visits, ordered = TRUE))
  
  if (nrow(pair_data) == 0) next
  
  sanitized_name <- str_replace_all(target$harmonized_annotation, "[^[:alnum:]]", "_")
  out_filename   <- paste0("DualSpaghetti_Log2_", target$ion_mode, "_", sanitized_name, ".png")
  plot_title <- wrap_annotation(
    target$harmonized_annotation,
    target$ion_mode
  )
  
  plot_dual_axis_spaghetti(
    pair_data,
    fims_col        = "fims_log2",
    lcms_col        = "lcms_log2",
    fims_axis_label = "FI-MS [log2]",
    lcms_axis_label = "LC-MS [log2]",
    title           = plot_title,
    out_path        = file.path(dual_spaghetti_log2_dir, out_filename)
  )
}

message("✅ Dual-axis spaghetti plots (log2 intensity, mean) saved to: ", normalizePath(dual_spaghetti_log2_dir))

# ── Log2-transformed version, median center line ────────────────────────────
dual_spaghetti_log2_median_dir <- file.path(output_dir, "dual_axis_spaghetti_log2_median")
dir.create(dual_spaghetti_log2_median_dir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(nrow(dual_axis_targets_log2))) {
  target <- dual_axis_targets_log2[i, ]
  
  pair_data <- concordance_data %>%
    filter(ion_mode == target$ion_mode,
           fims_feature == target$fims_feature,
           lcms_feature == target$lcms_feature) %>%
    mutate(Visit = factor(Visit, levels = dual_axis_visits, ordered = TRUE))
  
  if (nrow(pair_data) == 0) next
  
  sanitized_name <- str_replace_all(target$harmonized_annotation, "[^[:alnum:]]", "_")
  out_filename   <- paste0("DualSpaghetti_Log2_Median_", target$ion_mode, "_", sanitized_name, ".png")
  plot_title <- wrap_annotation(
    target$harmonized_annotation,
    target$ion_mode
  )
  
  plot_dual_axis_spaghetti(
    pair_data,
    fims_col        = "fims_log2",
    lcms_col        = "lcms_log2",
    fims_axis_label = "FI-MS [log2]",
    lcms_axis_label = "LC-MS [log2]",
    title           = plot_title,
    out_path        = file.path(dual_spaghetti_log2_median_dir, out_filename),
    center_fn       = median
  )
}

message("✅ Dual-axis spaghetti plots (log2 intensity, median) saved to: ", normalizePath(dual_spaghetti_log2_median_dir))

message("\n✨ Pipeline complete! All output files written to: ", normalizePath(output_dir))


## ==============================================================================
## EXTRA OUTPUT: RAW-INTENSITY MEDIAN DUAL SPAGHETTI, LMM-SIGNIFICANT ONLY
## ==============================================================================
message("\n🧬 [EXTRA] Loading LMM significance results (FI-MS + LC-MS)... Needs all FIMS and LCMS data, maybe run whole script")

lmm_results_paths <- list(
  pos = list(
    fims = file.path(dirname(fims_pos_path), "LMM_Visit_Results_POS.csv"),
    lcms = file.path(dirname(lcms_pos_path), "LMM_Visit_Results_LCMS_POS.csv")
  ),
  neg = list(
    fims = file.path(dirname(fims_neg_path), "LMM_Visit_Results_NEG.csv"),
    lcms = file.path(dirname(lcms_neg_path), "LMM_Visit_Results_LCMS_NEG.csv")
  )
)

fdr_threshold <- 0.05

load_sig_metabolites <- function(path) {
  message("    Checking: ", path)
  if (!file.exists(path)) {
    warning("   ❌ File NOT FOUND — treating as zero significant metabolites: ", path)
    return(character(0))
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("Metabolite", "FDR") %in% colnames(df))) {
    warning("   ❌ Missing 'Metabolite' or 'FDR' column in: ", path,
            " (found: ", paste(colnames(df), collapse = ", "), ")")
    return(character(0))
  }
  sig <- df %>% dplyr::filter(FDR < fdr_threshold) %>% dplyr::pull(Metabolite)
  message("   ✅ Found ", nrow(df), " tested, ", length(sig), " significant (FDR < ", fdr_threshold, ")")
  sig
}

sig_metabolites <- list()
for (m in names(lmm_results_paths)) {
  message("  [", toupper(m), "] FI-MS results:")
  fims_sig <- load_sig_metabolites(lmm_results_paths[[m]]$fims)
  message("  [", toupper(m), "] LC-MS results:")
  lcms_sig <- load_sig_metabolites(lmm_results_paths[[m]]$lcms)
  sig_metabolites[[m]] <- list(fims = fims_sig, lcms = lcms_sig)
}

tag_lmm_significance <- function(targets_df) {
  targets_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      sig_fims   = ion_mode %in% names(sig_metabolites) && fims_feature %in% sig_metabolites[[ion_mode]]$fims,
      sig_lcms   = ion_mode %in% names(sig_metabolites) && lcms_feature %in% sig_metabolites[[ion_mode]]$lcms,
      sig_either = sig_fims | sig_lcms
    ) %>%
    dplyr::ungroup()
}

dual_axis_targets_lmm <- matched_data_raw %>%
  filter(Visit %in% dual_axis_visits) %>%
  distinct(ion_mode, harmonized_annotation, fims_feature, lcms_feature) %>%
  tag_lmm_significance()

write.csv(
  dual_axis_targets_lmm,
  file.path(output_dir, "dual_axis_targets_lmm_significance.csv"), row.names = FALSE
)
message("  [SIG SUMMARY] Matched pairs total: ", nrow(dual_axis_targets_lmm))
message("  [SIG SUMMARY] Significant in FI-MS: ", sum(dual_axis_targets_lmm$sig_fims))
message("  [SIG SUMMARY] Significant in either platform: ", sum(dual_axis_targets_lmm$sig_either))

dual_axis_targets_lmm <- dual_axis_targets_lmm %>% filter(sig_fims)

# ── Output folder: separate from the unfiltered spaghetti folders ──────────
dual_spaghetti_lmm_dir <- file.path(output_dir, "dual_axis_spaghetti_LMM_significant")
dir.create(dual_spaghetti_lmm_dir, showWarnings = FALSE, recursive = TRUE)

if (nrow(dual_axis_targets_lmm) == 0) {
  warning("No FI-MS-significant matched metabolite pairs found — no plots generated. ",
          "Check dual_axis_targets_lmm_significance.csv and the file-not-found warnings above.")
} else {
  for (i in seq_len(nrow(dual_axis_targets_lmm))) {
    target <- dual_axis_targets_lmm[i, ]
    
    pair_data <- matched_data_raw %>%
      filter(ion_mode == target$ion_mode,
             fims_feature == target$fims_feature,
             lcms_feature == target$lcms_feature,
             Visit %in% dual_axis_visits) %>%
      mutate(Visit = factor(Visit, levels = dual_axis_visits, ordered = TRUE))
    
    if (nrow(pair_data) == 0) next
    
    sanitized_name <- str_replace_all(target$harmonized_annotation, "[^[:alnum:]]", "_")
    out_filename   <- paste0(target$ion_mode, "_", sanitized_name, ".png")
    plot_title <- wrap_annotation(
      target$harmonized_annotation,
      target$ion_mode
    )
    
    plot_dual_axis_spaghetti(
      pair_data,
      fims_col        = "fims_intensity",
      lcms_col        = "lcms_intensity",
      fims_axis_label = "FI-MS Intensity",
      lcms_axis_label = "LC-MS Intensity",
      title           = plot_title,
      out_path        = file.path(dual_spaghetti_lmm_dir, out_filename),
      center_fn       = median
    )
  }
  message("✅ [EXTRA] FI-MS-significant raw-intensity (median) dual spaghetti plots saved to: ",
          normalizePath(dual_spaghetti_lmm_dir))
}