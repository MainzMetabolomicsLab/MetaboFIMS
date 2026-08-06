#####
##### Raw data import and data rearrangement of LINEAR REGRESSION FILES
#####
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
library(stringr)


if (.Platform$OS.type == "unix") {
  base_ms_path   <- "/Volumes/T7/Arbeit/FIA/EmDia_NEU/10ppm annotation"
  metadata_path  <- "/Users/fabianschmitt/Seafile/Meine Bibliothek/FIA/EmDia_MS_Kohorte_FIMS_2025 batch and metadata/AG_Schmidlin_20251120_EmDia_V3.csv"
} else {
  base_ms_path   <- "D:/Arbeit/FIA/EmDia_NEU/10ppm annotation"
  metadata_path  <- "C:/Users/fabia/Seafile/Meine Bibliothek/FIA/EmDia_MS_Kohorte_FIMS_2025 batch and metadata/AG_Schmidlin_20251120_EmDia_V3.csv"
}

modes <- list(
  Positive = file.path(base_ms_path, "CSV/feature_matrix_transposed_EmDia_P1_5_with_metadata.csv"),
  Negative = file.path(base_ms_path, "CSV_NEG/feature_matrix_transposed_EmDia_P1_5_with_metadata_neg.csv")
)



for (mode_name in names(modes)) {
  
  emdia_data_path <- modes[[mode_name]]
  
  
  if (!file.exists(emdia_data_path)) {
    warning("❌ File not found, skipping: ", emdia_data_path)
    next
  }
  
  emdia_data_raw <- read.csv2(emdia_data_path, header = TRUE, sep = ",", dec = ".")
  
  existing_meta <- emdia_data_raw[1:3, ]       
  emdia_data    <- emdia_data_raw[-c(1:3), ]    
  
  sample_cols <- colnames(emdia_data)[-1]       
  
  tube_codes <- str_extract(sample_cols, "_\\d{9,}$")
  tube_codes <- sub("^_", "", tube_codes)       
  
  sample_info <- data.frame(
    column_name = sample_cols,
    tubeCode = tube_codes,
    stringsAsFactors = FALSE
  )
  
  if (!file.exists(metadata_path)) {
    stop("❌ External metadata file missing at: ", metadata_path)
  }
  emdia_metadata <- read.csv2(metadata_path, header = TRUE, sep = ",", dec = ".")
  
  sample_info$tubeCode    <- as.character(sample_info$tubeCode)
  emdia_metadata$tubeCode <- as.character(emdia_metadata$tubeCode)
  
  sample_info <- sample_info %>%
    left_join(emdia_metadata, by = "tubeCode")
  
  sample_info <- sample_info[match(sample_cols, sample_info$column_name), ]
  
  sex_row   <- data.frame(SampleID = "Sex",   t(sample_info$gender), stringsAsFactors = FALSE)
  group_row <- data.frame(SampleID = "Group", t(sample_info$group),  stringsAsFactors = FALSE)
  visit_row <- data.frame(SampleID = "Visit", t(sample_info$Visit),  stringsAsFactors = FALSE)
  age_row   <- data.frame(SampleID = "Age",   t(sample_info$age),    stringsAsFactors = FALSE)
  
  colnames(sex_row)   <- colnames(emdia_data)
  colnames(group_row) <- colnames(emdia_data)
  colnames(visit_row) <- colnames(emdia_data)
  colnames(age_row)   <- colnames(emdia_data)
  
  meta_keep <- c("class", "index", "batch")
  
  data_part <- emdia_data[!emdia_data$SampleID %in% c(meta_keep, "Sex", "Group", "Visit", "Age"), ]
  

  emdia_data_new <- rbind(
    existing_meta,   
    sex_row,
    group_row,
    visit_row,
    age_row,         
    data_part        
  )
  
  emdia_data_new[is.na(emdia_data_new)] <- ""
  
  unmatched <- sample_info %>% filter(!is.na(tubeCode) & (is.na(Visit) | is.na(group)))
  if (nrow(unmatched) > 0) {
    message("⚠️ Warning: Some tubeCodes did not match metadata in ", mode_name, ":")
    print(unmatched$column_name)
  }
  
  row_headers <- emdia_data_new$SampleID
  emdia_data_transposed <- as.data.frame(t(emdia_data_new[, -1]))
  colnames(emdia_data_transposed) <- row_headers
  
  emdia_data_transposed <- cbind(SampleID = rownames(emdia_data_transposed), emdia_data_transposed)
  rownames(emdia_data_transposed) <- NULL
  
  emdia_data_dir          <- dirname(emdia_data_path)
  emdia_data_filename     <- basename(emdia_data_path)
  emdia_data_new_filename <- sub("\\.csv$", "_with_age.csv", emdia_data_filename)
  emdia_data_new_path     <- file.path(emdia_data_dir, emdia_data_new_filename)
  
  write.csv(emdia_data_transposed, file = emdia_data_new_path, row.names = FALSE, quote = TRUE)
  
}



library(dplyr)
library(stringr)



file_paths <- list(
  pos = file.path(base_ms_path, "CSV/feature_matrix_transposed_EmDia_P1_5_with_metadata_with_age.csv"),
  neg = file.path(base_ms_path, "CSV_NEG/feature_matrix_transposed_EmDia_P1_5_with_metadata_neg_with_age.csv")
)

meta_cols <- c("class", "index", "batch", "Sex", "Group", "Visit", "Age")


for (mode in names(file_paths)) {
  
  input_file <- file_paths[[mode]]
  
  
  if (!file.exists(input_file)) {
    warning("File not found, skipping: ", input_file)
    next
  }
  

  emdia_data_transposed <- read.csv(
    input_file,
    header           = TRUE,
    sep              = ",",
    dec              = ".",
    check.names      = FALSE,   
    stringsAsFactors = FALSE,
    colClasses       = "character"
  )
  
  rownames(emdia_data_transposed) <- emdia_data_transposed[, 1]
  emdia_data_transposed           <- emdia_data_transposed[, -1, drop = FALSE]
  
  colnames(emdia_data_transposed) <- trimws(colnames(emdia_data_transposed))
  
  feature_cols <- setdiff(colnames(emdia_data_transposed), meta_cols)
  message("Initial feature space count: ", length(feature_cols))
  
  qc_data    <- emdia_data_transposed %>% filter(class == "QC")
  blank_data <- emdia_data_transposed %>% filter(class == "Blank_processing")
  
  if (nrow(qc_data) == 0) {
    warning("No 'QC' samples found in the 'class' column. Skipping mode: ", toupper(mode))
    next
  }
  
  qc_cvs <- sapply(feature_cols, function(f) {
    values <- as.numeric(qc_data[[f]])
    m <- mean(values, na.rm = TRUE)
    s <- sd(values,   na.rm = TRUE)
    if (is.na(m) || m == 0) return(NA_real_)
    (s / m) * 100
  })
  
  features_pass_cv <- names(qc_cvs)[!is.na(qc_cvs) & qc_cvs < 25]
  
  
  qc_means    <- sapply(feature_cols, function(f) mean(as.numeric(qc_data[[f]]),    na.rm = TRUE))
  blank_means <- sapply(feature_cols, function(f) mean(as.numeric(blank_data[[f]]), na.rm = TRUE))
  blank_means[blank_means == 0 | is.na(blank_means)] <- NA_real_
  
  qc_to_blank_ratio   <- qc_means / blank_means
  features_pass_blank <- names(qc_to_blank_ratio)[!is.na(qc_to_blank_ratio) & qc_to_blank_ratio >= 4]
  
  
  features_to_keep <- intersect(features_pass_cv, features_pass_blank)
  
  
  valid_cols          <- intersect(c(meta_cols, features_to_keep), colnames(emdia_data_transposed))
  emdia_data_filtered <- emdia_data_transposed[, valid_cols, drop = FALSE]
  
  emdia_data_filtered[is.na(emdia_data_filtered) | emdia_data_filtered == "NA"] <- ""
  
  
  
  input_basename <- tools::file_path_sans_ext(basename(input_file))
  out_file       <- file.path(dirname(input_file),
                              paste0(input_basename, "_CV_Blank_filtered.csv"))
  

  write.table(
    cbind(SampleID = rownames(emdia_data_filtered), emdia_data_filtered),
    file      = out_file,
    sep       = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote     = TRUE
  )
  
  
}








###
###
###Load in data
###
###
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

is_unix <- .Platform$OS.type == "unix"

file_paths <- list(
  pos = if (is_unix) {
    "/Volumes/T7/Arbeit/FIA/EmDia_NEU/CSV/feature_matrix_transposed_EmDia_P1_5_with_metadata_with_age_CV_Blank_filtered.csv"
  } else {
    "D:/Arbeit/FIA/EmDia_NEU/CSV/feature_matrix_transposed_EmDia_P1_5_with_metadata_with_age_CV_Blank_filtered.csv"
  },
  neg = if (is_unix) {
    "/Volumes/T7/Arbeit/FIA/EmDia_NEU/CSV_NEG/feature_matrix_transposed_EmDia_P1_5_with_metadata_neg_with_age_CV_Blank_filtered.csv"
  } else {
    "D:/Arbeit/FIA/EmDia_NEU/CSV_NEG/feature_matrix_transposed_EmDia_P1_5_with_metadata_neg_with_age_CV_Blank_filtered.csv"
  }
)

file_paths <- list(
  pos = file.path(base_ms_path, "CSV/feature_matrix_transposed_EmDia_P1_5_with_metadata_with_age_CV_Blank_filtered.csv"),
  neg = file.path(base_ms_path, "CSV_NEG/feature_matrix_transposed_EmDia_P1_5_with_metadata_neg_with_age_CV_Blank_filtered.csv")
)

meta_cols <- c("class", "index", "batch", "Sex", "Group", "Visit", "Age")

target_metabolites <- list(
  pos = c(
    "Uric acid (169.0349;M+H)",
    "1,5-Anhydroglucitol (187.0571;M+Na)"
  ),
  neg = c(
    "Uric acid (167.0203;M-H)",
    "1,5-Anhydroglucitol (209.067;M+HCOO-)"
  )
)

x_axis_labels <- list(
  pos = c(
    "Uric acid (169.0349;M+H)"            = "Uric acid\n(169.0349; M+H)",
    "1,5-Anhydroglucitol (187.0571;M+Na)" = "1,5-Anhydroglucitol\n(187.0571; M+Na)"
  ),
  neg = c(
    "Uric acid (167.0203;M-H)"                  = "Uric acid\n(167.0203; M-H)",
    "1,5-Anhydroglucitol (209.067;M+HCOO-)" = "1,5-Anhydroglucitol\n(209.067; M+HCOO-)"
  )
)

meta_cols <- c("class", "index", "batch", "Sex", "Group", "Visit", "Age")

for (mode in names(file_paths)) {
  
  input_file <- file_paths[[mode]]
  out_dir    <- dirname(input_file)
  
  message("\n", strrep("─", 60))
  message("▶  Processing EmDia ionisation mode: ", toupper(mode))
  message("   File: ", input_file)
  
  if (!file.exists(input_file)) {
    warning("File not found, skipping: ", input_file)
    next
  }
  
  dat <- read.csv(
    input_file,
    header           = TRUE,
    sep              = ",",
    dec              = ".",
    check.names      = FALSE,
    colClasses       = "character"
  )
  
  rownames(dat) <- dat[, 1]
  dat           <- dat[, -1, drop = FALSE]
  colnames(dat) <- trimws(colnames(dat))
  
  message("   Loaded: ", nrow(dat), " samples × ",
          length(setdiff(colnames(dat), meta_cols)), " features")
  
  pca_features <- dat[, !colnames(dat) %in% meta_cols]
  pca_features <- apply(pca_features, 2, as.numeric)
  
  pca_res  <- prcomp(pca_features, center = TRUE, scale. = TRUE)
  var_expl <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
  pc1_lab  <- paste0("PC1 (", var_expl[1], "%)")
  pc2_lab  <- paste0("PC2 (", var_expl[2], "%)")
  
  pca_df <- data.frame(
    dat[, meta_cols, drop = FALSE],
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2]
  )
  
  n_groups    <- length(unique(pca_df$class))
  brewer_cols <- colorRampPalette(brewer.pal(min(n_groups, 8), "Dark2"))(n_groups)
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = class)) +
    geom_point(size = 1.5, alpha = 0.8) +
    stat_ellipse(type = "norm", level = 0.95, linetype = 1, linewidth = 0.5) +
    scale_color_manual(values = brewer_cols) +
    labs(title = paste0("PCA EmDia [ESI", toupper(mode), "]"),
         x = pc1_lab, y = pc2_lab) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.major  = element_line(color = "grey85"),
      panel.grid.minor  = element_blank(),
      axis.text         = element_text(color = "black", size = 8),
      axis.title        = element_text(size = 8),
      plot.title        = element_blank(),
      legend.title      = element_blank(),
      legend.key        = element_blank(),
      legend.text       = element_text(size = 6),
      legend.key.size   = unit(0.3, "lines"),
      legend.spacing.y  = unit(0.1, "lines"),
      legend.position   = "bottom"
    )
  
  print(p_pca)
  ggsave(file.path(out_dir, paste0("PCA_EmDia_full_data_", mode, ".png")),
         p_pca, width = 2.3, height = 2, dpi = 300)
  
  samples_only        <- dat[dat$class == "Sample", ]
  samples_only$batch  <- factor(samples_only$batch)
  
  pca_feat_s <- samples_only[, !colnames(samples_only) %in% meta_cols]
  pca_feat_s <- apply(pca_feat_s, 2, as.numeric)
  
  pca_res_s  <- prcomp(pca_feat_s, center = TRUE, scale. = TRUE)
  var_expl_s <- round(100 * (pca_res_s$sdev^2 / sum(pca_res_s$sdev^2)), 1)
  pc1_lab_s  <- paste0("PC1 (", var_expl_s[1], "%)")
  pc2_lab_s  <- paste0("PC2 (", var_expl_s[2], "%)")
  
  pca_df_s <- data.frame(
    samples_only[, meta_cols, drop = FALSE],
    PC1 = pca_res_s$x[, 1],
    PC2 = pca_res_s$x[, 2]
  )
  
  batch_levels <- levels(pca_df_s$batch)
  n_batches    <- length(batch_levels)
  batch_cols   <- setNames(
    colorRampPalette(brewer.pal(min(n_batches, 8), "Dark2"))(n_batches),
    batch_levels
  )
  
  p_pca_batch <- ggplot(pca_df_s, aes(x = PC1, y = PC2, color = batch)) +
    geom_point(size = 1.5, alpha = 0.8) +
    stat_ellipse(type = "norm", level = 0.95, linetype = 1, linewidth = 0.5) +
    scale_color_manual(values = batch_cols) +
    labs(x = pc1_lab_s, y = pc2_lab_s) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.major  = element_line(color = "grey85"),
      panel.grid.minor  = element_blank(),
      axis.text         = element_text(color = "black", size = 8),
      axis.title        = element_text(size = 8),
      plot.title        = element_blank(),
      legend.title      = element_blank(),
      legend.key        = element_blank(),
      legend.text       = element_text(size = 6),
      legend.key.size   = unit(0.3, "lines"),
      legend.spacing.y  = unit(0.1, "lines"),
      legend.position   = "right"
    )
  
  print(p_pca_batch)
  ggsave(file.path(out_dir, paste0("PCA_EmDia_samples_by_batch_", mode, ".png")),
         p_pca_batch, width = 2.3, height = 2, dpi = 300)
  
  dat_cv <- dat %>%
    filter(xor(class == "QC", class == "Sample")) %>%
    dplyr::select(-c(index, Sex, Group, Visit, Age)) %>%
    rownames_to_column(var = "SampleID")
  
  data_long_cv <- pivot_longer(dat_cv,
                               cols      = -c(SampleID, batch, class),
                               names_to  = "Metabolite",
                               values_to = "Value") %>%
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
      title    = paste0("CV of Metabolites by Batch and Class (", toupper(mode), ")"),
      subtitle = "Violin and boxplot (CV capped at 200%)",
      x = "Class", y = "Coefficient of Variation (%)", fill = "Batch"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x   = element_text(size = 16, angle = 45, hjust = 1),
      axis.text.y   = element_text(size = 16),
      axis.title.x  = element_text(size = 18, margin = margin(t = 10)),
      axis.title.y  = element_text(size = 18, margin = margin(r = 10)),
      legend.title  = element_text(size = 16),
      legend.text   = element_text(size = 14),
      plot.title    = element_text(size = 20, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 10))
    )
  
  print(p_cv)
  ggsave(file.path(out_dir, paste0("CVs_EmDia_", mode, ".png")),
         p_cv, width = 9, height = 7, dpi = 300)
  
  dat_stats <- dat %>%
    filter(xor(Group == "Empagliflozin", Group == "Placebo")) %>%
    filter(xor(Visit == "V2",            Visit == "V3"))  %>%
    dplyr::select(-c(class, index, batch, Sex, Visit, Age)) %>%
    as.data.frame()
  
  dat_stats[, -1] <- sapply(dat_stats[, -1], as.numeric)
  dat_stats$Group <- as.factor(dat_stats$Group)
  
  feature_cols_stats <- colnames(dat_stats)[-1]
  
  results <- data.frame(Metabolite       = character(),
                        p.value          = numeric(),
                        log2_fold_change = numeric(),
                        stringsAsFactors = FALSE)
  
  for (metabolite in feature_cols_stats) {
    lvls    <- levels(dat_stats$Group)
    raw_g1  <- dat_stats[[metabolite]][dat_stats$Group == lvls[1]]
    raw_g2  <- dat_stats[[metabolite]][dat_stats$Group == lvls[2]]
    
    log_g1  <- log2(raw_g1 + 1)
    log_g2  <- log2(raw_g2 + 1)
    
    tt      <- t.test(log_g2, log_g1)
    
    mean_g1 <- mean(raw_g1, na.rm = TRUE)
    mean_g2 <- mean(raw_g2, na.rm = TRUE)
    
    if (is.na(mean_g1) || is.na(mean_g2) || mean_g1 == 0 || mean_g2 == 0) {
      log2_fc <- NA_real_
    } else {
      log2_fc <- log2(mean_g2 / mean_g1)
    }
    
    results <- rbind(results,
                     data.frame(Metabolite       = metabolite,
                                p.value          = tt$p.value,
                                log2_fold_change = log2_fc))
  }
  
  results <- results %>%
    mutate(
      neg_log10_p  = -log10(p.value),
      p.adj        = p.adjust(p.value, method = "fdr"),
      significance = case_when(
        p.adj < 0.0001 ~ "****",
        p.adj < 0.001  ~ "***",
        p.adj < 0.01   ~ "**",
        p.adj < 0.05   ~ "*",
        TRUE           ~ ""
      )
    )
  
  sig_metabolites <- results %>% filter(p.adj < 0.05) %>% pull(Metabolite)
  message("   Significant metabolites (FDR < 0.05 on Log2 scale): ", length(sig_metabolites))
  
  write.csv(results, file.path(out_dir, paste0("ttest_results_EmDia_Log2Transformed_", mode, ".csv")),
            row.names = FALSE)
  
  targets  <- target_metabolites[[mode]]
  x_labels <- x_axis_labels[[mode]]
  
  missing <- setdiff(targets, colnames(dat_stats))
  if (length(missing) > 0) {
    warning("Target metabolite(s) not found in ", toupper(mode), " data, skipping boxplot: ",
            paste(missing, collapse = ", "))
  } else {
    
    box_df <- dat_stats %>%
      dplyr::select(Group, all_of(targets)) %>%
      pivot_longer(cols = targets, names_to = "Metabolite", values_to = "Intensity") %>%
      mutate(
        Metabolite = factor(Metabolite, levels = targets),
        Intensity  = as.numeric(Intensity)
      ) %>%
      group_by(Metabolite) %>%
      mutate(Relative_Intensity = Intensity / max(Intensity, na.rm = TRUE) * 100) %>%
      ungroup()
    
    annotation_box <- results %>%
      filter(Metabolite %in% targets) %>%
      dplyr::select(Metabolite, significance) %>%
      mutate(
        Metabolite = factor(Metabolite, levels = targets),
        x      = as.numeric(Metabolite),
        xstart = x - 0.25,
        xend   = x + 0.25,
        yline  = 105,
        y      = 110
      )
    
    custom_colors <- c("#1f78b4", "#d95f02")
    
    p_box <- ggplot(box_df, aes(x = Metabolite, y = Relative_Intensity, fill = Group)) +
      geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
                   outlier.shape = NA) +
      scale_fill_manual(values = custom_colors) +
      scale_x_discrete(labels = x_labels) +
      geom_segment(data = annotation_box,
                   aes(x = xstart, xend = xend, y = yline, yend = yline),
                   inherit.aes = FALSE, linewidth = 0.8) +
      geom_text(data = annotation_box,
                aes(x = x, y = y, label = significance),
                inherit.aes = FALSE, vjust = 0.5, size = 3.5) +
      labs(y = "Relative Intensity (%)", x = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      theme_minimal(base_size = 8) +
      theme(
        plot.title     = element_blank(),
        plot.subtitle  = element_blank(),
        axis.title     = element_text(size = 8),
        axis.text.x    = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y    = element_text(size = 8),
        legend.title   = element_text(size = 8),
        legend.text    = element_text(size = 8),
        legend.key.size   = unit(0.3, "lines"),
        legend.spacing    = unit(0.1, "lines"),
        legend.margin     = ggplot2::margin(0, 0, 0, 0),
        plot.margin       = ggplot2::margin(2, 2, 2, 2, unit = "pt")
      )
    
    print(p_box)
    ggsave(file.path(out_dir, paste0("Boxplot_UricAcid_15AG_", mode, ".png")),
           p_box, width = 2.3, height = 3, dpi = 300)
  }
  
  message("✅  Done: ", toupper(mode))
}










###
### MetaboFIMS paper revision Elastic Net analysis
###








library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)


if (.Platform$OS.type == "unix") {
  Elastic_Net <- "/Volumes/T7/Arbeit/FIA/EmDia_NEU/Elastic_Net.csv"
} else {
  Elastic_Net <- "D:/Arbeit/FIA/EmDia_NEU/Elastic_Net.csv"
}

import_dir <- dirname(Elastic_Net)


df <- read.csv(
  Elastic_Net,
  check.names = TRUE,
  stringsAsFactors = FALSE
)

r2_cutoff <- 0.2

df_filtered <- df %>%
  filter(`LC.MS` > r2_cutoff | `FI.MS` > r2_cutoff) %>%
  mutate(
    Max_R2 = pmax(`LC.MS`, `FI.MS`, na.rm = TRUE)
  ) %>%
  arrange(Max_R2) %>%
  mutate(
    Clinical_Trait = factor(
      Clinical_Trait,
      levels = Clinical_Trait
    )
  )

df_long <- df_filtered %>%
  pivot_longer(
    cols = c("LC.MS", "FI.MS"),
    names_to = "Method",
    values_to = "R2"
  )


p <- ggplot(
  df_long,
  aes(
    x = R2,
    y = Clinical_Trait,
    fill = Method
  )
) +

  geom_vline(
    xintercept = r2_cutoff,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +

  geom_col(
    position = position_dodge(width = 0.85), 
    width = 0.55,                            
    colour = "black",
    linewidth = 0.3
  ) +

  geom_text(
    aes(label = sprintf("%.2f", R2)),
    position = position_dodge(width = 0.85), 
    hjust = -0.2,
    size = 2, 
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "LC.MS" = "#7570B3",  
      "FI.MS" = "#1B9E77"   
    ),
    labels = c(
      "LC.MS" = "LC-MS",
      "FI.MS" = "FI-MS"
    )
  ) +

  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = sprintf("%.1f", seq(0, 1, 0.2)),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = expression("10-Fold Cross-Validated "*R^2),
    y = NULL,
    fill = NULL
  ) +

  theme_bw(base_size = 8) +
  theme(

    plot.title         = element_blank(),
    plot.subtitle      = element_blank(),
    axis.title         = element_text(size = 8, color = "black"),

    axis.text.x        = element_text(size = 8, color = "black"), 
    axis.text.y        = element_text(size = 8, color = "black"),
    legend.title       = element_text(size = 8, color = "black"),
    legend.text        = element_text(size = 8, color = "black"),
    legend.key.size    = unit(0.3, "lines"),
    legend.spacing.y   = unit(0.1, "lines"), 
    legend.margin      = ggplot2::margin(0, 0, 0, 0),
    plot.margin        = ggplot2::margin(2, 2, 2, 2, unit = "pt"),
    

    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line          = element_blank(),
    axis.ticks         = element_line(linewidth = 0.4, color = "black"),
    panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.4), 
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position    = "top",
    legend.justification = "center",
    legend.key         = element_blank()
  )

print(p)


plot_height <- max(3, nrow(df_filtered) * 0.20) 


ggsave(
  filename = file.path(import_dir, "ElasticNet_R2_Barplot.png"),
  plot = p,
  width = 4.7, 
  height = plot_height,
  units = "in",
  dpi = 600
)
