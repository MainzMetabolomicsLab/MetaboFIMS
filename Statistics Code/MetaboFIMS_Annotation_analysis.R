
####
#### MetaboFIMS Paper Revision Annotation analysis of case study 1
####



library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(tidyr)

# =============================================================================
# Path & Configuration Setup
# =============================================================================
if (.Platform$OS.type == "unix") {
  parent_ms_path <- "/Volumes/T7/Arbeit/FIA/MSKohorte_NEU"
} else {
  parent_ms_path <- "D:/Arbeit/FIA/MSKohorte_NEU"
}

ppm_settings <- c("10ppm annotation", "20ppm annotation")

# Metadata column definitions for extracting labels
cv_meta_cols    <- c("SampleID", "class", "index", "batch", "Sex", "Group", "Visit", "Age")
nf_meta_labels  <- c("SampleID", "class", "index", "batch")

# Mass tolerance for matching feature labels back to row-1 m/z values
mz_match_tol    <- 0.001

# =============================================================================
# Helper & Parsing Functions
# =============================================================================

# Count annotations in string columns based on delimiter pattern
count_annotations_regex <- function(column, pattern) {
  ifelse(
    is.na(column) | column == "" | column == "NA", 
    0, 
    str_count(column, pattern) + 1
  )
}

# Generate distribution string (e.g. "2 annos: 15; 3 annos: 4")
get_distribution_string <- function(counts_vec) {
  multi_counts <- counts_vec[counts_vec > 1]
  if (length(multi_counts) == 0) return("None")
  
  dist_table <- table(multi_counts)
  paste(sapply(names(dist_table), function(x) paste0(x, " annos: ", dist_table[x])), collapse = "; ")
}

# Parse mass from "Name (mass;Adduct)" label in CV files
parse_formula_entry <- function(entry) {
  m <- str_match(entry, "^(.*?)\\s*\\(\\s*([-0-9.]+)\\s*;\\s*(.+?)\\s*\\)$")
  if (is.na(m[1, 1])) return(list(name = NA_character_, mass = NA_real_, adduct = NA_character_))
  list(name = str_trim(m[1, 2]), mass = as.numeric(m[1, 3]), adduct = str_trim(m[1, 4]))
}

parse_feature_label <- function(x) {
  parse_formula_entry(str_remove(x, "\\.\\.\\.\\d+$"))
}

# Load and transpose a raw matrix, then calculate all formula/mz annotation counts
process_matrix_file <- function(file_path) {
  if (!file.exists(file_path)) return(NULL)
  
  df_raw <- read_csv(file_path, col_names = FALSE, show_col_types = FALSE)
  
  # Extract true m/z values from row 1 for matching later
  peak_cols <- setdiff(colnames(df_raw), "X1")
  mz_lookup <- tibble(
    peak = peak_cols,
    mz_value = as.numeric(unlist(df_raw[1, peak_cols]))
  )
  
  # Transpose and merge
  df_transposed <- df_raw %>%
    pivot_longer(cols = -X1, names_to = "peak", values_to = "value") %>%
    pivot_wider(names_from = X1, values_from = value) %>%
    left_join(mz_lookup, by = "peak")
  
  # Calculate annotation depths and sources
  df_transposed %>%
    mutate(
      Formula_Annotation_Count = count_annotations_regex(Formula_Matched_NAME, ", "),
      MZ_Annotation_Count      = count_annotations_regex(mz_matched_NAME, "\\), "),
      Final_Annotation_Count   = ifelse(MZ_Annotation_Count == 0, Formula_Annotation_Count, MZ_Annotation_Count)
    )
}

# Map CV Blank feature headers to the base annotation df based on mass
get_matched_subset <- function(file_path, df_base) {
  if (!file.exists(file_path) || is.null(df_base)) return(NULL)
  
  headers <- colnames(read_csv(file_path, n_max = 0, show_col_types = FALSE))
  feature_labels <- setdiff(headers, cv_meta_cols)
  
  if (length(feature_labels) == 0) return(df_base[0, ])
  
  masses <- vapply(feature_labels, function(lbl) parse_feature_label(lbl)$mass, double(1))
  
  # For each mass, map to the matching peak row
  matched_rows <- map_dfr(masses, function(m) {
    if (is.na(m)) return(NULL)
    df_base %>% filter(abs(mz_value - m) < mz_match_tol) %>% head(1)
  })
  
  return(matched_rows)
}

# Calculate comprehensive metrics for a given dataframe
calc_metrics <- function(df, prefix) {
  if (is.null(df) || nrow(df) == 0) {
    res <- tibble(Total_Features = 0, Features_with_0_Annos = 0, Features_with_1_Anno = 0, 
                  Features_with_Multi_Annos = 0, Exact_Multi_Distribution = "None",
                  Annotated_by_BOTH = 0, Annotated_by_mz_ONLY = 0, Annotated_by_FORMULA_ONLY = 0)
  } else {
    res <- tibble(
      Total_Features            = nrow(df),
      Features_with_0_Annos     = sum(df$Final_Annotation_Count == 0),
      Features_with_1_Anno      = sum(df$Final_Annotation_Count == 1),
      Features_with_Multi_Annos = sum(df$Final_Annotation_Count > 1),
      Exact_Multi_Distribution  = get_distribution_string(df$Final_Annotation_Count),
      
      Annotated_by_BOTH         = sum(df$MZ_Annotation_Count > 0 & df$Formula_Annotation_Count > 0),
      Annotated_by_mz_ONLY      = sum(df$MZ_Annotation_Count > 0 & df$Formula_Annotation_Count == 0),
      Annotated_by_FORMULA_ONLY = sum(df$Formula_Annotation_Count > 0 & df$MZ_Annotation_Count == 0)
    )
  }
  
  colnames(res) <- paste0(prefix, "_", colnames(res))
  return(res)
}

# =============================================================================
# Main Analysis Loop
# =============================================================================
master_summary_data <- map_dfr(ppm_settings, function(ppm) {
  
  base_ms_path <- file.path(parent_ms_path, ppm)
  
  # Define file paths dynamically
  paths_NON_Filtered <- list(
    Positive = file.path(base_ms_path, "CSV/feature_matrix_with_formula_and_mz_annotations.csv"),
    Negative = file.path(base_ms_path, "CSV_NEG/feature_matrix_with_formula_and_mz_annotations_neg.csv")
  )
  paths_ANNOTATION_Filtered <- list(
    Positive = file.path(base_ms_path, "CSV/feature_matrix_with_formula_and_mz_annotations_filtered.csv"),
    Negative = file.path(base_ms_path, "CSV_NEG/feature_matrix_with_formula_and_mz_annotations_filtered_neg.csv")
  )
  paths_CV_BLANK_Filtered <- list(
    Positive = file.path(base_ms_path, "CSV/feature_matrix_transposed_MSKohorte_with_metadata_with_age_CV_Blank_filtered.csv"),
    Negative = file.path(base_ms_path, "CSV_NEG/feature_matrix_transposed_MSKohorte_with_metadata_neg_with_age_CV_Blank_filtered.csv")
  )
  
  # Loop through Positive/Negative modes
  mode_summary_data <- imap_dfr(paths_ANNOTATION_Filtered, function(af_path, mode_name) {
    if (!file.exists(af_path)) {
      warning(sprintf("Missing base annotation file for %s [%s]: %s", ppm, mode_name, af_path))
      return(NULL)
    }
    
    # -------------------------------------------------------------------------
    # 1. Process ANNOTATION_Filtered (Base Tracker)
    # -------------------------------------------------------------------------
    df_af <- process_matrix_file(af_path)
    metrics_af <- calc_metrics(df_af, "ANNOTATION_Filtered")
    
    # Export unannotated true errors
    unannotated_features <- df_af %>% filter(Final_Annotation_Count == 0)
    if (nrow(unannotated_features) > 0) {
      unannotated_csv_path <- file.path(base_ms_path, paste0("ERRORS_true_unannotated_", mode_name, ".csv"))
      write_csv(unannotated_features, unannotated_csv_path)
    }
    
    # -------------------------------------------------------------------------
    # 2. Process NON-Filtered File
    # -------------------------------------------------------------------------
    df_nf <- process_matrix_file(paths_NON_Filtered[[mode_name]])
    metrics_nf <- calc_metrics(df_nf, "NON-Filtered")
    
    # -------------------------------------------------------------------------
    # 3. Process CV_BLANK_Filtered via Mass Mapping
    # -------------------------------------------------------------------------
    df_cv <- get_matched_subset(paths_CV_BLANK_Filtered[[mode_name]], df_af)
    metrics_cv <- calc_metrics(df_cv, "CV_BLANK_Filtered")
    
    # -------------------------------------------------------------------------
    # Console Output Summary
    # -------------------------------------------------------------------------
    cat(sprintf("\n=======================================================\n"))
    cat(sprintf("=== [%s] -> %s Mode Detailed Master Summary ===\n", ppm, mode_name))
    cat(sprintf("=======================================================\n"))
    cat(sprintf("NON-Filtered Tier:\n"))
    cat(sprintf("  -> Total Features: %d | Annotated: %d | Unannotated: %d\n", metrics_nf[[1]], metrics_nf[[1]] - metrics_nf[[2]], metrics_nf[[2]]))
    cat(sprintf("  -> Distribution: Single: %d | Multi: %d (%s)\n", metrics_nf[[3]], metrics_nf[[4]], metrics_nf[[5]]))
    
    cat(sprintf("\nANNOTATION_Filtered Tier:\n"))
    cat(sprintf("  -> Total Features: %d | Annotated: %d | Unannotated: %d\n", metrics_af[[1]], metrics_af[[1]] - metrics_af[[2]], metrics_af[[2]]))
    cat(sprintf("  -> BOTH: %d | mz_ONLY: %d | FORMULA_ONLY: %d\n", metrics_af[[6]], metrics_af[[7]], metrics_af[[8]]))
    if (nrow(unannotated_features) > 0) cat(sprintf("  -> ⚠️ ERROR FILE SAVED: %d true unannotated features found!\n", nrow(unannotated_features)))
    
    cat(sprintf("\nCV_BLANK_Filtered Tier (Mapped):\n"))
    cat(sprintf("  -> Total Features: %d | Annotated: %d | Unannotated: %d\n", metrics_cv[[1]], metrics_cv[[1]] - metrics_cv[[2]], metrics_cv[[2]]))
    cat(sprintf("  -> BOTH: %d | mz_ONLY: %d | FORMULA_ONLY: %d\n", metrics_cv[[6]], metrics_cv[[7]], metrics_cv[[8]]))
    cat(sprintf("=======================================================\n"))
    
    # Combine everything into one ultra-wide row
    bind_cols(
      tibble(PPM_Setting = ppm, Experimental_Mode = mode_name),
      metrics_nf,
      metrics_af,
      metrics_cv
    )
  })
  
  return(mode_summary_data)
})

# =============================================================================
# Save Combined Master Output
# =============================================================================
master_output_csv <- file.path(parent_ms_path, "annotation_analysis_comprehensive_summary.csv")
write_csv(master_summary_data, master_output_csv)

cat(sprintf("\n[SUCCESS] Unified analysis completed seamlessly!\n"))
cat(sprintf("Master Summary File saved to: %s\n", master_output_csv))









