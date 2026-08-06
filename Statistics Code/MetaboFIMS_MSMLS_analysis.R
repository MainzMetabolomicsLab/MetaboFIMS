### =========================================================================
### FIA-MS MS1 DETECTION ANALYSIS - WITH BLANK FILTERING
### =========================================================================
### Uses RT-window accumulation, smoothing, and GLOBAL centroiding
### All files share a single centroid list (like document 6)
### Applies 3x blank average filter to reduce false positives
### =========================================================================
### COMPARISON VERSION: Use with independent centroiding version to compare
### =========================================================================

library(MSnbase)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(cowplot)
library(gridExtra)
library(grid)
library(ggpubr)
library(Spectra)
library(S4Vectors)
library(mzR)
library(data.table)
library(signal)  # For Savitzky-Golay smoothing

# --- 1. CONFIGURATION ---
if (Sys.info()['sysname'] == "Darwin") {
  data_dir_MSMLS <- "/Volumes/T7/Arbeit/FIA/MSMLS/MSMLS P1-5/mzML/MSMLS.csv"
  data_dir_mzML  <- "/Volumes/T7/Arbeit/FIA/MSMLS/MSMLS P1-5/mzML"  
  output_dir     <- "/Volumes/T7/Arbeit/FIA/MSMLS/MSMLS P1-5/mzML"
} else {
  data_dir_MSMLS <- "D:/Arbeit/FIA/MSMLS/MSMLS P1-5/mzML/MSMLS.csv"
  data_dir_mzML  <- "D:/Arbeit/FIA/MSMLS/MSMLS P1-5/mzML"
  output_dir     <- "D:/Arbeit/FIA/MSMLS/MSMLS P1-5/mzML"
}

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Parameters
ppm_tolerance_centroiding <- 20  # PPM tolerance for centroiding and peak grouping
ppm_tolerance <- 10              # PPM tolerance for m/z matching to MSMLS targets
intensity_threshold <- 5000      # Minimum intensity to consider a peak detected
blank_multiplier <- 3            # Require sample to be 3x blank average
rt_min <- 5                      # Start of RT window (seconds) for FIA integration
rt_max <- 25                     # End of RT window (seconds) for FIA integration
EXPORT_SPECTRA <- T           # Export centroided spectra to mzML files
EXPORT_SMOOTHED <- T          # Export smoothed accumulated profile spectra

# --- FIA ACCUMULATION / CENTROIDING PARAMETERS ---
common_mz_min <- 70
common_mz_max <- 1000
common_mz_step <- 0.0005

apply_smoothing <- TRUE
smoothing_n <- 7
smoothing_p <- 2

# --- CENTROIDING SETTINGS ---
local_max_window_size <- 2   # Window for finding local maxima

# Create export directories
centroided_export_dir <- file.path(output_dir, "centroided_mzML_GLOBAL")
smoothed_export_dir <- file.path(output_dir, "smoothed_accumulated_mzML_GLOBAL")
dir.create(centroided_export_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(smoothed_export_dir, showWarnings = FALSE, recursive = TRUE)

# Figure parameters
max_figure_width <- 30
top_n_compounds_fig2 <- 30
top_n_compounds_suppl <- 100

# --- PUBLICATION-QUALITY THEME ---
theme_publication <- function(base_size = 11, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5, margin = margin(0,0,10,0)),
      plot.subtitle = element_text(size = rel(0.9), hjust = 0.5, margin = margin(0,0,10,0)),
      plot.caption = element_text(size = rel(0.8), hjust = 1, margin = margin(10,0,0,0)),
      axis.title = element_text(face = "bold", size = rel(1.0)),
      axis.title.x = element_text(margin = margin(10,0,0,0)),
      axis.title.y = element_text(margin = margin(0,10,0,0)),
      axis.text = element_text(size = rel(0.9), color = "black"),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      legend.title = element_text(face = "bold", size = rel(1.0)),
      legend.text = element_text(size = rel(0.9)),
      legend.key.size = unit(1, "lines"),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = NA),
      legend.margin = margin(0,0,0,0),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "black", size = 1),
      panel.grid.major = element_line(color = "grey90", size = 0.25),
      panel.grid.minor = element_line(color = "grey95", size = 0.15),
      strip.background = element_rect(fill = "grey90", color = "black", size = 0.5),
      strip.text = element_text(face = "bold", size = rel(1.0)),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ACS-appropriate color palettes
acs_colors <- list(
  primary = c("#0077BB", "#CC3311", "#009988", "#EE7733", "#33BBEE", "#EE3377"),
  polarity = c("Positive" = "#0077BB", "Negative" = "#CC3311"),
  sequential = c("#FFFFCC", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#0C2C84"),
  diverging = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
)

# --- 2. LOAD LIBRARY ---
cat("Loading MSMLS library...\n")
msmls_lib <- read.csv(data_dir_MSMLS, stringsAsFactors = FALSE) %>%
  dplyr::mutate(
    Position = str_trim(as.character(Position)),
    SMILES = if("SMILES" %in% names(.)) as.character(SMILES) else NA_character_
  )

cat(paste("Library loaded:", nrow(msmls_lib), "entries\n"))

# --- 3. PREPARE TARGET LIST ---
prepare_target_list <- function(library_df) {
  cat("\nPreparing target list from library...\n")
  
  # Positive mode adducts
  pos_targets <- library_df %>%
    dplyr::select(Position, MSMLS_NAME, MSMLS_H, MSMLS_Na, MSMLS_K) %>%
    tidyr::pivot_longer(cols = c(MSMLS_H, MSMLS_Na, MSMLS_K),
                        names_to = "adduct_type",
                        values_to = "target_mz") %>%
    dplyr::filter(!is.na(target_mz)) %>%
    dplyr::mutate(
      polarity = "Positive",
      adduct = case_when(
        adduct_type == "MSMLS_H" ~ "[M+H]+",
        adduct_type == "MSMLS_Na" ~ "[M+Na]+",
        adduct_type == "MSMLS_K" ~ "[M+K]+",
        TRUE ~ adduct_type
      )
    )
  
  # Negative mode adducts
  neg_targets <- library_df %>%
    dplyr::select(Position, MSMLS_NAME, MSMLS_Hneg, MSMLS_Cl, MSMLS_FA) %>%
    tidyr::pivot_longer(cols = c(MSMLS_Hneg, MSMLS_Cl, MSMLS_FA),
                        names_to = "adduct_type",
                        values_to = "target_mz") %>%
    dplyr::filter(!is.na(target_mz)) %>%
    dplyr::mutate(
      polarity = "Negative",
      adduct = case_when(
        adduct_type == "MSMLS_Hneg" ~ "[M-H]-",
        adduct_type == "MSMLS_Cl" ~ "[M+Cl]-",
        adduct_type == "MSMLS_FA" ~ "[M+FA-H]-",
        TRUE ~ adduct_type
      )
    )
  
  # Combine
  all_targets <- dplyr::bind_rows(pos_targets, neg_targets) %>%
    dplyr::select(Position, MSMLS_NAME, target_mz, adduct, polarity) %>%
    dplyr::arrange(Position, MSMLS_NAME, adduct)
  
  cat(paste("Total targets prepared:", nrow(all_targets), "\n"))
  cat(paste("  - Positive mode:", sum(all_targets$polarity == "Positive"), "\n"))
  cat(paste("  - Negative mode:", sum(all_targets$polarity == "Negative"), "\n"))
  cat(paste("Unique compounds:", n_distinct(all_targets$MSMLS_NAME), "\n"))
  
  return(all_targets)
}

target_list <- prepare_target_list(msmls_lib)


# ============================================================================
# FIA RT-WINDOW ACCUMULATION + GLOBAL CENTROIDING FUNCTIONS
# (Following document 6 approach - global centroid list)
# ============================================================================

# --- RT-WINDOW ACCUMULATION ---
accumulate_rt_window <- function(file, mz_grid, rt_min, rt_max) {
  raw_data <- readMSData(files = file, mode = "onDisk")
  ms1_data <- filterMsLevel(raw_data, msLevel = 1)
  spectra_list <- spectra(ms1_data)
  retention_times <- rtime(ms1_data)
  
  valid_indices <- which(retention_times >= rt_min & retention_times <= rt_max)
  selected_spectra <- spectra_list[valid_indices]
  
  if (length(selected_spectra) == 0) {
    return(rep(0, length(mz_grid)))
  }
  
  resampled_intensities <- Reduce("+", lapply(selected_spectra, function(spectrum) {
    approx(
      x = mz(spectrum),
      y = intensity(spectrum),
      xout = mz_grid,
      method = "linear",
      yleft = 0,
      yright = 0
    )$y
  }))
  
  return(resampled_intensities)
}

# --- SMOOTHING ---
smooth_intensity <- function(intensity, n = smoothing_n, p = smoothing_p) {
  signal::sgolayfilt(intensity, p = p, n = n)
}

# --- PPM TOLERANCE ---
get_mz_tolerance <- function(mz, ppm = ppm_tolerance_centroiding) {
  (mz * ppm) / 1e6
}

# --- FIND LOCAL MAXIMA ---
find_local_maxima <- function(mz, intensity, window_size = 5) {
  n <- length(intensity)
  local_maxima <- rep(FALSE, n)
  for(i in 2:(n - 1)) {
    if(intensity[i] > intensity[i - 1] && intensity[i] > intensity[i + 1]) {
      window_start <- max(1, i - window_size)
      window_end <- min(n, i + window_size)
      if(intensity[i] == max(intensity[window_start:window_end])) {
        local_maxima[i] <- TRUE
      }
    }
  }
  idx <- which(local_maxima)
  data.frame(mz = mz[idx], intensity = intensity[idx])
}


# --- 4. ACCUMULATE AND SMOOTH ALL FILES (FIRST PASS) ---
accumulate_all_files <- function(file_paths, rt_window_min, rt_window_max) {
  cat("\n=== ACCUMULATING AND SMOOTHING ALL FILES ===\n")
  
  mz_grid <- seq(common_mz_min, common_mz_max, by = common_mz_step)
  
  all_data <- purrr::map(file_paths, function(file_path) {
    cat(paste("Accumulating:", basename(file_path), "\n"))
    
    is_blank <- grepl("Blank", basename(file_path), ignore.case = TRUE)
    
    pool_id <- if (is_blank) {
      "Blank"
    } else {
      pid <- str_trim(str_extract(basename(file_path), "P[0-9]+[A-Z]"))
      if (is.na(pid)) pid <- "Unknown"
      pid
    }
    
    ms_data <- readMSData(file_path, mode = "onDisk")
    hdr <- fData(ms_data)
    pol <- unique(hdr$polarity[hdr$msLevel == 1])[1]
    polarity_name <- ifelse(pol == 1, "Positive", "Negative")
    
    # Accumulate
    accumulated_intensity <- accumulate_rt_window(
      file_path,
      mz_grid,
      rt_window_min,
      rt_window_max
    )
    
    if (sum(accumulated_intensity) == 0) {
      cat("  Warning: Empty spectrum\n")
      return(NULL)
    }
    
    # Smooth
    if (apply_smoothing) {
      accumulated_intensity <- smooth_intensity(accumulated_intensity)
    }
    
    # Store smoothed profile
    accumulated_smoothed <- data.frame(
      mz = mz_grid,
      intensity = accumulated_intensity
    ) %>%
      dplyr::filter(intensity > 0)
    
    cat(paste("  Accumulated peaks:", nrow(accumulated_smoothed), "\n"))
    
    return(list(
      pool_id = pool_id,
      polarity = polarity_name,
      mz = mz_grid,
      intensity = accumulated_intensity,
      spectrum_smoothed = accumulated_smoothed,
      n_scans = sum(hdr$msLevel == 1 &
                      hdr$retentionTime >= rt_window_min &
                      hdr$retentionTime <= rt_window_max),
      source_file = basename(file_path),
      is_blank = is_blank
    ))
  })
  
  # Remove NULL entries
  all_data <- all_data[!sapply(all_data, is.null)]
  
  return(all_data)
}


# --- 5. GLOBAL CENTROIDING (DOCUMENT 6 APPROACH) ---
create_global_centroids <- function(all_accumulated_data, ppm_tol_centroid, int_threshold) {
  cat("\n=== CREATING GLOBAL CENTROID LIST ===\n")
  
  # Combine all peaks above threshold from ALL samples
  all_peaks <- rbindlist(lapply(all_accumulated_data, function(data) {
    df <- data.frame(
      mz = data$mz,
      intensity = data$intensity,
      sample = data$source_file
    ) %>%
      dplyr::filter(intensity > int_threshold)
    return(df)
  }), fill = TRUE)
  
  cat(paste("Total peaks above threshold (all files):", nrow(all_peaks), "\n"))
  
  # Find local maxima across ALL data
  local_maxima <- find_local_maxima(all_peaks$mz, all_peaks$intensity, 
                                    window_size = local_max_window_size)
  
  cat(paste("Local maxima detected:", nrow(local_maxima), "\n"))
  
  if (nrow(local_maxima) == 0) {
    return(data.frame(mz_centroid = numeric()))
  }
  
  # Sort and group by PPM tolerance
  mz_sorted <- sort(local_maxima$mz)
  int_sorted <- local_maxima$intensity[order(local_maxima$mz)]
  
  mz_diffs <- diff(mz_sorted)
  mz_tolerances <- get_mz_tolerance(head(mz_sorted, -1), ppm_tol_centroid)
  group <- c(1, 1 + cumsum(mz_diffs > mz_tolerances))
  
  # Calculate intensity-weighted centroids
  centroids <- data.frame(mz = mz_sorted, intensity = int_sorted, group = group) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      mz_centroid = sum(mz * intensity) / sum(intensity),
      .groups = "drop"
    )
  
  cat(paste("Global centroids created:", nrow(centroids), "\n"))
  
  return(centroids)
}


# --- 6. INTEGRATE EACH FILE USING GLOBAL CENTROIDS ---
integrate_file_with_global_centroids <- function(file_data, global_centroids, ppm_tol_centroid) {
  
  # Create data.table of full spectrum
  spec_dt <- data.table(
    mz = file_data$mz,
    intensity = file_data$intensity
  )
  
  if (nrow(global_centroids) == 0) {
    return(data.frame(mz = numeric(), intensity = numeric()))
  }
  
  # Set up integration windows
  centroids_dt <- as.data.table(global_centroids)
  centroids_dt[, ppm_tol := get_mz_tolerance(mz_centroid, ppm_tol_centroid)]
  centroids_dt[, mz_min := mz_centroid - ppm_tol]
  centroids_dt[, mz_max := mz_centroid + ppm_tol]
  
  spec_dt[, mz_start := mz]
  spec_dt[, mz_end := mz]
  
  setkey(spec_dt, mz_start, mz_end)
  setkey(centroids_dt, mz_min, mz_max)
  
  # Integrate using foverlaps
  matched <- foverlaps(
    spec_dt,
    centroids_dt,
    by.x = c("mz_start", "mz_end"),
    by.y = c("mz_min", "mz_max"),
    nomatch = 0
  )
  
  # Sum intensities for each centroid
  final_centroids <- matched[, .(intensity = sum(intensity)), by = .(mz = mz_centroid)]
  
  return(as.data.frame(final_centroids[order(mz)]))
}


# --- 7. MATCH TARGETS TO DETECTED PEAKS ---
match_targets_to_spectrum <- function(spectrum_data, target_df, ppm_tol, int_threshold) {
  
  spectrum <- spectrum_data$spectrum
  pool_id <- spectrum_data$pool_id
  polarity <- spectrum_data$polarity
  is_blank <- spectrum_data$is_blank
  
  # For blank files, match ALL targets of the correct polarity
  # For sample files, match only targets for this pool
  if (is_blank) {
    pool_targets <- target_df %>%
      dplyr::filter(polarity == !!polarity)
    cat(paste("  Matching", nrow(pool_targets), "targets in", polarity, "BLANK\n"))
  } else {
    pool_targets <- target_df %>%
      dplyr::filter(Position == pool_id, polarity == !!polarity)
    
    if (nrow(pool_targets) == 0) {
      cat(paste("  No targets for pool", pool_id, "in", polarity, "mode\n"))
      return(NULL)
    }
    cat(paste("  Matching", nrow(pool_targets), "targets in", polarity, "mode\n"))
  }
  
  # Match each target to spectrum
  results <- purrr::map_df(seq_len(nrow(pool_targets)), function(i) {
    target <- pool_targets[i, ]
    target_mz <- target$target_mz
    
    # Calculate PPM window (using MATCHING tolerance, not centroiding)
    mz_window <- target_mz * ppm_tol / 1e6
    
    # Find matching peaks
    matches <- spectrum %>%
      dplyr::filter(
        mz >= (target_mz - mz_window),
        mz <= (target_mz + mz_window),
        intensity >= int_threshold
      )
    
    if (nrow(matches) == 0) {
      return(data.frame(
        pool_id = pool_id,
        compound = target$MSMLS_NAME,
        adduct = target$adduct,
        target_mz = target_mz,
        detected = FALSE,
        detected_mz = NA,
        detected_intensity = NA,
        ppm_error = NA,
        n_matches = 0,
        polarity = polarity,
        source_file = spectrum_data$source_file,
        stringsAsFactors = FALSE
      ))
    } else {
      best_match <- matches %>%
        dplyr::arrange(desc(intensity)) %>%
        dplyr::slice(1)
      
      ppm_err <- abs(best_match$mz - target_mz) / target_mz * 1e6
      
      return(data.frame(
        pool_id = pool_id,
        compound = target$MSMLS_NAME,
        adduct = target$adduct,
        target_mz = target_mz,
        detected = TRUE,
        detected_mz = best_match$mz,
        detected_intensity = best_match$intensity,
        ppm_error = ppm_err,
        n_matches = nrow(matches),
        polarity = polarity,
        source_file = spectrum_data$source_file,
        stringsAsFactors = FALSE
      ))
    }
  })
  
  cat(paste("  Detected:", sum(results$detected), "of", nrow(results), "targets\n"))
  
  return(results)
}

# --- 8. MAIN EXECUTION ---
cat("\n=== FIA-MS MS1 Detection Analysis WITH BLANK FILTERING ===\n")
cat("=== GLOBAL CENTROIDING MODE (all files share centroids) ===\n\n")
cat(paste("RT window:", rt_min, "-", rt_max, "seconds\n"))
cat(paste("PPM tolerance (centroiding):", ppm_tolerance_centroiding, "\n"))
cat(paste("PPM tolerance (MSMLS matching):", ppm_tolerance, "\n"))
cat(paste("Intensity threshold:", intensity_threshold, "\n"))
cat(paste("Blank filter:", blank_multiplier, "x blank average\n"))
cat(paste("Smoothing:", ifelse(apply_smoothing, "ENABLED", "DISABLED"), "\n"))
cat(paste("Centroiding: GLOBAL across all files\n\n"))

# Get all mzML files
mzml_files <- list.files(data_dir_mzML, 
                         pattern = "\\.mz(ML|XML)$", 
                         full.names = TRUE,
                         ignore.case = TRUE)

if (length(mzml_files) == 0) {
  stop("No mzML files found in directory: ", data_dir_mzML)
}

cat(paste("Found", length(mzml_files), "mzML files\n"))

# Separate blank and sample files
blank_files <- mzml_files[grepl("Blank", basename(mzml_files), ignore.case = TRUE)]
sample_files <- mzml_files[!grepl("Blank", basename(mzml_files), ignore.case = TRUE)]

cat(paste("  Blank files:", length(blank_files), "\n"))
cat(paste("  Sample files:", length(sample_files), "\n"))

# --- 9. ACCUMULATE ALL FILES ---
all_accumulated <- accumulate_all_files(mzml_files, rt_min, rt_max)

# --- 10. CREATE GLOBAL CENTROIDS ---
global_centroids <- create_global_centroids(
  all_accumulated,
  ppm_tolerance_centroiding,
  intensity_threshold
)

# --- 11. INTEGRATE EACH FILE WITH GLOBAL CENTROIDS ---
cat("\n=== INTEGRATING FILES WITH GLOBAL CENTROIDS ===\n")

all_spectra <- list()
blank_spectra <- list()

for (i in seq_along(all_accumulated)) {
  file_data <- all_accumulated[[i]]
  
  cat(paste("\nIntegrating:", file_data$source_file, "\n"))
  
  # Integrate using global centroids
  centroided <- integrate_file_with_global_centroids(
    file_data,
    global_centroids,
    ppm_tolerance_centroiding
  ) %>%
    dplyr::filter(intensity > 0) %>%
    dplyr::arrange(mz)
  
  cat(paste("  Centroided peaks:", nrow(centroided), "\n"))
  
  result <- list(
    pool_id = file_data$pool_id,
    polarity = file_data$polarity,
    spectrum = centroided,
    spectrum_smoothed = file_data$spectrum_smoothed,
    n_scans = file_data$n_scans,
    source_file = file_data$source_file,
    is_blank = file_data$is_blank
  )
  
  if (file_data$is_blank) {
    blank_spectra[[length(blank_spectra) + 1]] <- result
  } else {
    all_spectra[[length(all_spectra) + 1]] <- result
  }
}

cat(paste("\n=== Successfully processed", length(all_spectra), "sample files ===\n"))
cat(paste("=== Successfully processed", length(blank_spectra), "blank files ===\n"))

# --- 12. EXPORT SMOOTHED ACCUMULATED SPECTRA (PROFILE MODE) ---
if (EXPORT_SMOOTHED) {
  cat("\n=== EXPORTING SMOOTHED ACCUMULATED SPECTRA (PROFILE) ===\n")
  
  all_files_for_export <- c(blank_spectra, all_spectra)
  
  for (spec in all_files_for_export) {
    sample_name <- gsub("\\.mzML$|\\.mzXML$", "", spec$source_file, ignore.case = TRUE)
    
    spectrum_df <- spec$spectrum_smoothed
    mzs <- as.numeric(spectrum_df$mz)
    ints <- as.numeric(spectrum_df$intensity)
    
    if (length(mzs) == 0) next
    
    tic_value <- sum(ints)
    
    spectra_data <- S4Vectors::DataFrame(
      mz = S4Vectors::I(list(mzs)),
      intensity = S4Vectors::I(list(ints)),
      rtime = mean(c(rt_min, rt_max)),
      msLevel = 1L,
      centroided = FALSE,
      fromFile = 1L,
      scanWindowLowerLimit = common_mz_min,
      scanWindowUpperLimit = common_mz_max,
      tic = tic_value,
      scanIndex = 1L,
      polarity = ifelse(spec$polarity == "Positive", 1L, 0L)
    )
    
    tryCatch({
      sps <- Spectra(spectra_data, backend = MsBackendMemory())
      outfile <- file.path(smoothed_export_dir, paste0(sample_name, "_smoothed_profile.mzML"))
      export(sps, MsBackendMzR(), file = outfile)
      cat(paste("  ✔ Exported:", sample_name, "(TIC =", scales::comma(tic_value), ")\n"))
    }, error = function(e) {
      cat(paste("  ✖ Export failed for", sample_name, ":", e$message, "\n"))
    })
  }
  
  cat(paste("\n✅ Smoothed profile spectra exported to:", smoothed_export_dir, "\n"))
} else {
  cat("\n⛔ Skipping smoothed export (EXPORT_SMOOTHED = FALSE)\n")
}

# --- 13. EXPORT CENTROIDED SPECTRA ---
if (EXPORT_SPECTRA) {
  cat("\n=== EXPORTING CENTROIDED SPECTRA ===\n")
  
  all_files_for_export <- c(blank_spectra, all_spectra)
  
  for (spec in all_files_for_export) {
    sample_name <- gsub("\\.mzML$|\\.mzXML$", "", spec$source_file, ignore.case = TRUE)
    
    mzs <- as.numeric(spec$spectrum$mz)
    ints <- as.numeric(spec$spectrum$intensity)
    
    if (length(mzs) == 0) next
    
    spectra_data <- S4Vectors::DataFrame(
      mz = S4Vectors::I(list(mzs)),
      intensity = S4Vectors::I(list(ints)),
      rtime = mean(c(rt_min, rt_max)),
      msLevel = 1L,
      centroided = TRUE,
      fromFile = 1L
    )
    
    tryCatch({
      sps <- Spectra(spectra_data, backend = MsBackendMemory())
      outfile <- file.path(centroided_export_dir, paste0(sample_name, "_centroided.mzML"))
      export(sps, MsBackendMzR(), file = outfile)
      cat(paste("  ✔ Exported:", sample_name, "(", length(mzs), "peaks )\n"))
    }, error = function(e) {
      cat(paste("  ✖ Export failed for", sample_name, ":", e$message, "\n"))
    })
  }
  
  cat(paste("\n✅ Centroided spectra exported to:", centroided_export_dir, "\n"))
} else {
  cat("\n⛔ Skipping centroid export (EXPORT_SPECTRA = FALSE)\n")
}

# --- 14. PROCESS BLANKS ---
cat("\n=== PROCESSING BLANK FILES ===\n")

if (length(blank_spectra) > 0) {
  blank_detections <- purrr::map_df(blank_spectra, function(spec_data) {
    match_targets_to_spectrum(spec_data, target_list, ppm_tolerance, intensity_threshold)
  })
  
  blank_averages <- blank_detections %>%
    dplyr::filter(detected) %>%
    dplyr::group_by(compound, adduct, polarity, target_mz) %>%
    dplyr::summarise(
      blank_avg_intensity = mean(detected_intensity),
      blank_n_detections = n(),
      .groups = "drop"
    )
  
  cat(paste("\n=== BLANK SUMMARY ===\n"))
  cat(paste("Total blank files processed:", length(blank_spectra), "\n"))
  cat(paste("Unique targets detected in blanks:", nrow(blank_averages), "\n"))
  cat(paste("  - Positive mode:", sum(blank_averages$polarity == "Positive"), "\n"))
  cat(paste("  - Negative mode:", sum(blank_averages$polarity == "Negative"), "\n"))
  
} else {
  cat("*** WARNING: No blank files found! Proceeding without blank filtering. ***\n")
  blank_averages <- data.frame(
    compound = character(),
    adduct = character(),
    polarity = character(),
    target_mz = numeric(),
    blank_avg_intensity = numeric(),
    blank_n_detections = integer(),
    stringsAsFactors = FALSE
  )
}

# --- 15. MATCH TARGETS IN SAMPLES ---
cat("\n=== MATCHING TARGETS IN SAMPLES ===\n")

all_detections <- purrr::map_df(all_spectra, function(spec_data) {
  match_targets_to_spectrum(spec_data, target_list, ppm_tolerance, intensity_threshold)
})

# --- 16. APPLY BLANK FILTER ---
cat("\n=== APPLYING BLANK FILTER ===\n")

all_detections <- all_detections %>%
  dplyr::left_join(blank_averages, 
                   by = c("compound", "adduct", "polarity", "target_mz"),
                   suffix = c("", "_blank")) %>%
  dplyr::mutate(
    blank_avg_intensity = ifelse(is.na(blank_avg_intensity), 0, blank_avg_intensity),
    blank_threshold = blank_avg_intensity * blank_multiplier,
    passes_blank_filter = ifelse(detected, 
                                 detected_intensity >= blank_threshold, 
                                 NA),
    detected_final = detected & (detected_intensity >= blank_threshold)
  )

# --- 17. CALCULATE STATISTICS ---
cat("\n=== DETECTION STATISTICS (BEFORE BLANK FILTER) ===\n")

total_targets <- nrow(all_detections)
total_detected_raw <- sum(all_detections$detected)
detection_rate_raw <- total_detected_raw / total_targets * 100

cat(paste("Total targets tested:", total_targets, "\n"))
cat(paste("Total detected (raw):", total_detected_raw, "\n"))
cat(sprintf("Overall detection rate (raw): %.1f%%\n", detection_rate_raw))

cat("\n=== DETECTION STATISTICS (AFTER BLANK FILTER) ===\n")

total_detected <- sum(all_detections$detected_final)
detection_rate <- total_detected / total_targets * 100
n_filtered_by_blank <- total_detected_raw - total_detected

cat(paste("Total detected (after blank filter):", total_detected, "\n"))
cat(sprintf("Overall detection rate (after blank filter): %.1f%%\n", detection_rate))
cat(paste("Filtered out by blank:", n_filtered_by_blank, 
          sprintf("(%.1f%% of raw detections)\n", 
                  ifelse(total_detected_raw > 0, n_filtered_by_blank/total_detected_raw*100, 0))))

# By polarity
polarity_stats <- all_detections %>%
  dplyr::group_by(polarity) %>%
  dplyr::summarise(
    n_targets = n(),
    n_detected_raw = sum(detected),
    n_detected_final = sum(detected_final),
    detection_rate_raw = n_detected_raw / n_targets * 100,
    detection_rate_final = n_detected_final / n_targets * 100,
    .groups = "drop"
  )

cat("\nDetection by polarity:\n")
print(polarity_stats)

# By compound
compound_stats <- all_detections %>%
  dplyr::group_by(compound) %>%
  dplyr::summarise(
    n_adducts_tested = n(),
    n_adducts_detected_raw = sum(detected),
    n_adducts_detected_final = sum(detected_final),
    detection_rate = n_adducts_detected_final / n_adducts_tested * 100,
    detected_in_files = paste(unique(source_file[detected_final]), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(detection_rate), compound)

cat("\nTop 20 most detected compounds (after blank filter):\n")
print(head(compound_stats, 20))

# By pool
pool_stats <- all_detections %>%
  dplyr::group_by(pool_id, source_file) %>%
  dplyr::summarise(
    n_targets = n(),
    n_detected_raw = sum(detected),
    n_detected_final = sum(detected_final),
    detection_rate_final = n_detected_final / n_targets * 100,
    .groups = "drop"
  ) %>%
  dplyr::arrange(pool_id)

cat("\nDetection by pool (after blank filter):\n")
print(pool_stats)

# By adduct
adduct_stats <- all_detections %>%
  dplyr::group_by(adduct, polarity) %>%
  dplyr::summarise(
    n_targets = n(),
    n_detected_raw = sum(detected),
    n_detected_final = sum(detected_final),
    detection_rate_raw = n_detected_raw / n_targets * 100,
    detection_rate_final = n_detected_final / n_targets * 100,
    mean_intensity = mean(detected_intensity[detected_final], na.rm = TRUE),
    median_intensity = median(detected_intensity[detected_final], na.rm = TRUE),
    sd_intensity = sd(detected_intensity[detected_final], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(polarity, desc(detection_rate_final))

cat("\nDetection by adduct type (after blank filter):\n")
print(adduct_stats)

# --- 18. METABOLITE-LEVEL DETECTION ANALYSIS ---
cat("\n=== METABOLITE-LEVEL DETECTION (AFTER BLANK FILTER) ===\n")

compound_detection <- all_detections %>%
  dplyr::group_by(compound, polarity) %>%
  dplyr::summarise(
    detected_in_polarity = any(detected_final),
    n_adducts_tested = n(),
    n_adducts_detected = sum(detected_final),
    detected_adducts = paste(adduct[detected_final], collapse = ", "),
    .groups = "drop"
  )

pos_detected <- compound_detection %>%
  dplyr::filter(polarity == "Positive", detected_in_polarity == TRUE) %>%
  dplyr::pull(compound) %>%
  unique()

neg_detected <- compound_detection %>%
  dplyr::filter(polarity == "Negative", detected_in_polarity == TRUE) %>%
  dplyr::pull(compound) %>%
  unique()

combined_detected <- unique(c(pos_detected, neg_detected))
total_compounds <- n_distinct(all_detections$compound)

pos_detection_rate <- length(pos_detected) / total_compounds * 100
neg_detection_rate <- length(neg_detected) / total_compounds * 100
combined_detection_rate <- length(combined_detected) / total_compounds * 100

cat(paste("\nTotal unique metabolites tested:", total_compounds, "\n"))
cat(paste("\n--- POSITIVE MODE ---\n"))
cat(paste("Metabolites detected:", length(pos_detected), 
          sprintf("(%.1f%%)\n", pos_detection_rate)))

cat(paste("\n--- NEGATIVE MODE ---\n"))
cat(paste("Metabolites detected:", length(neg_detected), 
          sprintf("(%.1f%%)\n", neg_detection_rate)))

cat(paste("\n--- COMBINED (Both Modes) ---\n"))
cat(paste("Metabolites detected:", length(combined_detected), 
          sprintf("(%.1f%%)\n", combined_detection_rate)))

only_pos <- setdiff(pos_detected, neg_detected)
only_neg <- setdiff(neg_detected, pos_detected)
both_modes <- intersect(pos_detected, neg_detected)

cat(paste("\n--- DETECTION BREAKDOWN ---\n"))
cat(paste("Detected only in positive mode:", length(only_pos), "\n"))
cat(paste("Detected only in negative mode:", length(only_neg), "\n"))
cat(paste("Detected in BOTH modes:", length(both_modes), "\n"))

metabolite_summary <- data.frame(
  Category = c("Total metabolites tested",
               "Detected in positive mode",
               "Detected in negative mode", 
               "Detected in combined modes",
               "Detected only in positive",
               "Detected only in negative",
               "Detected in both modes",
               "NOT detected in either mode"),
  Count = c(total_compounds,
            length(pos_detected),
            length(neg_detected),
            length(combined_detected),
            length(only_pos),
            length(only_neg),
            length(both_modes),
            total_compounds - length(combined_detected)),
  Percentage = c(100,
                 pos_detection_rate,
                 neg_detection_rate,
                 combined_detection_rate,
                 length(only_pos) / total_compounds * 100,
                 length(only_neg) / total_compounds * 100,
                 length(both_modes) / total_compounds * 100,
                 (total_compounds - length(combined_detected)) / total_compounds * 100)
)

cat("\n--- SUMMARY TABLE ---\n")
print(metabolite_summary, row.names = FALSE)

compound_detection_wide <- compound_detection %>%
  dplyr::select(compound, polarity, detected_in_polarity, detected_adducts) %>%
  tidyr::pivot_wider(
    names_from = polarity,
    values_from = c(detected_in_polarity, detected_adducts),
    values_fill = list(detected_in_polarity = FALSE, detected_adducts = "")
  ) %>%
  dplyr::mutate(
    detected_in_any_mode = detected_in_polarity_Positive | detected_in_polarity_Negative
  ) %>%
  dplyr::arrange(desc(detected_in_any_mode), compound)

# --- 19. WRITE OUTPUTS ---
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

get_pol_stat <- function(df, pol, col) {
  val <- df[[col]][df$polarity == pol]
  if (length(val) == 0) return(0) else return(val)
}

summary_stats <- data.frame(
  Parameter = c("Total files processed", "Blank files", "Total targets tested", 
                "Detected (before blank filter)", "Detected (after blank filter)",
                "Filtered by blank", "Overall detection rate (%)", 
                "Positive mode targets", "Positive mode detected",
                "Negative mode targets", "Negative mode detected",
                "Unique compounds tested", "PPM tolerance (centroiding)", "PPM tolerance (matching)", 
                "Intensity threshold", "Blank filter multiplier", "RT window (s)", 
                "Smoothing enabled", "Centroiding mode", "Global centroids count"),
  Value = c(
    length(all_spectra),
    length(blank_files),
    total_targets, 
    total_detected_raw,
    total_detected,
    n_filtered_by_blank,
    round(detection_rate, 1),
    get_pol_stat(polarity_stats, "Positive", "n_targets"),
    get_pol_stat(polarity_stats, "Positive", "n_detected_final"),
    get_pol_stat(polarity_stats, "Negative", "n_targets"),
    get_pol_stat(polarity_stats, "Negative", "n_detected_final"),
    n_distinct(all_detections$compound), 
    ppm_tolerance_centroiding,
    ppm_tolerance,
    intensity_threshold,
    blank_multiplier,
    paste(rt_min, "-", rt_max),
    ifelse(apply_smoothing, "Yes", "No"),
    "Global (all files share centroids)",
    nrow(global_centroids)
  )
)

detailed_file <- file.path(output_dir, paste0("DetailedDetections_BlankFiltered_GLOBAL_", timestamp, ".csv"))
write.csv(all_detections, detailed_file, row.names = FALSE)
cat(paste("\nDetailed results saved to:", detailed_file, "\n"))

compound_file <- file.path(output_dir, paste0("CompoundSummary_BlankFiltered_GLOBAL_", timestamp, ".csv"))
write.csv(compound_stats, compound_file, row.names = FALSE)
cat(paste("Compound summary saved to:", compound_file, "\n"))

pool_file <- file.path(output_dir, paste0("PoolSummary_BlankFiltered_GLOBAL_", timestamp, ".csv"))
write.csv(pool_stats, pool_file, row.names = FALSE)
cat(paste("Pool summary saved to:", pool_file, "\n"))

summary_file <- file.path(output_dir, paste0("OverallSummary_BlankFiltered_GLOBAL_", timestamp, ".csv"))
write.csv(summary_stats, summary_file, row.names = FALSE)
cat(paste("Overall summary saved to:", summary_file, "\n"))

metabolite_file <- file.path(output_dir, paste0("MetaboliteDetection_BlankFiltered_GLOBAL_", timestamp, ".csv"))
write.csv(metabolite_summary, metabolite_file, row.names = FALSE)
cat(paste("Metabolite detection summary saved to:", metabolite_file, "\n"))

compound_detection_file <- file.path(output_dir, paste0("CompoundDetectionByPolarity_BlankFiltered_GLOBAL_", timestamp, ".csv"))
write.csv(compound_detection_wide, compound_detection_file, row.names = FALSE)
cat(paste("Compound detection by polarity saved to:", compound_detection_file, "\n"))

blank_file <- file.path(output_dir, paste0("BlankAverages_GLOBAL_", timestamp, ".csv"))
write.csv(blank_averages, blank_file, row.names = FALSE)
cat(paste("Blank averages saved to:", blank_file, "\n"))

centroids_file <- file.path(output_dir, paste0("GlobalCentroids_", timestamp, ".csv"))
write.csv(global_centroids, centroids_file, row.names = FALSE)
cat(paste("Global centroids saved to:", centroids_file, "\n"))

# --- 20. PUBLICATION-QUALITY PLOTS ---
cat("\n=== Generating publication-quality plots ===\n")

# Create Figure 1 panels
# Create Figure 1 panels
create_figure1_panels <- function(detections_df, detection_column, subtitle_suffix = "") {
  
  compound_det <- detections_df %>%
    dplyr::group_by(compound, polarity) %>%
    dplyr::summarise(
      detected_in_polarity = any(!!sym(detection_column)),
      .groups = "drop"
    )
  
  pos_det <- compound_det %>%
    dplyr::filter(polarity == "Positive", detected_in_polarity == TRUE) %>%
    dplyr::pull(compound) %>%
    unique()
  
  neg_det <- compound_det %>%
    dplyr::filter(polarity == "Negative", detected_in_polarity == TRUE) %>%
    dplyr::pull(compound) %>%
    unique()
  
  combined_det <- unique(c(pos_det, neg_det))
  total_comp <- n_distinct(detections_df$compound)
  
  only_pos <- setdiff(pos_det, neg_det)
  only_neg <- setdiff(neg_det, pos_det)
  both_modes <- intersect(pos_det, neg_det)
  
  pos_rate <- length(pos_det) / total_comp * 100
  neg_rate <- length(neg_det) / total_comp * 100
  combined_rate <- length(combined_det) / total_comp * 100
  
  # Panel A: Overall detection pie chart
  molecule_summary <- data.frame(
    Category = c("Detected", "Not Detected"),
    Count = c(length(combined_det), total_comp - length(combined_det)),
    Percentage = c(combined_rate, 100 - combined_rate)
  )
  
  p1a <- ggplot(molecule_summary, aes(x = "", y = Count, fill = Category)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c("Detected" = acs_colors$primary[1], 
                                 "Not Detected" = "grey80")) +
    geom_text(aes(label = paste0(round(Percentage, 1), "%\n(n=", Count, ")")),
              position = position_stack(vjust = 0.5),
              size = 2, fontface = "bold", color = "white") +
    labs(title = "Overall Metabolite Detection Rate",
         subtitle = sprintf("Total: %d metabolites%s", total_comp, subtitle_suffix)) +
    theme_void(base_size = 6) +
    theme(
      plot.title = element_text(face = "bold", size = 7, hjust = 0.5),
      plot.subtitle = element_text(size = 6, hjust = 0.5, margin = margin(3,0,0,0)),
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      legend.position = "bottom"
    )
  
  # Panel B: Detection by polarity
  polarity_mol_stats <- data.frame(
    polarity = c("Positive", "Negative"),
    n_molecules = c(total_comp, total_comp),
    n_detected = c(length(pos_det), length(neg_det)),
    detection_rate = c(pos_rate, neg_rate)
  )
  
  p1b <- ggplot(polarity_mol_stats, aes(x = polarity, y = detection_rate, fill = polarity)) +
    geom_col(width = 0.7, color = "black", size = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d/%d)", 
                                  detection_rate, n_detected, n_molecules)),
              vjust = -0.5, size = 2, fontface = "bold") +
    scale_fill_manual(values = acs_colors$polarity) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                       expand = expansion(mult = c(0, 0.1))) +
    labs(title = "Metabolite Detection by Polarity",
         subtitle = paste0("Molecule detected (all adducts)", subtitle_suffix),
         x = NULL,
         y = "Detection Rate (%)") +
    theme_publication(base_size = 6) +
    theme(legend.position = "none",
          axis.text.x = element_text(face = "bold", size = 7))
  
  # Panel C: Detection by adduct type
  adduct_stats_plot <- detections_df %>%
    dplyr::group_by(adduct, polarity) %>%
    dplyr::summarise(
      n_targets = n(),
      n_detected = sum(!!sym(detection_column)),
      detection_rate = n_detected / n_targets * 100,
      .groups = "drop"
    )
  
  p1c <- ggplot(adduct_stats_plot, aes(x = reorder(adduct, detection_rate), 
                                       y = detection_rate, fill = polarity)) +
    geom_col(width = 0.7, color = "black", size = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", detection_rate)),
              hjust = -0.2, size = 2, fontface = "bold") +
    scale_fill_manual(values = acs_colors$polarity) +
    scale_y_continuous(limits = c(0, max(adduct_stats_plot$detection_rate) * 1.15),
                       expand = expansion(mult = c(0, 0))) +
    coord_flip() +
    labs(title = "Detection Rate by Adduct Type",
         x = NULL,
         y = "Detection Rate (%)",
         fill = "Polarity") +
    theme_publication(base_size = 6) +
    theme(legend.position = c(0.83, 0.23),
          legend.background = element_rect(fill = "white", color = "black"))
  
  # Panel D: Mode overlap
  venn_data <- data.frame(
    Category = factor(c("Positive\nOnly", "Both\nModes", "Negative\nOnly"),
                      levels = c("Positive\nOnly", "Both\nModes", "Negative\nOnly")),
    Count = c(length(only_pos), length(both_modes), length(only_neg)),
    Percentage = c(length(only_pos)/total_comp*100,
                   length(both_modes)/total_comp*100,
                   length(only_neg)/total_comp*100)
  )
  
  p1d <- ggplot(venn_data, aes(x = Category, y = Count, fill = Category)) +
    geom_col(width = 0.7, color = "black", size = 0.5) +
    geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
              vjust = -0.5, size = 2, fontface = "bold") +
    scale_fill_manual(values = c("Positive\nOnly" = acs_colors$primary[1], 
                                 "Both\nModes" = acs_colors$primary[3], 
                                 "Negative\nOnly" = acs_colors$primary[2])) +
    scale_y_continuous(limits = c(0, max(venn_data$Count) * 1.4),
                       expand = expansion(mult = c(0, 0))) +
    labs(title = "Metabolite Detection by Ion Mode",
         subtitle = sprintf("Total detected: %d of %d (%.1f%%)", 
                            length(combined_det), total_comp, combined_rate),
         x = NULL,
         y = "Number of Metabolites") +
    theme_publication(base_size = 6) +
    theme(legend.position = "none")
  
  return(list(p1a = p1a, p1b = p1b, p1c = p1c, p1d = p1d))
}

# --- CREATE RAW (UNFILTERED) FIGURE 1 ---
cat("\n=== Creating Figure 1 (RAW - Before Blank Filter) ===\n")
raw_panels <- create_figure1_panels(all_detections, "detected", " (before blank filter)")

fig1_raw <- plot_grid(
  raw_panels$p1a, raw_panels$p1b, raw_panels$p1c, raw_panels$p1d,
  ncol = 2, nrow = 2,
  labels = NULL,
  rel_widths = c(1, 1),
  rel_heights = c(1, 1)
)

fig1_raw_file <- file.path(output_dir, paste0("Figure1_Overview_RAW_GLOBAL_", timestamp, ".png"))
ggsave(fig1_raw_file, fig1_raw, width = 12, height = 10, dpi = 600, bg = "white")
cat(paste("Figure 1 (RAW) saved to:", fig1_raw_file, "\n"))

fig1_raw_tiff <- file.path(output_dir, paste0("Figure1_Overview_RAW_GLOBAL_", timestamp, ".tiff"))
ggsave(fig1_raw_tiff, fig1_raw, width = 12, height = 10, dpi = 600, 
       compression = "lzw", bg = "white")

# --- CREATE FILTERED FIGURE 1 ---
cat("\n=== Creating Figure 1 (FILTERED - After Blank Filter) ===\n")
filtered_panels <- create_figure1_panels(all_detections, "detected_final", 
                                         sprintf(" (after %dx blank filter)", blank_multiplier))

fig1_filtered <- plot_grid(
  filtered_panels$p1a, filtered_panels$p1b, filtered_panels$p1c, filtered_panels$p1d,
  ncol = 2, nrow = 2,
  labels = NULL,
  rel_widths = c(1, 1),
  rel_heights = c(1, 1)
)

fig1_filtered_file <- file.path(output_dir, paste0("Figure1_Overview_FILTERED_GLOBAL_", timestamp, ".png"))
ggsave(fig1_filtered_file, fig1_filtered, width = 12, height = 10, dpi = 600, bg = "white")
cat(paste("Figure 1 (FILTERED) saved to:", fig1_filtered_file, "\n"))

fig1_filtered_tiff <- file.path(output_dir, paste0("Figure1_Overview_FILTERED_GLOBAL_", timestamp, ".tiff"))
ggsave(fig1_filtered_tiff, fig1_filtered, width = 12, height = 10, dpi = 600, 
       compression = "lzw", bg = "white")

# Save individual panels (FILTERED version only)
fig1a_file <- file.path(output_dir, paste0("Figure1A_MetaboliteDetection_FILTERED_GLOBAL_", timestamp, ".png"))
ggsave(fig1a_file, filtered_panels$p1a, width = 2.3, height = 2, dpi = 600, bg = "white")

fig1b_file <- file.path(output_dir, paste0("Figure1B_DetectionByPolarity_FILTERED_GLOBAL_", timestamp, ".png"))
ggsave(fig1b_file, filtered_panels$p1b, width = 2.3, height = 2, dpi = 600, bg = "white")

fig1c_file <- file.path(output_dir, paste0("Figure1C_DetectionByAdduct_FILTERED_GLOBAL_", timestamp, ".png"))
ggsave(fig1c_file, filtered_panels$p1c, width = 2.3, height = 2, dpi = 600, bg = "white")

fig1d_file <- file.path(output_dir, paste0("Figure1D_ModeOverlap_FILTERED_GLOBAL_", timestamp, ".png"))
ggsave(fig1d_file, filtered_panels$p1d, width = 2.3, height = 2, dpi = 600, bg = "white")

cat("\nFigure 1 individual panels (filtered) saved\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(paste("Output directory:", output_dir, "\n\n"))

cat("=== OUTPUTS GENERATED ===\n")

cat("\n📊 Figure 1 - TWO VERSIONS:\n")
cat("  1. Figure1_Overview_RAW_GLOBAL - Before blank filtering\n")
cat("  2. Figure1_Overview_FILTERED_GLOBAL - After blank filtering\n")
cat("  Plus individual panels (A, B, C, D) for filtered version\n")

cat("\n📁 Export Directories:\n")
cat(paste("  Smoothed:", smoothed_export_dir, "\n"))
cat(paste("  Centroided:", centroided_export_dir, "\n"))

if (EXPORT_SMOOTHED) {
  cat("\n📁 Smoothed Accumulated Profile Spectra:\n")
  cat(paste("  Total files:", length(blank_spectra) + length(all_spectra), "\n"))
  cat("  Format: Profile mode (Savitzky-Golay smoothed)\n")
}

if (EXPORT_SPECTRA) {
  cat("\n📁 Centroided Spectra:\n")
  cat(paste("  Total files:", length(blank_spectra) + length(all_spectra), "\n"))
  cat("  Format: Centroided (GLOBAL centroids with noise integration)\n")
}

cat("\n📊 Global Centroiding:\n")
cat(paste("  Total global centroids:", nrow(global_centroids), "\n"))
cat("  All files share the same centroid list\n")

cat("\n📈 Detection Statistics:\n")
cat(sprintf("  Before blank filter: %.1f%% (%d detected)\n", 
            detection_rate_raw, total_detected_raw))
cat(sprintf("  After %dx blank filter: %.1f%% (%d detected)\n", 
            blank_multiplier, detection_rate, total_detected))
cat(sprintf("  Filtered by blank: %d detections (%.1f%%)\n",
            n_filtered_by_blank, 
            ifelse(total_detected_raw > 0, n_filtered_by_blank/total_detected_raw*100, 0)))

cat("\n🧪 Metabolite-Level Results:\n")
cat(sprintf("  Total metabolites: %d\n", total_compounds))
cat(sprintf("  Detected (filtered): %d (%.1f%%)\n", 
            length(combined_detected), combined_detection_rate))
cat(sprintf("    - Positive only: %d\n", length(only_pos)))
cat(sprintf("    - Negative only: %d\n", length(only_neg)))
cat(sprintf("    - Both modes: %d\n", length(both_modes)))

cat("\n✅ All outputs saved with timestamp:", timestamp, "\n")
cat("\n📝 Processing mode: GLOBAL centroiding (all files share centroids)\n")
cat("💡 Compare with INDEPENDENT centroiding version to assess differences\n")


# ============================================================================
# DROP-IN ADDON: Combined Detectability for [M+H]+ and [M-H]- only
# Paste at the END of your existing script
# ============================================================================

cat("\n=== ADDON: [M+H]+ / [M-H]- COMBINED DETECTABILITY ===\n")

# Filter to primary adducts only
mh_detections <- all_detections %>%
  dplyr::filter(
    (polarity == "Positive" & adduct == "[M+H]+") |
      (polarity == "Negative" & adduct == "[M-H]-")
  )

# Per-compound: detected in [M+H]+ and/or [M-H]-
mh_compound <- mh_detections %>%
  dplyr::group_by(compound, polarity, adduct) %>%
  dplyr::summarise(
    detected_final = any(detected_final),
    best_intensity = max(detected_intensity[detected_final], na.rm = TRUE),
    .groups = "drop"
  )

mh_pos <- mh_compound %>% dplyr::filter(polarity == "Positive", detected_final) %>% dplyr::pull(compound)
mh_neg <- mh_compound %>% dplyr::filter(polarity == "Negative", detected_final) %>% dplyr::pull(compound)
mh_both <- intersect(mh_pos, mh_neg)
mh_only_pos <- setdiff(mh_pos, mh_neg)
mh_only_neg <- setdiff(mh_neg, mh_pos)
mh_combined <- unique(c(mh_pos, mh_neg))
total_comp <- n_distinct(all_detections$compound)

cat(sprintf("Total compounds tested:          %d\n", total_comp))
cat(sprintf("Detected via [M+H]+ only:        %d (%.1f%%)\n", length(mh_pos),      length(mh_pos)      / total_comp * 100))
cat(sprintf("Detected via [M-H]- only:        %d (%.1f%%)\n", length(mh_neg),      length(mh_neg)      / total_comp * 100))
cat(sprintf("Detected via both [M+H]+/[M-H]-: %d (%.1f%%)\n", length(mh_both),     length(mh_both)     / total_comp * 100))
cat(sprintf("Detected in either mode:         %d (%.1f%%)\n", length(mh_combined), length(mh_combined) / total_comp * 100))
cat(sprintf("Not detected in either mode:     %d (%.1f%%)\n",
            total_comp - length(mh_combined),
            (total_comp - length(mh_combined)) / total_comp * 100))

# Summary table
mh_summary <- data.frame(
  Category = c(
    "Total compounds tested",
    "Detected via [M+H]+ (pos only)",
    "Detected via [M-H]- (neg only)",
    "Detected via both [M+H]+ and [M-H]-",
    "Detected in [M+H]+ OR [M-H]- (combined)",
    "Not detected in either"
  ),
  Count = c(
    total_comp,
    length(mh_only_pos),
    length(mh_only_neg),
    length(mh_both),
    length(mh_combined),
    total_comp - length(mh_combined)
  ),
  Percentage = round(c(
    100,
    length(mh_only_pos)  / total_comp * 100,
    length(mh_only_neg)  / total_comp * 100,
    length(mh_both)      / total_comp * 100,
    length(mh_combined)  / total_comp * 100,
    (total_comp - length(mh_combined)) / total_comp * 100
  ), 1)
)

print(mh_summary, row.names = FALSE)

# --- Figure ---
set2 <- brewer.pal(8, "Set2")

mh_plot_df <- data.frame(
  Category = factor(
    c("[M+H]+\nOnly", "Both", "[M-H]-\nOnly"),
    levels = c("[M+H]+\nOnly", "Both", "[M-H]-\nOnly")
  ),
  Count      = c(length(mh_only_pos), length(mh_both), length(mh_only_neg)),
  Percentage = round(c(
    length(mh_only_pos) / total_comp * 100,
    length(mh_both)     / total_comp * 100,
    length(mh_only_neg) / total_comp * 100
  ), 1)
)

p_mh <- ggplot(mh_plot_df, aes(x = Category, y = Count, fill = Category)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
            vjust = -0.3, size = 2) +
  scale_fill_manual(values = c(
    "[M+H]+\nOnly" = set2[1],
    "Both"         = set2[3],
    "[M-H]-\nOnly" = set2[2]
  )) +
  scale_y_continuous(
    limits = c(0, max(mh_plot_df$Count) * 1.28),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title    = "Detectability via [M+H]+ and [M-H]-",
    subtitle = sprintf("Combined: %d/%d (%.1f%%)", length(mh_combined), total_comp,
                       length(mh_combined) / total_comp * 100),
    x = NULL, y = "Number of Compounds"
  ) +
  theme_bw(base_size = 6) +
  theme(legend.position = "none")

ggsave(
  file.path(output_dir, paste0("MH_only_Detectability_", timestamp, ".png")),
  p_mh, width = 2.3, height = 2, dpi = 300, bg = "white"
)

# Save CSV
write.csv(
  mh_summary,
  file.path(output_dir, paste0("MH_only_Detectability_Summary_", timestamp, ".csv")),
  row.names = FALSE
)

cat(sprintf("\n✅ [M+H]+/[M-H]- analysis saved (figure + CSV)\n"))







# ============================================================================
# ADDON: PPM ERROR ANALYSIS OF MSMLS-ANNOTATED FEATURES
# Paste at the END of your existing script
# ============================================================================

cat("\n=== ADDON: PPM ERROR ANALYSIS (MSMLS-ANNOTATED FEATURES) ===\n")

# Restrict to features that passed detection AND blank filtering
ppm_data <- all_detections %>%
  dplyr::filter(detected_final, !is.na(ppm_error))

cat(paste("Annotated features (detected_final) with ppm error:", nrow(ppm_data), "\n"))

if (nrow(ppm_data) == 0) {
  cat("*** WARNING: No detected_final features with ppm_error found. Skipping ppm error analysis. ***\n")
} else {
  
  # --- Overall summary statistics ---
  ppm_overall_stats <- ppm_data %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_ppm = mean(ppm_error),
      median_ppm = median(ppm_error),
      sd_ppm = sd(ppm_error),
      min_ppm = min(ppm_error),
      max_ppm = max(ppm_error),
      q25_ppm = quantile(ppm_error, 0.25),
      q75_ppm = quantile(ppm_error, 0.75),
      q95_ppm = quantile(ppm_error, 0.95),
      pct_within_5ppm = mean(ppm_error <= 5) * 100,
      pct_within_10ppm = mean(ppm_error <= 10) * 100
    )
  
  cat("\n--- OVERALL PPM ERROR STATISTICS ---\n")
  print(ppm_overall_stats)
  
  # --- By polarity ---
  ppm_polarity_stats <- ppm_data %>%
    dplyr::group_by(polarity) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_ppm = mean(ppm_error),
      median_ppm = median(ppm_error),
      sd_ppm = sd(ppm_error),
      q95_ppm = quantile(ppm_error, 0.95),
      pct_within_5ppm = mean(ppm_error <= 5) * 100,
      pct_within_10ppm = mean(ppm_error <= 10) * 100,
      .groups = "drop"
    )
  
  cat("\n--- PPM ERROR BY POLARITY ---\n")
  print(ppm_polarity_stats)
  
  # --- By adduct ---
  ppm_adduct_stats <- ppm_data %>%
    dplyr::group_by(adduct, polarity) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_ppm = mean(ppm_error),
      median_ppm = median(ppm_error),
      sd_ppm = sd(ppm_error),
      q95_ppm = quantile(ppm_error, 0.95),
      .groups = "drop"
    ) %>%
    dplyr::arrange(polarity, desc(n))
  
  cat("\n--- PPM ERROR BY ADDUCT ---\n")
  print(ppm_adduct_stats)
  
  # --- By pool (to check for systematic per-pool mass calibration drift) ---
  ppm_pool_stats <- ppm_data %>%
    dplyr::group_by(pool_id, source_file) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_ppm = mean(ppm_error),
      median_ppm = median(ppm_error),
      sd_ppm = sd(ppm_error),
      .groups = "drop"
    ) %>%
    dplyr::arrange(pool_id)
  
  cat("\n--- PPM ERROR BY POOL ---\n")
  print(ppm_pool_stats)
  
  # --- Write CSVs ---
  ppm_detail_file <- file.path(output_dir, paste0("PPMError_MSMLS_Detailed_", timestamp, ".csv"))
  write.csv(ppm_data %>%
              dplyr::select(pool_id, source_file, compound, adduct, polarity,
                            target_mz, detected_mz, ppm_error, detected_intensity),
            ppm_detail_file, row.names = FALSE)
  cat(paste("\nDetailed ppm error data saved to:", ppm_detail_file, "\n"))
  
  ppm_overall_file <- file.path(output_dir, paste0("PPMError_MSMLS_OverallSummary_", timestamp, ".csv"))
  write.csv(ppm_overall_stats, ppm_overall_file, row.names = FALSE)
  cat(paste("Overall ppm error summary saved to:", ppm_overall_file, "\n"))
  
  ppm_polarity_file <- file.path(output_dir, paste0("PPMError_MSMLS_ByPolarity_", timestamp, ".csv"))
  write.csv(ppm_polarity_stats, ppm_polarity_file, row.names = FALSE)
  
  ppm_adduct_file <- file.path(output_dir, paste0("PPMError_MSMLS_ByAdduct_", timestamp, ".csv"))
  write.csv(ppm_adduct_stats, ppm_adduct_file, row.names = FALSE)
  
  ppm_pool_file <- file.path(output_dir, paste0("PPMError_MSMLS_ByPool_", timestamp, ".csv"))
  write.csv(ppm_pool_stats, ppm_pool_file, row.names = FALSE)
  cat(paste("By-polarity / by-adduct / by-pool ppm error summaries saved\n"))
  
  # --- Histogram: overall ppm error distribution, faceted by polarity ---
  p_ppm_hist <- ggplot(ppm_data, aes(x = ppm_error, fill = polarity)) +
    geom_histogram(binwidth = 0.5, color = "black", size = 0.2, alpha = 0.85,
                   boundary = 0) +
    geom_vline(data = ppm_polarity_stats,
               aes(xintercept = median_ppm),
               linetype = "dashed", color = "black", size = 0.5) +
    facet_wrap(~polarity, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = acs_colors$polarity) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    labs(title = "Mass Accuracy of MSMLS-Annotated Features",
         subtitle = sprintf("n = %d detections (after blank filter) | dashed line = median",
                            nrow(ppm_data)),
         x = "ppm error",
         y = "Count") +
    theme_publication(base_size = 9) +
    theme(legend.position = "none")
  
  ppm_hist_file <- file.path(output_dir, paste0("Figure_PPMError_MSMLS_Histogram_", timestamp, ".png"))
  ggsave(ppm_hist_file, p_ppm_hist, width = 6, height = 6, dpi = 600, bg = "white")
  cat(paste("\nPPM error histogram saved to:", ppm_hist_file, "\n"))
  
  ppm_hist_tiff <- file.path(output_dir, paste0("Figure_PPMError_MSMLS_Histogram_", timestamp, ".tiff"))
  ggsave(ppm_hist_tiff, p_ppm_hist, width = 6, height = 6, dpi = 600,
         compression = "lzw", bg = "white")
  
  # --- Boxplot by adduct (helps spot adduct-specific mass calibration issues) ---
  p_ppm_box <- ggplot(ppm_data, aes(x = reorder(adduct, ppm_error, FUN = median),
                                    y = ppm_error, fill = polarity)) +
    geom_boxplot(outlier.size = 0.6, size = 0.3) +
    scale_fill_manual(values = acs_colors$polarity) +
    coord_flip() +
    labs(title = "PPM Error by Adduct Type",
         subtitle = "MSMLS-annotated features, after blank filter",
         x = NULL, y = "ppm error", fill = "Polarity") +
    theme_publication(base_size = 9)
  
  ppm_box_file <- file.path(output_dir, paste0("Figure_PPMError_MSMLS_ByAdduct_Boxplot_", timestamp, ".png"))
  ggsave(ppm_box_file, p_ppm_box, width = 6, height = 5, dpi = 600, bg = "white")
  cat(paste("PPM error by-adduct boxplot saved to:", ppm_box_file, "\n"))
  
  cat(sprintf("\n✅ PPM error analysis complete: mean = %.2f ppm, median = %.2f ppm, 95th pct = %.2f ppm\n",
              ppm_overall_stats$mean_ppm, ppm_overall_stats$median_ppm, ppm_overall_stats$q95_ppm))
}


# --- Boxplot by adduct, ordered by MEAN (comparison to median-ordered version) ---
p_ppm_box_mean <- ggplot(ppm_data, aes(x = reorder(adduct, ppm_error, FUN = mean),
                                       y = ppm_error, fill = polarity)) +
  geom_boxplot(outlier.size = 0.6, size = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 1.8,
               fill = "white", color = "black", stroke = 0.5) +
  scale_fill_manual(values = acs_colors$polarity) +
  coord_flip() +
  labs(title = "PPM Error by Adduct Type (ordered by mean)",
       subtitle = "MSMLS-annotated features, after blank filter | \u25c7 = mean",
       x = NULL, y = "ppm error", fill = "Polarity") +
  theme_publication(base_size = 9)

ppm_box_mean_file <- file.path(output_dir, paste0("Figure_PPMError_MSMLS_ByAdduct_Boxplot_MeanOrdered_", timestamp, ".png"))
ggsave(ppm_box_mean_file, p_ppm_box_mean, width = 6, height = 5, dpi = 600, bg = "white")
cat(paste("PPM error by-adduct boxplot (mean-ordered) saved to:", ppm_box_mean_file, "\n"))