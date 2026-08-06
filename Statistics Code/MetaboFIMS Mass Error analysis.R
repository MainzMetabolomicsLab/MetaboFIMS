
library(MSnbase)
library(Rdisop)
library(R.utils)
library(webchem)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringdist)
library(BiocParallel)
library(Spectra)
library(S4Vectors)
library(mzR)

get_os <- function() {
  sys <- Sys.info()[["sysname"]]
  if (.Platform$OS.type == "windows") return("windows")
  if (sys == "Darwin") return("macos")
  if (sys == "Linux") return("linux")
  return("unknown")
}
OS_TYPE <- get_os()

default_path <- function(win, mac, linux) {
  if (OS_TYPE == "windows") return(win)
  if (OS_TYPE == "macos")   return(mac)
  if (OS_TYPE == "linux")   return(linux)
  return(mac)
}

file_path <- default_path(
  win   = "D:/Arbeit/FIA/HMDB_BLOOD_ENDO_FOOD_DRUG_DRUGMET.csv",
  mac   = "/Volumes/T7/Arbeit/FIA/HMDB_BLOOD_ENDO_FOOD_DRUG_DRUGMET.csv",
  linux = "~/HMDB_BLOOD.csv"
)

data_dir_raw_files_pos <- default_path(
  win   = "D:/Arbeit/FIA/processing/raw_mzML",
  mac   = "/Volumes/T7/Arbeit/FIA/processing/raw_mzML",
  linux = "~/processing/raw_mzML"
)

data_dir_raw_files_neg <- default_path(
  win   = "D:/Arbeit/FIA/processing_neg/raw_mzML",
  mac   = "/Volumes/T7/Arbeit/FIA/processing_neg/raw_mzML",
  linux = "~/processing_neg/raw_mzML"
)

output_dir <- default_path(
  win   = "D:/Arbeit/FIA",
  mac   = "/Volumes/T7/Arbeit/FIA",
  linux = "~"
)

rt_min <- 5           
rt_max <- 25         

# Positive-mode targets
target_mz_labels_pos <- c(
  "Glucose [M+Na]+" = 203.052606,
  "PC 34:2 [M+H]+"  = 758.569476
)

# Negative-mode targets
target_mz_labels_neg <- c(
  "Glucose [M-H]-"  = 179.056112,
  "PC 34:2 [M+Cl]-" = 792.531602
)

apply_smoothing     <- TRUE
smoothing_n          <- 7      
smoothing_p          <- 2      
ppm_tolerance        <- 20     
intensity_threshold  <- 5000   
target_window_da     <- 0.02   



process_by_rt_window <- function(file, common_mz_grid, rt_min, rt_max) {
  raw_data <- readMSData(files = file, mode = "onDisk")
  ms1_data <- filterMsLevel(raw_data, msLevel = 1)
  spectra_list <- spectra(ms1_data)
  retention_times <- rtime(ms1_data)
  
  valid_indices <- which(retention_times >= rt_min & retention_times <= rt_max)
  selected_spectra <- spectra_list[valid_indices]
  
  if (length(selected_spectra) == 0)
    return(rep(0, length(common_mz_grid)))
  
  Reduce("+", lapply(selected_spectra, function(spectrum) {
    approx(
      x = mz(spectrum),
      y = intensity(spectrum),
      xout = common_mz_grid,
      method = "linear",
      yleft = 0,
      yright = 0
    )$y
  }))
}

smooth_intensity <- function(intensity, n = smoothing_n, p = smoothing_p) {
  if (!requireNamespace("signal", quietly = TRUE)) install.packages("signal")
  return(signal::sgolayfilt(intensity, p = p, n = n))
}

get_mz_tolerance <- function(mz, ppm = ppm_tolerance) {
  (mz * ppm) / 1e6
}

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

targeted_centroid_one_sample <- function(sample_name, spectrum, theoretical_mz,
                                         target_window_da, intensity_threshold,
                                         ppm_tolerance, window_size = 2, centroid_frac = 0.5) {
  local_spectrum <- spectrum %>%
    filter(mz >= theoretical_mz - target_window_da & mz <= theoretical_mz + target_window_da)
  
  empty_result <- data.frame(
    sample = sample_name, theoretical_mz = theoretical_mz,
    mz_centroid = NA_real_, centroid_intensity = 0
  )
  
  if (nrow(local_spectrum) < 3) return(empty_result)
  
  peaks <- local_spectrum %>% filter(intensity > intensity_threshold)
  if (nrow(peaks) < 3) return(empty_result)
  
  local_maxima <- find_local_maxima(peaks$mz, peaks$intensity, window_size = window_size)
  if (nrow(local_maxima) == 0) return(empty_result)
  
  apex_idx <- which.max(local_maxima$intensity)
  apex_mz  <- local_maxima$mz[apex_idx]
  apex_int <- local_maxima$intensity[apex_idx]
  tol      <- get_mz_tolerance(apex_mz, ppm_tolerance)
  
  window <- local_spectrum %>% filter(mz >= apex_mz - tol & mz <= apex_mz + tol)
  keep   <- window$intensity >= (centroid_frac * apex_int)
  if (sum(keep) == 0) return(empty_result)
  
  data.frame(
    sample              = sample_name,
    theoretical_mz      = theoretical_mz,
    mz_centroid         = sum(window$mz[keep] * window$intensity[keep]) / sum(window$intensity[keep]),
    centroid_intensity  = sum(window$intensity[keep])
  )
}


run_polarity_pipeline <- function(polarity_label, data_dir_raw_files, target_mz_labels,
                                  rt_min, rt_max, apply_smoothing,
                                  intensity_threshold, ppm_tolerance, target_window_da) {
  
  message("\n############################################################")
  message("### Processing polarity: ", polarity_label, "  (", data_dir_raw_files, ")")
  message("############################################################")
  
  target_mzs <- unname(target_mz_labels)
  
  # --- Resampling (load from cache if present, else build from mzML/mzXML) ---
  mzxml_files <- list.files(data_dir_raw_files, pattern = "\\.(mzXML|mzML)$", full.names = TRUE)
  common_mz_grid <- seq(from = 70, to = 1000, by = 0.0005)
  resampled_data_file <- file.path(data_dir_raw_files, "resampled_data.rds")
  
  if (file.exists(resampled_data_file)) {
    message("[", polarity_label, "] Resampled data RDS already exists. Loading from file...")
    resampled_data_list <- readRDS(resampled_data_file)
  } else {
    message("[", polarity_label, "] No RDS found. Processing mzML/mzXML files using RT-window...")
    resampled_data_list <- list()
    total_files <- length(mzxml_files)
    
    for (i in seq_along(mzxml_files)) {
      file <- mzxml_files[i]
      message(sprintf("[%s] Processing file %d/%d: %s", polarity_label, i, total_files, basename(file)))
      
      intensities <- process_by_rt_window(file, common_mz_grid, rt_min = rt_min, rt_max = rt_max)
      sample_name <- basename(tools::file_path_sans_ext(file))
      resampled_data_list[[sample_name]] <- data.frame(mz = common_mz_grid, intensity = intensities)
    }
    
    saveRDS(resampled_data_list, file = resampled_data_file)
    message("[", polarity_label, "] Saved RT-window resampled data to RDS.")
  }
  
  if (length(resampled_data_list) == 0) {
    warning("[", polarity_label, "] No samples found in ", data_dir_raw_files, " -- skipping.")
    return(data.table())
  }
  
  if (apply_smoothing) {
    message("[", polarity_label, "] Smoothing ENABLED.")
    for (sample_name in names(resampled_data_list)) {
      spectrum <- resampled_data_list[[sample_name]]
      spectrum$intensity <- smooth_intensity(spectrum$intensity)
      resampled_data_list[[sample_name]] <- spectrum
    }
  } else {
    message("[", polarity_label, "] Smoothing DISABLED.")
  }
  gc()
  
  centroided_peaks_dt <- rbindlist(lapply(names(resampled_data_list), function(sample_name) {
    message("[", polarity_label, "] Targeted centroiding sample: ", sample_name)
    rbindlist(lapply(target_mzs, function(theoretical_mz) {
      targeted_centroid_one_sample(
        sample_name          = sample_name,
        spectrum             = resampled_data_list[[sample_name]],
        theoretical_mz       = theoretical_mz,
        target_window_da     = target_window_da,
        intensity_threshold  = intensity_threshold,
        ppm_tolerance        = ppm_tolerance,
        window_size          = 2
      )
    }), fill = TRUE)
  }), fill = TRUE)
  
  message("[", polarity_label, "] Targeted centroiding complete: ", nrow(centroided_peaks_dt),
          " sample x target rows across ", length(unique(centroided_peaks_dt$sample)), " samples.")
  
  report <- as.data.table(centroided_peaks_dt)
  report[, Abs_Error := mz_centroid - theoretical_mz]
  report[, PPM_Error := (Abs_Error / theoretical_mz) * 1e6]
  
  label_lookup <- data.table(Metabolite = names(target_mz_labels), theoretical_mz = unname(target_mz_labels))
  report <- merge(report, label_lookup, by = "theoretical_mz", all.x = TRUE)
  
  setnames(report,
           old = c("sample", "theoretical_mz", "mz_centroid", "centroid_intensity"),
           new = c("File", "Theoretical_mz", "Observed_mz", "Intensity"))
  
  report[, Polarity := polarity_label]
  
  setcolorder(report,
              c("Polarity", "File", "Metabolite", "Theoretical_mz", "Observed_mz", "Intensity", "Abs_Error", "PPM_Error"))
  
  report
}


mass_error_report_pos <- run_polarity_pipeline(
  polarity_label       = "positive",
  data_dir_raw_files   = data_dir_raw_files_pos,
  target_mz_labels     = target_mz_labels_pos,
  rt_min               = rt_min,
  rt_max               = rt_max,
  apply_smoothing      = apply_smoothing,
  intensity_threshold  = intensity_threshold,
  ppm_tolerance        = ppm_tolerance,
  target_window_da     = target_window_da
)

mass_error_report_neg <- run_polarity_pipeline(
  polarity_label       = "negative",
  data_dir_raw_files   = data_dir_raw_files_neg,
  target_mz_labels     = target_mz_labels_neg,
  rt_min               = rt_min,
  rt_max               = rt_max,
  apply_smoothing      = apply_smoothing,
  intensity_threshold  = intensity_threshold,
  ppm_tolerance        = ppm_tolerance,
  target_window_da     = target_window_da
)

mass_error_report <- rbindlist(list(mass_error_report_pos, mass_error_report_neg), fill = TRUE)
setorder(mass_error_report, Polarity, Metabolite, File)


cat("\n===== Consolidated Mass Error Analysis Report (per file, both polarities) =====\n")
print(as.data.frame(mass_error_report))
cat("==================================================================================\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
fwrite(mass_error_report, file.path(output_dir, "mass_error_analysis_summary.csv"))



mass_error_summary_stats <- mass_error_report[, .(
  N_Files            = .N,
  N_Detected         = sum(!is.na(PPM_Error)),
  Detection_Rate     = round(sum(!is.na(PPM_Error)) / .N, 3),
  Mean_PPM_Error     = round(mean(PPM_Error, na.rm = TRUE), 3),
  Median_PPM_Error   = round(median(PPM_Error, na.rm = TRUE), 3),
  SD_PPM_Error       = round(sd(PPM_Error, na.rm = TRUE), 3),
  Min_PPM_Error      = round(min(PPM_Error, na.rm = TRUE), 3),
  Max_PPM_Error      = round(max(PPM_Error, na.rm = TRUE), 3),
  Mean_Abs_Error     = round(mean(Abs_Error, na.rm = TRUE), 6),
  Mean_Intensity     = round(mean(Intensity, na.rm = TRUE), 1)
), by = .(Polarity, Metabolite, Theoretical_mz)]

setorder(mass_error_summary_stats, Polarity, Metabolite)

cat("\n===== Per-Metabolite Mass Accuracy Summary (both polarities) =====\n")
print(as.data.frame(mass_error_summary_stats))
cat("=====================================================================\n")

fwrite(mass_error_summary_stats, file.path(output_dir, "mass_error_summary_stats.csv"))
message("✔ Exported: mass_error_analysis_summary.csv (per file) and mass_error_summary_stats.csv (per metabolite), both combining positive + negative mode, to: ", output_dir)