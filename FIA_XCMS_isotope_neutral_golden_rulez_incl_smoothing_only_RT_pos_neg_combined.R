
#load all required R packages
{
  #library(xcms)
  library(MSnbase)
  #library(CAMERA)
  library(Rdisop)  # For molecular formula prediction
  library(R.utils)
  library(webchem)
  library(data.table)  # More efficient than base R
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(stringdist)
  library(data.table)
  library(BiocParallel)
  library(Spectra)
  library(data.table)
  library(S4Vectors)
  library(mzR)
}

{
  # Read HMDB Annotation CSV and mz(X)ML from processing folder direction provided by the SHINY App.
  {
    HMDB_Metabolites <- fread(file_path)
    
    # List mzXML / mzML files in processing folder
    mzxml_files <- list.files(data_dir_raw_files, pattern = "\\.(mzXML|mzML)$", full.names = TRUE)
    
    # Define common m/z grid
    common_mz_grid <- seq(from = 70, to = 1000, by = 0.0005)
    
    total_files <- length(mzxml_files)
    processed_files <- 0
    
    message("Top-N disabled. All files will be processed using RT-window resampling.")
    
    # ---- PROCESSING FUNCTION for RT-WINDOW defined in the SHINY App ----
    process_by_rt_window <- function(file, common_mz_grid, rt_min, rt_max) {
      raw_data <- readMSData(files = file, mode = "onDisk")
      ms1_data <- filterMsLevel(raw_data, msLevel = 1)
      spectra_list <- spectra(ms1_data)
      retention_times <- rtime(ms1_data)
      
      valid_indices <- which(retention_times >= rt_min & retention_times <= rt_max)
      selected_spectra <- spectra_list[valid_indices]
      
      if (length(selected_spectra) == 0)
        return(rep(0, length(common_mz_grid)))
      
      resampled_intensities <- Reduce("+", lapply(selected_spectra, function(spectrum) {
        approx(
          x = mz(spectrum),
          y = intensity(spectrum),
          xout = common_mz_grid,
          method = "linear",
          yleft = 0,
          yright = 0
        )$y
      }))
      
      return(resampled_intensities)
    }
    
    # Storage list
    resampled_data_list <- list()
    resampled_data_file <- file.path(data_dir_raw_files, "resampled_data.rds")
    
    if (file.exists(resampled_data_file)) {
      message("Resampled data RDS already exists. Loading from file...")
      resampled_data_list <- readRDS(resampled_data_file)
      
    } else {
      message("No RDS found. Processing mzML/mzXML files using RT-window...")
      
      for (file in mzxml_files) {
        processed_files <- processed_files + 1
        message(sprintf("Processing file %d/%d: %s",
                        processed_files, total_files, basename(file)))
        
        intensities <- process_by_rt_window(file, common_mz_grid, rt_min = rt_min, rt_max = rt_max)
        
        sample_name <- basename(tools::file_path_sans_ext(file))
        resampled_data_list[[sample_name]] <- data.frame(
          mz = common_mz_grid,
          intensity = intensities
        )
      }
      
      saveRDS(resampled_data_list,
              file = file.path(data_dir_raw_files, "resampled_data.rds"))
      
      message("Saved RT-window resampled data to RDS.")
    }
  }
}

#Export the averaged spectrum as mzml for usage of external programs like ToppView etc.
if (EXPORT_SPECTRA) {

  
  export_dir <- file.path(data_dir, "accumulated_mzML")
  dir.create(export_dir, showWarnings = FALSE)
  
  for (sample_name in names(resampled_data_list)) {
    spectrum_df <- resampled_data_list[[sample_name]]
    
    # Ensure numeric
    spectrum_df <- as.data.table(spectrum_df)
    spectrum_df[, mz := as.numeric(mz)]
    spectrum_df[, intensity := as.numeric(intensity)]
    
    # Remove zero intensities
    spectrum_df <- spectrum_df[intensity > 0, ]
    
    if (nrow(spectrum_df) == 0) {
      message("Skipping empty spectrum for sample: ", sample_name)
      next
    }
    
    mzs <- spectrum_df$mz
    ints <- spectrum_df$intensity
    tic_value <- sum(ints)
    
    spectra_data <- DataFrame(
      mz = I(list(mzs)),
      intensity = I(list(ints)),
      rtime = 10,
      msLevel = 1L,
      centroided = FALSE,
      fromFile = 1L,
      scanWindowLowerLimit = 70,
      scanWindowUpperLimit = 1000,
      tic = tic_value,
      scanIndex = 1L,
      polarity = 1L
    )
    
    sps <- Spectra(spectra_data, backend = MsBackendMemory())
    
    outfile <- file.path(export_dir, paste0(sample_name, "_accumulated_profile_data.mzML"))
    export(sps, MsBackendMzR(), file = outfile)
    message("Exported: ", outfile, " with TIC = ", tic_value)
  }
} else {
  message("mzML export disabled (EXPORT_SPECTRA == FALSE)")
}


  
####Optional smoothing section of profile mass spectrum before centroiding. Uses the smoothing options provided by the SHINY App.
  

# =============================
# STEP 1: OPTIONAL SMOOTHING
# =============================

message("===== STEP 1: (Optional) Smoothing =====")

smooth_intensity <- function(intensity, n = smoothing_n, p = smoothing_p) {
  if (!requireNamespace("signal", quietly = TRUE)) install.packages("signal")
  return(signal::sgolayfilt(intensity, p = p, n = n))
}

if (apply_smoothing) {
  message("✔ Smoothing ENABLED — spectra will be updated and used for centroiding.")
  
  for (sample_name in names(resampled_data_list)) {
    spectrum <- resampled_data_list[[sample_name]]
    
    spectrum$intensity <- smooth_intensity(spectrum$intensity)
    
    # *** CRITICAL: overwrite so centroiding uses smoothed data ***
    resampled_data_list[[sample_name]] <- spectrum
  }
  
} else {
  message("❌ Smoothing DISABLED — centroiding will use raw intensities.")
}

gc()


#Centroiding step with settings provided by the SHINY app.

# =============================
# STEP 2: CENTROIDING
# =============================

message("===== STEP 2: Centroiding =====")

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


# === Build combined peak list from (raw OR smoothed) data ===
all_peaks <- rbindlist(lapply(names(resampled_data_list), function(sample_name) {
  spectrum <- resampled_data_list[[sample_name]]
  spectrum %>%
    filter(intensity > intensity_threshold) %>%
    mutate(sample = sample_name)
}), fill = TRUE)


# === Detect peaks === #window_size =2 was use for most analysis, but increase could reduce peak spliting!
local_maxima <- find_local_maxima(all_peaks$mz, all_peaks$intensity, window_size = 2)

mz_sorted <- sort(local_maxima$mz)
int_sorted <- local_maxima$intensity[order(local_maxima$mz)]

# Grouping by ppm tolerance
mz_diffs <- diff(mz_sorted)
mz_tolerances <- get_mz_tolerance(head(mz_sorted, -1), ppm_tolerance)
group <- c(1, 1 + cumsum(mz_diffs > mz_tolerances))

centroided_peaks <- data.frame(mz = mz_sorted, intensity = int_sorted) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise(mz_centroid = sum(mz * intensity) / sum(intensity))


# =============================
# STEP 3: BUILD FEATURE MATRIX
# =============================

message("===== STEP 3: Feature Matrix Generation =====")

resampled_data_dt <- rbindlist(lapply(names(resampled_data_list), function(sample_name) {
  spectrum <- as.data.table(resampled_data_list[[sample_name]])
  spectrum[, sample := sample_name]
}), fill = TRUE)


filtered_peaks <- resampled_data_dt[intensity > intensity_threshold]
centroided_peaks_dt <- as.data.table(centroided_peaks)

centroided_peaks_dt[, ppm_tol := get_mz_tolerance(mz_centroid, ppm_tolerance)]
centroided_peaks_dt[, mz_min := mz_centroid - ppm_tol]
centroided_peaks_dt[, mz_max := mz_centroid + ppm_tol]

filtered_peaks[, mz_start := mz]
filtered_peaks[, mz_end := mz]

setkey(filtered_peaks, mz_start, mz_end)
setkey(centroided_peaks_dt, mz_min, mz_max)

matched_peaks <- foverlaps(
  filtered_peaks,
  centroided_peaks_dt,
  by.x = c("mz_start", "mz_end"),
  by.y = c("mz_min", "mz_max"),
  nomatch = 0
)

feature_matrix_dt_detected <- matched_peaks[, .(intensity = sum(intensity)), 
                                            by = .(sample, mz_centroid)]


# Fill missing (noise) values by forced reintegration to ensure no NA or 0 are present in the data matrix.
all_combinations <- CJ(
  sample = unique(resampled_data_dt$sample),
  mz_centroid = centroided_peaks_dt$mz_centroid
)

all_combinations <- merge(
  all_combinations,
  centroided_peaks_dt[, .(mz_centroid, mz_min, mz_max)],
  by = "mz_centroid",
  all.x = TRUE
)

integrated_noise <- resampled_data_dt[
  all_combinations,
  on = .(sample, mz >= mz_min, mz <= mz_max),
  allow.cartesian = TRUE
][, .(intensity = sum(intensity, na.rm = TRUE)), by = .(sample, mz_centroid)]


combined_dt <- merge(
  integrated_noise,
  feature_matrix_dt_detected,
  by = c("sample", "mz_centroid"),
  all.x = TRUE,
  suffixes = c("_noise", "_detected")
)

combined_dt[, intensity := ifelse(!is.na(intensity_detected), intensity_detected, intensity_noise)]

feature_matrix_df <- dcast(
  combined_dt[, .(sample, mz_centroid, intensity)],
  sample ~ mz_centroid,
  value.var = "intensity",
  fill = NA
)

rownames(feature_matrix_df) <- feature_matrix_df$sample
feature_matrix_df$sample <- NULL

message("✔ Feature matrix complete!")

gc()

# =============================
# EXPORT: SMOOTHED / PROFILE DATA
# =============================
if (EXPORT_SPECTRA) {
  message("📤 Exporting smoothed/profile spectra (EXPORT_SPECTRA = TRUE)")
  
  library(Spectra)
  library(data.table)
  library(S4Vectors)
  library(mzR)
  
  export_dir <- file.path(data_dir, "accumulated_smoothed_mzML")
  dir.create(export_dir, showWarnings = FALSE)
  
  for (sample_name in names(resampled_data_list)) {
    spectrum_df <- as.data.table(resampled_data_list[[sample_name]])
    spectrum_df[, mz := as.numeric(mz)]
    spectrum_df[, intensity := as.numeric(intensity)]
    spectrum_df <- spectrum_df[intensity > 0]
    
    if (nrow(spectrum_df) == 0) next
    
    mzs <- spectrum_df$mz
    ints <- spectrum_df$intensity
    tic_value <- sum(ints)
    
    spectra_data <- DataFrame(
      mz = I(list(mzs)),
      intensity = I(list(ints)),
      rtime = 10,
      msLevel = 1L,
      centroided = FALSE,
      fromFile = 1L
    )
    
    sps <- Spectra(spectra_data, backend = MsBackendMemory())
    
    outfile <- file.path(export_dir, paste0(sample_name, "_smoothed_profile.mzML"))
    export(sps, MsBackendMzR(), file = outfile)
    message("✔ Exported smoothed profile: ", outfile)
  }
} else {
  message("⛔ Skipping profile export (EXPORT_SPECTRA = FALSE)")
}


# =============================
# EXPORT: CENTROIDED DATA
# =============================
if (EXPORT_SPECTRA) {
  message("📤 Exporting centroided spectra (EXPORT_SPECTRA = TRUE)")
  
  centroided_export_dir <- file.path(data_dir, "centroided_mzML")
  dir.create(centroided_export_dir, showWarnings = FALSE)
  
  for (sample_name in unique(combined_dt$sample)) {
    sample_peaks <- combined_dt[sample == sample_name & intensity > 0]
    if (nrow(sample_peaks) == 0) next
    
    mzs <- as.numeric(sample_peaks$mz_centroid)
    ints <- as.numeric(sample_peaks$intensity)
    tic_value <- sum(ints, na.rm = TRUE)
    
    spectra_data <- DataFrame(
      mz = I(list(mzs)),
      intensity = I(list(ints)),
      rtime = 10,
      msLevel = 1L,
      centroided = TRUE,
      fromFile = 1L
    )
    
    sps <- Spectra(spectra_data, backend = MsBackendMemory())
    outfile <- file.path(centroided_export_dir, paste0(sample_name, "_centroided_LM.mzML"))
    export(sps, MsBackendMzR(), file = outfile)
    
    message("✔ Exported centroided spectrum for ", sample_name, " (", length(mzs), " peaks)")
  }
} else {
  message("⛔ Skipping centroid export (EXPORT_SPECTRA = FALSE)")
}



gc()
  
  
###
##### M/z range limiter to reduce computational load. Molecules greater than 850Da require high computational load due to the increasing number of formulas being calculated.
###
  {
    mz_min_limit <- 70
    mz_max_limit <- 850
    
    sample_names <- rownames(feature_matrix_df)
    mz_values <- suppressWarnings(as.numeric(colnames(feature_matrix_df)))
    valid_cols <- which(!is.na(mz_values) & mz_values >= mz_min_limit & mz_values <= mz_max_limit)
    
    feature_matrix_df <- feature_matrix_df[, valid_cols, with = FALSE]
    rownames(feature_matrix_df) <- sample_names
    message("m/z range set to ", mz_min_limit, " - ", mz_max_limit)
}
  

  ###
  #### Isotope clustering based on 13C isotope pattern with a set tolerance of 10ppm.
  ###
  {
    mass_diff <- 1.0033548378
    ppm <- 10
    
    detect_isotopic_groups_with_numbers <- function(neutral_masses, intensities, mass_diff = 1.0033548378, ppm = 10) {
      neutral_masses <- sort(neutral_masses)
      intensities <- intensities[order(neutral_masses)]
      
      labels <- rep("monoisotopic", length(neutral_masses))
      group_numbers <- rep(NA, length(neutral_masses))
      
      group_id <- 1
      used_indices <- rep(FALSE, length(neutral_masses))
      
      for (i in seq_along(neutral_masses)) {
        if (used_indices[i]) next
        
        base_mass <- neutral_masses[i]
        base_intensity <- intensities[i]
        tol <- base_mass * ppm / 1e6
        
        expected_m1 <- base_mass + mass_diff
        m1_index <- which(abs(neutral_masses - expected_m1) <= tol & !used_indices)
        
        if (length(m1_index) > 0 && intensities[m1_index[1]] < base_intensity) {
          labels[i] <- "monoisotopic"
          labels[m1_index[1]] <- "M+1"
          group_numbers[c(i, m1_index[1])] <- group_id
          used_indices[c(i, m1_index[1])] <- TRUE
          
          expected_m2 <- base_mass + 2 * mass_diff
          m2_index <- which(abs(neutral_masses - expected_m2) <= tol & !used_indices)
          if (length(m2_index) > 0 &&
              intensities[m2_index[1]] < base_intensity &&
              intensities[m2_index[1]] < intensities[m1_index[1]]) {
            labels[m2_index[1]] <- "M+2"
            group_numbers[m2_index[1]] <- group_id
            used_indices[m2_index[1]] <- TRUE
          }
          
          group_id <- group_id + 1
          next
        }
        
        labels[i] <- "monoisotopic"
        group_numbers[i] <- group_id
        used_indices[i] <- TRUE
        group_id <- group_id + 1
      }
      
      return(list(labels = labels, group_numbers = group_numbers))
    }
    message("Isotope groups detected with mass difference set to ", mass_diff, " Da and mass tolerance set to ", ppm, " ppm")
  }
  

###
#### Adduct-aware neutral mass calculation. Assumes every peak to be a (de)-protonated species and searches for other possible matching adducts.
###

if (Ionization_method == "Positive") {
  
  message("Calculation of Molecular Forumla using RDISOP (Ionization_method = Positive)")
  
  proton_mass <- 1.007276  # [M+H]+
  na_mass     <- 22.989221 # [M+Na]+
  k_mass      <- 38.963158 # [M+K]+
  ppm <- 10
  
  detect_neutral_masses <- function(mz_values, ppm = 10) {
    adduct_masses <- list(H = proton_mass, Na = na_mass, K = k_mass)
    
    neutral_masses <- rep(NA, length(mz_values))
    used <- rep(FALSE, length(mz_values))
    
    for (i in seq_along(mz_values)) {
      if (used[i]) next
      
      # Neutral mass candidate for positive mode
      candidate_neutral <- mz_values[i] - proton_mass
      
      # ppm tolerance
      tol <- mz_values[i] * ppm / 1e6
      
      # expected adduct masses
      na_idx <- which(abs(mz_values - (candidate_neutral + na_mass)) <= tol & !used)
      k_idx  <- which(abs(mz_values - (candidate_neutral + k_mass)) <= tol & !used)
      
      neutral_masses[i] <- candidate_neutral
      used[i] <- TRUE
      
      if (length(na_idx) > 0) { neutral_masses[na_idx[1]] <- candidate_neutral; used[na_idx[1]] <- TRUE }
      if (length(k_idx)  > 0) { neutral_masses[k_idx[1]]  <- candidate_neutral; used[k_idx[1]]  <- TRUE }
    }
    
    return(neutral_masses)
  }
  
} else if (Ionization_method == "Negative") {
  
  message("Calculation of Molecular Forumla using RDISOP (Ionization_method = Negative)")
  
  proton_mass <- 1.007276 #M-H
  cl_mass <- 34.969401 #M-Cl
  fa_mass <- 44.998201 #M-H+FA
  ppm <- 10 
  
  detect_neutral_masses <- function(mz_values, ppm = 10) {
    
    neutral_masses <- rep(NA, length(mz_values))
    used <- rep(FALSE, length(mz_values))
    
    for (i in seq_along(mz_values)) {
      if (used[i]) next
      
      # Assume [M-H]- to get candidate neutral
      candidate_neutral <- mz_values[i] + proton_mass
      
      tol <- mz_values[i] * ppm / 1e6
      
      # check [M+Cl]-
      expected_cl <- candidate_neutral + cl_mass
      cl_idx <- which(abs(mz_values - expected_cl) <= tol & !used)
      
      # check [M+FA]-
      expected_fa <- candidate_neutral + fa_mass
      fa_idx <- which(abs(mz_values - expected_fa) <= tol & !used)
      
      neutral_masses[i] <- candidate_neutral
      used[i] <- TRUE
      
      if (length(cl_idx) > 0) {
        neutral_masses[cl_idx[1]] <- candidate_neutral
        used[cl_idx[1]] <- TRUE
      }
      if (length(fa_idx) > 0) {
        neutral_masses[fa_idx[1]] <- candidate_neutral
        used[fa_idx[1]] <- TRUE
      }
    }
    
    return(neutral_masses)
  }
}
  
###
#### Apply to feature matrix
###
  {
    mz_values <- as.numeric(colnames(feature_matrix_df))
    intensities <- as.numeric(feature_matrix_df[1, ])
    
    # adduct-aware neutral mass calculation (internal only)
    neutral_mass_values <- detect_neutral_masses(mz_values, ppm = ppm)
    
    # isotopic grouping on neutral_mass_values
    result <- detect_isotopic_groups_with_numbers(neutral_mass_values, intensities)
    
    isotope_labels_row <- as.data.frame(t(result$labels), stringsAsFactors = FALSE)
    peak_group_numbers_row <- as.data.frame(t(result$group_numbers), stringsAsFactors = FALSE)
    neutral_masses_row <- as.data.frame(t(neutral_mass_values), stringsAsFactors = FALSE)
    
    colnames(neutral_masses_row) <- colnames(feature_matrix_df)
    colnames(isotope_labels_row) <- colnames(feature_matrix_df)
    colnames(peak_group_numbers_row) <- colnames(feature_matrix_df)
    
    rownames(neutral_masses_row) <- "Neutral Masses"
    rownames(isotope_labels_row) <- "Isotopic Labels"
    rownames(peak_group_numbers_row) <- "Peak Group Numbers"
    
    feature_matrix_df_isotope <- rbind(
      neutral_masses_row,
      isotope_labels_row,
      peak_group_numbers_row,
      feature_matrix_df
    )
    
    new_row_names <- c("Neutral Masses", "Isotopic Labels", "Peak Group Numbers", rownames(feature_matrix_df))
    rownames(feature_matrix_df_isotope) <- new_row_names
    
    # numeric conversion except labels
    for (i in 1:nrow(feature_matrix_df_isotope)) {
      if (!(rownames(feature_matrix_df_isotope)[i] %in% c("Isotopic Labels"))) {
        feature_matrix_df_isotope[i, ] <- as.numeric(feature_matrix_df_isotope[i, ])
      }
    }
    
    # save CSV (same name as before)
    write.csv(feature_matrix_df_isotope,
              file = file.path(data_dir, "feature_matrix_df_isotope.csv"),
              row.names = TRUE)
  }
  




###
####Setting up parameters for RDISOP formula calculation of neutral masses with the use of the (six of) seven golden rules for heuristic filtering.
###

{
  #Molecular formula calculation using RDISOP:
  
  # Define the elements you want to include
  getEssentialElements <- function(neutral_mass) {
    Da <- abs(neutral_mass)
    
    if (Da <= 140) {
      return(initializeElements(c("C", "H", "N", "O", "P", "S")))  # CHNO
    } else if (Da <= 180) { 
      return(initializeElements(c("C", "H", "N", "O", "P", "S")))  # CHNOS
    } else if (Da <= 250) { 
      return(initializeElements(c("C", "H", "N", "O", "P", "S")))  # CHNOS 
    } else if (Da <= 1000) {
      return(initializeElements(c("C", "H", "N", "O", "P", "S")))  # CHNOPS 
    }
  }
  
  #In the analysis of the molecular formula of organic molecules, the degree of unsaturation (DU) (also known as the index of hydrogen deficiency (IHD), double bond equivalents (DBE), or unsaturation index[1]) is a calculation that determines the total number of rings and π bonds.
  #Must be POSITIVE and ideally a numeric integer. 
  
  #Defines mass accuracy based on precursor mass/m/z, allowing a wider ppm value for smaller molecules
  #and stricter ppm values for higher mass/m/z, as they get more computational expensive when allowing more
  #potential formulas. As our ZenoTOF 7600 provides high accuracy measurements, 5-10 ppm is typically more than enough. 
  getmassaccuracy <- function (neutral_mass) {
    Da <- abs(neutral_mass)
    
    if (Da <= 300) {
      return(10)  # CHNO
    } else if (Da <= 1000) { 
      return(5)
    }}
  message("RDISOP formula prediction using the seven golden rules for heuristic filtering, set to 10ppm mass tolerance for molecules < 300 Da and 5 ppm for molecules >300 Da")

  
  #Golden rule #1: Restriction for element numbers
  
  getMaxElementLimits <- function(neutral_mass) {
    Da <- abs(neutral_mass)  # Calculate Da value
    
    # More detail in the lower mass ranges
    if (Da <= 100) {
      return("C7H14N5O4P2S2")  # More elements in lower mass range
    } else if (Da <= 150) {
      return("C10H21N6O6P4S3")  
    } else if (Da <= 200) {
      return("C14H27N6O8P2S3")  
    } else if (Da <= 250) {
      return("C19H35N5O9P2S4")  
    } else if (Da <= 300) {
      return("C22H43N8O10P4S4")  
    } else if (Da <= 350) {
      return("C26H46N9O12P2S4")  
    } else if (Da <= 400) {
      return("C29H54N8O14P3S4")  
    } else if (Da <= 450) {
      return("C35H62N12O15P4S4")  
    } else if (Da <= 500) {
      return("C31H59N8O15P3S2")  
    } else if (Da <= 550) {
      return("C40H69N8O16P3S1")  
    } else if (Da <= 600) {
      return("C41H77N6O17P2S2") 
    } else if (Da <= 650) {
      return("C45H83N6O21P1S3")  
    } else if (Da <= 700) {
      return("C49H89N8O24P6S4")  
    } else if (Da <= 750) {
      return("C51H92N8O17P3S2")  
    } else if (Da <= 800) {
      return("C54H93N10O23P3")  
    } else if (Da <= 850) {
      return("C53H102N10O24P4")  
    } else if (Da <= 900) { #Not used as cutoff is <850 Da
      return("C59H106N7O23P3S1")  
    } else if (Da <= 950) { #Not used as cutoff is <850 Da
      return("C59H102N12O25P5S1")  
    } else if (Da <= 1000) { #Not used as cutoff is <850 Da
      return("C63H112N10O27P6")  
    } else { #Not used as cutoff is <850 Da
      return(NULL)  
    }
  }
  
  getMinElementLimits <- function(neutral_mass) {
    Da <- abs(neutral_mass) 
    
    
    if (Da <= 100) {
      return("C2H2")  
    } else if (Da <= 150) {
      return("C3H3")  
    } else if (Da <= 250) {
      return("C4H4")  
    } else if (Da <= 350) {
      return("C5H10")  
    } else if (Da <= 500) {
      return("C7H14")  
    } else if (Da <= 750) {
      return("C15H30")  
    } else if (Da <= 1000) { #Not used as cutoff is <850 Da
      return("C35H70")  
    } else { #Not used as cutoff is <850 Da
      return(NULL) 
    }
  }
  
  
  # Helper to parse formulas into element counts
  parseFormulaElements <- function(formula) {
    matches <- regmatches(formula, gregexpr("([A-Z][a-z]?)([0-9]*)", formula))[[1]]
    elements <- gsub("[0-9]", "", matches)
    counts <- as.numeric(gsub("[^0-9]", "", matches))
    counts[is.na(counts)] <- 1
    structure(tapply(counts, elements, sum), names = unique(elements))
  }
  
  
  # Lipid class assumption and scoring bonus
  detectLipidClass <- function(formula) {
    el <- parseFormulaElements(formula)
    
    # Extract counts (default to 0 if not present)
    N <- ifelse("N" %in% names(el), el["N"], 0)
    P <- ifelse("P" %in% names(el), el["P"], 0)
    O <- ifelse("O" %in% names(el), el["O"], 0)
    C <- ifelse("C" %in% names(el), el["C"], 0)
    H <- ifelse("H" %in% names(el), el["H"], 0)
    S <- ifelse("S" %in% names(el), el["S"], 0)
    
    lipid_score <- 0
    lipid_class <- "Unknown"
    
    # Phosphatidylcholines (PC) and Phosphatidylethanolamines (PE)
    if (N == 1 && P == 1 && O >= 7 && O <= 9 && S == 0) {
      lipid_score <- 10
      lipid_class <- "PC/PE"
    }
    # Phosphatidylserines (PS)
    else if (N == 1 && P == 1 && O >= 9 && O <= 11 && S == 0) {
      lipid_score <- 10
      lipid_class <- "PS"
    }
    # Phosphatidylglycerols (PG)
    else if (N == 0 && P == 1 && O >= 9 && O <= 11 && S == 0) {
      lipid_score <- 10
      lipid_class <- "PG"
    }
    # Phosphatidylinositols (PI)
    else if (N == 0 && P == 1 && O >= 12 && O <= 14 && S == 0) {
      lipid_score <- 10
      lipid_class <- "PI"
    }
    # Phosphatidic acids (PA)
    else if (N == 0 && P == 1 && O >= 7 && O <= 9 && S == 0) {
      lipid_score <- 9
      lipid_class <- "PA"
    }
    # Ceramides
    else if (N == 2 && P == 0 && O >= 3 && O <= 6 && S == 0) {
      lipid_score <- 8
      lipid_class <- "Cer"
    }
    # Sphingomyelins
    else if (N == 2 && P == 1 && O >= 5 && O <= 7 && S == 0) {
      lipid_score <- 10
      lipid_class <- "SM"
    }
    # Glycerolipids (DG/TG) - typically no N or P
    else if (N == 0 && P == 0 && O >= 4 && O <= 7 && C >= 20 && S == 0) {
      lipid_score <- 7
      lipid_class <- "DG/TG"
    }
    # Cardiolipins (CL) - double phosphate
    else if (N == 0 && P == 2 && O >= 15 && O <= 18 && S == 0) {
      lipid_score <- 10
      lipid_class <- "CL"
    }
    # Lysophospholipids (LPC, LPE, etc.) - similar to PC/PE but fewer carbons
    else if (N == 1 && P == 1 && O >= 6 && O <= 8 && C >= 15 && C <= 30 && S == 0) {
      lipid_score <- 8
      lipid_class <- "LysoPC/LysoPE"
    }
    # Sulfatides (sulfated galactoceramides)
    else if (N == 1 && P == 0 && S == 1 && O >= 8 && O <= 12) {
      lipid_score <- 9
      lipid_class <- "Sulfatide"
    }
    
    # Additional heuristic: lipids typically have high C and H counts
    # and characteristic H/C ratios (usually 1.5-2.0)
    if (C > 0) {
      hc_ratio <- H / C
      if (C >= 20 && hc_ratio >= 1.5 && hc_ratio <= 2.0) {
        lipid_score <- lipid_score + 2
      }
    }
    
    return(list(score = lipid_score, class = lipid_class))
  }
  
  # Extract neutral mass values and intensities from feature_matrix_df_isotope
  neutral_mass_values <- as.numeric(feature_matrix_df_isotope[1, ])
  intensity_values <- as.numeric(feature_matrix_df_isotope[4, ])  # Assuming 4th row is the first sample
  
  # Call the clustering function with both arguments
  peak_groups <- detect_isotopic_groups_with_numbers(
    neutral_masses = neutral_mass_values,
    intensities = intensity_values
  )$group_numbers
  
  # Get unique peak groups
  unique_groups <- unique(na.omit(peak_groups))
  
  
  # Create additional rows to store the best molecular information
  molecular_info_row_formula <- rep(NA, ncol(feature_matrix_df))
  molecular_info_row_score <- rep(NA, ncol(feature_matrix_df))
  molecular_info_row_valid <- rep(NA, ncol(feature_matrix_df))
  molecular_info_row_hcratio <- rep(NA, ncol(feature_matrix_df))
  molecular_info_row_ntoc <- rep(NA, ncol(feature_matrix_df))
  molecular_info_row_otoc <- rep(NA, ncol(feature_matrix_df))
  molecular_info_row_ptoc <- rep(NA, ncol(feature_matrix_df))
  molecular_info_row_stoc <- rep(NA, ncol(feature_matrix_df))
  molecular_info_row_heteroprob <- rep(NA, length(mz_values))
  molecular_info_row_dbe <- rep(NA, ncol(feature_matrix_df))
  molecular_info_row_valence <- rep(NA, ncol(feature_matrix_df))
  molecular_info_row_lipid_class <- rep(NA, ncol(feature_matrix_df))
  molecular_info_row_lipid_score <- rep(NA, ncol(feature_matrix_df))
  
  
  # Initialize molecular data storage (to store all calculated formulas in a separate CSV)
  molecular_data <- data.frame(
    Peak_Group = integer(),
    Monoisotopic_Mass = numeric(),
    Mono_Intensity = numeric(),
    M1_Mass = numeric(),
    M1_Intensity = numeric(),
    M2_Mass = numeric(),
    M2_Intensity = numeric(),
    Molecular_Formula = character(),
    Score = numeric(),
    Valid = logical(),
    stringsAsFactors = FALSE
  )
}


for (group in unique_groups) {
  indices <- which(peak_groups == group)
  mono_index <- indices[1]
  mono_mass <- neutral_mass_values[mono_index]
  
  max_elements <- getMaxElementLimits(mono_mass)
  min_elements <- getMinElementLimits(mono_mass)
  essential_elements <- getEssentialElements(mono_mass)  
  mz_dependent_mass_accuracy <- getmassaccuracy(mono_mass)
  
  
  if (!is.null(max_elements)) {
    masses <- neutral_mass_values[indices]
    intensities <- colMeans(feature_matrix_df[, ..indices, drop = FALSE], na.rm = TRUE)
    
    molecules <- tryCatch(
      decomposeIsotopes(masses, intensities, ppm = mz_dependent_mass_accuracy, z = 0, 
                        elements = essential_elements, 
                        minElements = min_elements,
                        maxElements = max_elements,
                        maxisotopes = 3),
      error = function(e) return(NULL)
    )
    
    if (!is.null(molecules)) {
      formulas <- getFormula(molecules)
      scores <- getScore(molecules)
      validities <- getValid(molecules)
      
      scores <- as.numeric(scores)
      scores[is.na(scores)] <- -Inf  
      
      #Golden Rule #4:Hydrogen/Carbon element ratios for all formulas
      hc_ratios <- sapply(formulas, function(f) {
        el <- parseFormulaElements(f)
        if (!("C" %in% names(el)) || el["C"] == 0) return(NA)
        el["H"] / el["C"]
      })
      
      #Golden Rule #5: Heteroatom check
      
      # Calculate heteroatom-to-carbon ratios
      hetero_ratios <- lapply(formulas, function(f) {
        el <- parseFormulaElements(f)
        c_count <- ifelse("C" %in% names(el), el["C"], NA)
        if (is.na(c_count) || c_count == 0) return(rep(NA, 4))
        c(
          NtoC = ifelse("N" %in% names(el), el["N"] / c_count, 0),
          OtoC = ifelse("O" %in% names(el), el["O"] / c_count, 0),
          PtoC = ifelse("P" %in% names(el), el["P"] / c_count, 0),
          StoC = ifelse("S" %in% names(el), el["S"] / c_count, 0)
        )
      })
      hetero_ratios <- do.call(rbind, hetero_ratios)
      
      
      #Golden rule 6. Element probability
      checkElementProbability <- function(elements) {
        # Extract only N, O, P, S and make sure there are no missing values
        required_elements <- c("N", "O", "P", "S")
        # Set missing elements to 0 if they are not present in the formula
        for (el in required_elements) {
          if (!(el %in% names(elements))) {
            elements[el] <- 0  # Assign a default value of 0 if the element is missing
          }
        }
        
        # Filter to only retain required elements (N, O, P, S)
        elements <- elements[names(elements) %in% required_elements]
        
        # Remove NA values and elements with count <= 1
        r1 <- na.omit(elements)
        r2 <- r1[r1 > 1]
        
        # If there are two or fewer heteroatoms with a high count, we don't apply the restriction
        if (length(r2) <= 2) return(TRUE)
        
        # If all four heteroatoms are present with high counts, apply a strict rule
        if (length(r2) == 4 && 
            (r2["N"] >= 10 || r2["O"] >= 20 || r2["P"] >= 4 || r2["S"] >= 3)) {
          return(FALSE)
        }
        
        # Apply rules based on combinations of N, O, P, S
        combo <- paste(sort(names(r2)), collapse = "")
        switch(combo,
               NOP = {
                 if (min(r2) > 3 && (r2["N"] >= 11 || r2["O"] >= 22 || r2["P"] >= 7)) return(FALSE)
               },
               OPS = {
                 if (min(r2) > 1 && (r2["S"] >= 3 || r2["O"] >= 14 || r2["P"] >= 3)) return(FALSE)
               },
               NPS = {
                 if (min(r2) > 1 && (r2["N"] >= 4 || r2["S"] >= 3 || r2["P"] >= 3)) return(FALSE)
               },
               NOS = {
                 if (min(r2) > 6 && (r2["S"] >= 8 || r2["O"] >= 14 || r2["N"] >= 19)) return(FALSE)
               }
        )
        
        return(TRUE)
      }
      
      
      # In the loop where you process each formula
      element_probabilities <- sapply(formulas, function(f) {
        el <- parseFormulaElements(f)
        checkElementProbability(el)
      })
      
      
      # DBE calculator based on extended CHNOPS formula
      calculateDBE <- function(formula) {
        el <- parseFormulaElements(formula)
        C <- ifelse("C" %in% names(el), el["C"], 0)
        H <- ifelse("H" %in% names(el), el["H"], 0)
        N <- ifelse("N" %in% names(el), el["N"], 0)
        X <- ifelse("Halogen" %in% names(el), el["Halogen"], 0)  # Optional, if halogens included
        P <- ifelse("P" %in% names(el), el["P"], 0)
        S <- ifelse("S" %in% names(el), el["S"], 0)

        dbe <- C - H / 2 + N / 2 + P / 2 + 1
        return(dbe)
      }
      # Compute DBE values for each formula
      dbe_values <- sapply(formulas, calculateDBE)
      
      # Only allow positive or 0 (0 are no DB or rings, and ideally integer or near-integer) DBEs. 0.5 is usually a protonated species if used on mz values instead of neutral masses.
      valid_dbe <- !is.na(dbe_values) & dbe_values >= 0
      
      
      
      # Golden Rule 7: Valence check
      checkValenceCompliance <- function(formula) {
        el <- parseFormulaElements(formula)
        
        # Define typical valences for each element
        valences <- c(
          H = 1, C = 4, N = 3, O = 2, P = 3, S = 2, Na = 1
        )
        
        # Sum of electrons used for bonding based on valence and atom count
        total_bonds <- 0
        
        for (element in names(el)) {
          if (!element %in% names(valences)) {
            return(FALSE)  # Unknown element, reject
          }
          total_bonds <- total_bonds + valences[element] * el[element]
        }
        
        return(total_bonds %% 2 == 0)
      }
      
      valid_valence <- sapply(formulas, checkValenceCompliance)
      
      # Lipid detection and scoring
      lipid_info <- lapply(formulas, detectLipidClass)
      lipid_scores <- sapply(lipid_info, function(x) x$score)
      lipid_classes <- sapply(lipid_info, function(x) x$class)
      
      # Adjust final scores by adding lipid bonus
      adjusted_scores <- scores + lipid_scores
      
    
      # Use these probabilities in your molecular_data dataframe
      temp_df_all <- data.frame(
        Peak_Group = rep(group, length(formulas)),
        Monoisotopic_Mass = rep(mono_mass, length(formulas)),
        Mono_Intensity = rep(intensities[1], length(formulas)),
        M1_Mass = ifelse(length(masses) > 1, rep(masses[2], length(formulas)), NA),
        M1_Intensity = ifelse(length(intensities) > 1, rep(intensities[2], length(formulas)), NA),
        M2_Mass = ifelse(length(masses) > 2, rep(masses[3], length(formulas)), NA),
        M2_Intensity = ifelse(length(intensities) > 2, rep(intensities[3], length(formulas)), NA),
        Molecular_Formula = formulas,
        Score = scores,
        Adjusted_Score = adjusted_scores,
        Valid = validities,
        DBE = dbe_values,
        HC_Ratio = hc_ratios,
        NtoC = hetero_ratios[, "NtoC"],
        OtoC = hetero_ratios[, "OtoC"],
        PtoC = hetero_ratios[, "PtoC"],
        StoC = hetero_ratios[, "StoC"],
        Heteroatom_probability = element_probabilities,
        Valence_Valid = valid_valence,
        Lipid_Score = lipid_scores,
        Lipid_Class = lipid_classes,
        stringsAsFactors = FALSE
      )
      
      molecular_data <- rbind(molecular_data, temp_df_all)
      
      #Golden Rule #4: Hydrogen/Carbon element ratio check. Adapted from calcMF package
      # Filter by H/C ratio
      valid_hc <- is.na(hc_ratios) | (hc_ratios >= 0.2 & hc_ratios <= 3.1)
      valid_indices <- which(valid_hc)
      
      # Golden Rule #5: Heteroatom filter
      # NtoC < 1.3; OtoC < 1.2; PtoC < 0.3; StoC < 0.8))
      valid_hetero <- (
        (is.na(hetero_ratios[, "NtoC"]) | hetero_ratios[, "NtoC"] < 1.3) &
          (is.na(hetero_ratios[, "OtoC"]) | hetero_ratios[, "OtoC"] < 1.2) &
          (is.na(hetero_ratios[, "PtoC"]) | hetero_ratios[, "PtoC"] < 0.3) &
          (is.na(hetero_ratios[, "StoC"]) | hetero_ratios[, "StoC"] < 0.8)
      )
      
      valid_indices <- which(valid_hc & valid_hetero & valid_dbe & valid_valence)
      
      if (length(valid_indices) > 0) {
        # Use adjusted_scores instead of scores for best formula selection
        best_index <- valid_indices[which.max(adjusted_scores[valid_indices])]
        molecular_info_row_formula[mono_index] <- formulas[best_index]
        molecular_info_row_score[mono_index] <- format(adjusted_scores[best_index], scientific = TRUE)
        molecular_info_row_valid[mono_index] <- validities[best_index]
        molecular_info_row_hcratio[mono_index] <- hc_ratios[best_index]
        molecular_info_row_ntoc[mono_index] <- hetero_ratios[best_index, "NtoC"]
        molecular_info_row_otoc[mono_index] <- hetero_ratios[best_index, "OtoC"]
        molecular_info_row_ptoc[mono_index] <- hetero_ratios[best_index, "PtoC"]
        molecular_info_row_stoc[mono_index] <- hetero_ratios[best_index, "StoC"]
        molecular_info_row_heteroprob[mono_index] <- element_probabilities[best_index]
        molecular_info_row_dbe[mono_index] <- dbe_values[best_index]
        molecular_info_row_valence[mono_index] <- valid_valence[best_index]
        molecular_info_row_lipid_class[mono_index] <- lipid_classes[best_index]
        molecular_info_row_lipid_score[mono_index] <- lipid_scores[best_index]
      }
    }
  }
}

###
####Convert peak group RDISOP information to df for csv export
###
{
  # Convert to data frames
  molecular_info_row_formula <- as.data.frame(t(molecular_info_row_formula), stringsAsFactors = FALSE)
  molecular_info_row_score <- as.data.frame(t(molecular_info_row_score), stringsAsFactors = FALSE)
  molecular_info_row_valid <- as.data.frame(t(molecular_info_row_valid), stringsAsFactors = FALSE)
  molecular_info_row_dbe <- as.data.frame(t(molecular_info_row_dbe), stringAsFactors = FALSE)
  molecular_info_row_hcratio <- as.data.frame(t(molecular_info_row_hcratio), stringsAsFactors = FALSE)  # <-- Convert H/C ratio row
  molecular_info_row_ntoc <- as.data.frame(t(molecular_info_row_ntoc), stringsAsFactors = FALSE)
  molecular_info_row_otoc <- as.data.frame(t(molecular_info_row_otoc), stringsAsFactors = FALSE)
  molecular_info_row_ptoc <- as.data.frame(t(molecular_info_row_ptoc), stringsAsFactors = FALSE)
  molecular_info_row_stoc <- as.data.frame(t(molecular_info_row_stoc), stringsAsFactors = FALSE)
  molecular_info_row_heteroprob <- as.data.frame(t(molecular_info_row_heteroprob), stringsAsFactors = FALSE)
  molecular_info_row_valence <- as.data.frame(t(molecular_info_row_valence), stringsAsFactors = FALSE)
  
  
  # Ensure column names match
  colnames(molecular_info_row_formula) <- colnames(feature_matrix_df)
  colnames(molecular_info_row_score) <- colnames(feature_matrix_df)
  colnames(molecular_info_row_valid) <- colnames(feature_matrix_df)
  colnames(molecular_info_row_dbe) <- colnames(feature_matrix_df)
  colnames(molecular_info_row_hcratio) <- colnames(feature_matrix_df)  # <-- Set column names for H/C row
  colnames(molecular_info_row_ntoc) <- colnames(feature_matrix_df)
  colnames(molecular_info_row_otoc) <- colnames(feature_matrix_df)
  colnames(molecular_info_row_ptoc) <- colnames(feature_matrix_df)
  colnames(molecular_info_row_stoc) <- colnames(feature_matrix_df)
  colnames(molecular_info_row_heteroprob) <- colnames(feature_matrix_df)
  colnames(molecular_info_row_valence) <- colnames(feature_matrix_df)
  
  
  feature_matrix_df_final <- rbind(
    feature_matrix_df_isotope[1:3, ],
    molecular_info_row_formula,
    molecular_info_row_score,
    molecular_info_row_valid,
    molecular_info_row_dbe,
    molecular_info_row_valence,
    molecular_info_row_hcratio,
    molecular_info_row_ntoc,
    molecular_info_row_otoc,
    molecular_info_row_ptoc,
    molecular_info_row_stoc,
    molecular_info_row_heteroprob,  # <-- Add this line
    feature_matrix_df_isotope[-(1:3), ]
  )
  
  
  rownames(feature_matrix_df_final)[4:14] <- c(
    "Molecular_Formula", "Score", "Validity", "DBE", "Valence_Valid", "H/C ratio",
    "N/C ratio", "O/C ratio", "P/C ratio", "S/C ratio", "Heteroatom Probability"  # <-- Added label
  )
  
  
  # Save both dataframes as CSV files
  write.csv(feature_matrix_df_final, file = file.path(data_dir, "feature_matrix_with_molecular_info.csv"), row.names = T)
  write.csv(molecular_data, file = file.path(data_dir, "molecular_data.csv"), row.names = F)
}

message("RDISOP calculations complete. Export of feature_matrix_with_molecular_info.csv and molecular_data.csv successful")

message("Performing annotation of features by matching mz values and predicted formulas using the HMDB database and a mass tolerance of ", ppm_tolerance_annotation, " ppm" )


###############################################################
### ANNOTATION BASED ON MOLECULAR FORMULA using the provided HMDB CSV defined in the SHINY APP
###############################################################

feature_matrix_df_final_annotated <- feature_matrix_df_final

library(dplyr)

# Normalize molecular formula
feature_matrix_df_final_annotated["Molecular_Formula", ] <- tolower(trimws(as.character(unlist(feature_matrix_df_final_annotated["Molecular_Formula", ]))))
HMDB_Metabolites$CHEMICAL_FORMULA <- trimws(tolower(HMDB_Metabolites$CHEMICAL_FORMULA))

# Match formula -> NAME
feature_matrix_df_final_annotated["Formula_Matched_NAME", ] <- sapply(feature_matrix_df_final_annotated["Molecular_Formula", ], function(MF) {
  if (is.na(MF) || MF == "") return(NA)
  match_idx <- which(HMDB_Metabolites$CHEMICAL_FORMULA == MF)
  if (length(match_idx) > 0) paste(HMDB_Metabolites$NAME[match_idx], collapse = ", ") else NA
})

# Match formula -> HMDB ID
feature_matrix_df_final_annotated["Formula_Matched_HMDB", ] <- sapply(feature_matrix_df_final_annotated["Molecular_Formula", ], function(MF) {
  if (is.na(MF) || MF == "") return(NA)
  match_idx <- which(HMDB_Metabolites$CHEMICAL_FORMULA == MF)
  if (length(match_idx) > 0) paste(HMDB_Metabolites$HMDB[match_idx], collapse = ", ") else NA
})

# Match formula -> SMILES
feature_matrix_df_final_annotated["Formula_Matched_SMILES", ] <- sapply(feature_matrix_df_final_annotated["Molecular_Formula", ], function(MF) {
  if (is.na(MF) || MF == "") return(NA)
  match_idx <- which(HMDB_Metabolites$CHEMICAL_FORMULA == MF)
  if (length(match_idx) > 0) paste(HMDB_Metabolites$SMILES[match_idx], collapse = ", ") else NA
})

# Convert any list columns to character
list_columns <- sapply(feature_matrix_df_final_annotated, is.list)
if (any(list_columns)) {
  feature_matrix_df_final_annotated[list_columns] <- lapply(feature_matrix_df_final_annotated[list_columns], function(x) {
    sapply(x, function(y) if (is.list(y)) paste(unlist(y), collapse = ", ") else y)
  })
}

write.csv(feature_matrix_df_final_annotated,
          file = file.path(data_dir, "feature_matrix_with_formula_annotations.csv"),
          row.names = TRUE)


###############################################################
### ANNOTATION BASED ON m/z based on used ionization mode using the provided HMDB CSV defined in the SHINY APP
###############################################################

# Convert column names (m/z) to numeric
mz_features <- as.numeric(colnames(feature_matrix_df_final_annotated))

# PPM error function
calculate_ppm_error <- function(mz_observed, mz_theoretical) {
  abs((mz_observed - mz_theoretical) / mz_theoretical) * 1e6
}

###############################################
### POSITIVE ION MODE
###############################################
if (Ionization_method == "Positive") {
  
  message("M/z feature matching against HMDB [M+H]+, [M+Na]+ and [M+K]+ (Ionization_method = Positive)")
  
  match_mz <- function(mz_value, hmdb_data, tolerance_ppm) {
    tolerance <- mz_value * tolerance_ppm / 1e6
    matches <- hmdb_data %>%
      filter(
        abs(M.H  - mz_value) <= tolerance |
          abs(M.Na - mz_value) <= tolerance |
          abs(M.K  - mz_value) <= tolerance
      )
    if (nrow(matches) > 0) {
      matches <- matches %>%
        mutate(
          ppm_error_M.H  = ifelse(abs(M.H  - mz_value) <= tolerance, calculate_ppm_error(mz_value, M.H),  NA),
          ppm_error_M.Na = ifelse(abs(M.Na - mz_value) <= tolerance, calculate_ppm_error(mz_value, M.Na), NA),
          ppm_error_M.K  = ifelse(abs(M.K  - mz_value) <= tolerance, calculate_ppm_error(mz_value, M.K),  NA),
          adduct_type = case_when(
            !is.na(ppm_error_M.H) & !is.na(ppm_error_M.Na) & !is.na(ppm_error_M.K) ~ "M+H/M+Na/M+K",
            !is.na(ppm_error_M.H)  ~ "M+H",
            !is.na(ppm_error_M.Na) ~ "M+Na",
            !is.na(ppm_error_M.K)  ~ "M+K",
            TRUE                   ~ NA_character_
          )
        )
    }
    return(matches)
  }
  
  ###############################################
  ### NEGATIVE ION MODE
  ###############################################
} else if (Ionization_method == "Negative") {
  
  message("M/z feature matching against HMDB [M-H]−, [M+HCOO]−, [M+Cl]− (Ionization_method = Negative)")
  
  match_mz <- function(mz_value, hmdb_data, tolerance_ppm) {
    tolerance <- mz_value * tolerance_ppm / 1e6
    matches <- hmdb_data %>%
      filter(
        abs(M..H   - mz_value) <= tolerance |
          abs(M.HCOO - mz_value) <= tolerance |
          abs(M.Cl   - mz_value) <= tolerance
      )
    if (nrow(matches) > 0) {
      matches <- matches %>%
        mutate(
          ppm_error_M..H   = ifelse(abs(M..H   - mz_value) <= tolerance, calculate_ppm_error(mz_value, M..H),   NA),
          ppm_error_M.HCOO = ifelse(abs(M.HCOO - mz_value) <= tolerance, calculate_ppm_error(mz_value, M.HCOO), NA),
          ppm_error_M.Cl   = ifelse(abs(M.Cl   - mz_value) <= tolerance, calculate_ppm_error(mz_value, M.Cl),   NA),
          adduct_type = case_when(
            !is.na(ppm_error_M..H) & !is.na(ppm_error_M.HCOO) & !is.na(ppm_error_M.Cl) ~ "M-H/M+HCOO-/M+Cl-",
            !is.na(ppm_error_M..H)   ~ "M-H",
            !is.na(ppm_error_M.HCOO) ~ "M+HCOO-",
            !is.na(ppm_error_M.Cl)   ~ "M+Cl-",
            TRUE                     ~ NA_character_
          )
        )
    }
    return(matches)
  }
}

###############################################################
### RUN MATCHING ON ALL FEATURES
###############################################################

mz_names        <- vector("character", length(mz_features))
mz_hmdb_ids     <- vector("character", length(mz_features))
mz_adduct_types <- vector("character", length(mz_features))

for (i in seq_along(mz_features)) {
  
  mz_value <- mz_features[i]
  matches  <- match_mz(mz_value, HMDB_Metabolites, ppm_tolerance_annotation)
  
  if (nrow(matches) > 0) {
    
    matched_names <- sapply(1:nrow(matches), function(j) {
      if (Ionization_method == "Positive") {
        ppm_error <- case_when(
          !is.na(matches$ppm_error_M.H[j])  ~ matches$ppm_error_M.H[j],
          !is.na(matches$ppm_error_M.Na[j]) ~ matches$ppm_error_M.Na[j],
          !is.na(matches$ppm_error_M.K[j])  ~ matches$ppm_error_M.K[j],
          TRUE ~ NA_real_
        )
      } else {
        ppm_error <- case_when(
          !is.na(matches$ppm_error_M..H[j])   ~ matches$ppm_error_M..H[j],
          !is.na(matches$ppm_error_M.HCOO[j]) ~ matches$ppm_error_M.HCOO[j],
          !is.na(matches$ppm_error_M.Cl[j])   ~ matches$ppm_error_M.Cl[j],
          TRUE ~ NA_real_
        )
      }
      paste0(matches$NAME[j], " (", sprintf("%.2f", ppm_error), " ppm, ", matches$adduct_type[j], ")")
    })
    
    mz_names[i]        <- paste(matched_names, collapse = ", ")
    mz_hmdb_ids[i]     <- paste(matches$HMDB_ID, collapse = ", ")
    mz_adduct_types[i] <- paste(matches$adduct_type, collapse = ", ")
    
  } else {
    mz_names[i]        <- NA
    mz_hmdb_ids[i]     <- NA
    mz_adduct_types[i] <- NA
  }
}

feature_matrix_df_final_annotated["mz_matched_NAME", ]   <- mz_names
feature_matrix_df_final_annotated["mz_matched_HMDB", ]   <- mz_hmdb_ids
feature_matrix_df_final_annotated["mz_matched_ADDUCT", ] <- mz_adduct_types


###############################################################
### PREDICTED ADDUCT LABELING to return the adduct prediction into the dataframe
###############################################################

# Use the same mz_features already derived from column names above

if (Ionization_method == "Positive") {
  
  predicted_adducts <- rep(NA, length(mz_features))
  used <- rep(FALSE, length(mz_features))
  
  for (i in seq_along(mz_features)) {
    if (used[i]) next
    
    candidate_neutral      <- mz_features[i] - proton_mass
    predicted_adducts[i]   <- "M+H"
    used[i]                <- TRUE
    
    tol <- mz_features[i] * ppm / 1e6
    
    expected_na <- candidate_neutral + na_mass
    na_idx <- which(abs(mz_features - expected_na) <= tol & !used)
    if (length(na_idx) > 0) {
      predicted_adducts[na_idx[1]] <- "M+Na"
      used[na_idx[1]]              <- TRUE
    }
    
    expected_k <- candidate_neutral + k_mass
    k_idx <- which(abs(mz_features - expected_k) <= tol & !used)
    if (length(k_idx) > 0) {
      predicted_adducts[k_idx[1]] <- "M+K"
      used[k_idx[1]]              <- TRUE
    }
  }
  
} else if (Ionization_method == "Negative") {
  
  predicted_adducts <- rep(NA, length(mz_features))
  used <- rep(FALSE, length(mz_features))
  
  for (i in seq_along(mz_features)) {
    if (used[i]) next
    
    candidate_neutral    <- mz_features[i] + proton_mass
    predicted_adducts[i] <- "M-H"
    used[i]              <- TRUE
    
    tol <- mz_features[i] * ppm / 1e6
    
    expected_cl <- candidate_neutral + cl_mass
    cl_idx <- which(abs(mz_features - expected_cl) <= tol & !used)
    if (length(cl_idx) > 0) {
      predicted_adducts[cl_idx[1]] <- "M+Cl"
      used[cl_idx[1]]              <- TRUE
    }
    
    expected_fa <- candidate_neutral + fa_mass
    fa_idx <- which(abs(mz_features - expected_fa) <= tol & !used)
    if (length(fa_idx) > 0) {
      predicted_adducts[fa_idx[1]] <- "M+HCOO"
      used[fa_idx[1]]              <- TRUE
    }
  }
}

adduct_labels_row <- as.data.frame(t(predicted_adducts), stringsAsFactors = FALSE)
colnames(adduct_labels_row) <- colnames(feature_matrix_df_final_annotated)
rownames(adduct_labels_row) <- "Predicted_Adduct"


###############################################################
### ASSEMBLE FINAL ANNOTATED MATRIX
###############################################################

# Collect annotation rows
annotations_df <- feature_matrix_df_final_annotated[
  c("Formula_Matched_NAME",
    "Formula_Matched_HMDB",
    "Formula_Matched_SMILES",
    "mz_matched_NAME",
    "mz_matched_HMDB",
    "mz_matched_ADDUCT"), , drop = FALSE]

# Remove annotation rows from main body
feature_matrix_df_final_annotated <- feature_matrix_df_final_annotated[
  !rownames(feature_matrix_df_final_annotated) %in% rownames(annotations_df), , drop = FALSE]

# Restore Molecular_Formula to uppercase
feature_matrix_df_final_annotated["Molecular_Formula", ] <-
  toupper(trimws(as.character(unlist(feature_matrix_df_final_annotated["Molecular_Formula", ]))))

# Insert Predicted_Adduct at row 3, then annotation rows after row 14
feature_matrix_df_final_annotated <- rbind(
  feature_matrix_df_final_annotated[1:2, , drop = FALSE],
  adduct_labels_row,
  feature_matrix_df_final_annotated[-(1:2), , drop = FALSE]
)

feature_matrix_df_final_annotated <- rbind(
  feature_matrix_df_final_annotated[1:15, ],   # rows 1-2 + Predicted_Adduct + rows 3-14 = 15 rows
  annotations_df,
  feature_matrix_df_final_annotated[-(1:15), ]
)

# Remove any residual list columns
list_columns <- sapply(feature_matrix_df_final_annotated, is.list)
if (any(list_columns)) {
  feature_matrix_df_final_annotated[list_columns] <-
    lapply(feature_matrix_df_final_annotated[list_columns], function(x) {
      sapply(x, function(y) if (is.list(y)) paste(unlist(y), collapse = ", ") else y)
    })
}

cat("Predicted_Adduct row added at row 3.\n")


###############################################################
### EXPORT FULL MATRIX of ANNOTATED AND NON-ANNOTATED FEATURES
###############################################################

write.csv(feature_matrix_df_final_annotated,
          file = file.path(data_dir, "feature_matrix_with_formula_and_mz_annotations.csv"),
          row.names = TRUE)


###############################################################
### EXPORT FILTERED MATRIX (ANNOTATED FEATURES ONLY)
###############################################################

na_cols <- which(
  is.na(feature_matrix_df_final_annotated["Formula_Matched_NAME", ]) &
    is.na(feature_matrix_df_final_annotated["mz_matched_NAME", ])
)

feature_matrix_df_final_annotated_filtered <-
  feature_matrix_df_final_annotated[, -na_cols, drop = FALSE]

write.csv(feature_matrix_df_final_annotated_filtered,
          file = file.path(data_dir, "feature_matrix_with_formula_and_mz_annotations_filtered.csv"),
          row.names = TRUE)
  




#########################################
# Parse annotations to replace m/z column names for METABOANALYST Styled output.
#########################################

library(stringr)

# ------------------------------------------------------------------
# Helper 1: depth-aware comma/semicolon split
# ------------------------------------------------------------------
split_top_level <- function(x) {
  depth   <- 0L
  result  <- character(0)
  current <- ""
  for (ch in strsplit(x, "")[[1]]) {
    if      (ch == "(")                          { depth <- depth + 1L; current <- paste0(current, ch) }
    else if (ch == ")")                          { depth <- depth - 1L; current <- paste0(current, ch) }
    else if (ch %in% c(",", ";") && depth == 0L) { result <- c(result, trimws(current)); current <- "" }
    else                                         { current <- paste0(current, ch) }
  }
  if (nchar(trimws(current)) > 0) result <- c(result, trimws(current))
  result
}

# ------------------------------------------------------------------
# Helper 2: lipid class detection
# ------------------------------------------------------------------
LIPID_CLASS_REGEX <- paste0(
  "^(TG|DG|MG|",
  "PC|PE|PS|PI|PG|PA|",
  "LysoPC|LysoPE|LysoPI|LysoPG|LysoPS|LysoPA|",
  "SM|Cer|GlcCer|GalCer|LacCer|HexCer|",
  "Trihexosylceramide|Tetrahexosylceramide|",
  "Ganglioside(\\s+\\S+)?|",   # allows "Ganglioside GD2("
  "CE)\\s*\\("
)

is_lipid_name <- function(name) {
  grepl(LIPID_CLASS_REGEX, name) && grepl("\\d+:\\d+", name)
}

# ------------------------------------------------------------------
# Helper 3: sum composition converter
# ------------------------------------------------------------------
convert_to_sum_composition <- function(lipid_name) {
  if (is.na(lipid_name) || lipid_name == "") return(lipid_name)
  
  # Normalise spaces before opening paren
  name_norm <- gsub("\\s+\\(", "(", trimws(lipid_name))
  
  # Collapse "Ganglioside GD2(" -> "Ganglioside(" (subtype dropped for sum composition)
  name_norm <- gsub("^(Ganglioside)\\s*\\S+\\s*\\(", "\\1(", name_norm)
  
  # Match CLASS(chain_content)
  m <- regmatches(name_norm, regexec("^([A-Za-z]+)\\((.+)\\)$", name_norm))[[1]]
  if (length(m) != 3) return(lipid_name)
  
  lipid_class   <- m[2]
  chain_content <- m[3]
  
  # Detect ether/plasmalogen prefix (O- or P-) on first chain
  ether_prefix <- ""
  if (grepl("^[OP]-", chain_content)) {
    ether_prefix <- sub("^([OP]-).*", "\\1", chain_content)
  }
  
  # Strip double-bond position annotations: (9Z), (9Z,12Z), etc.
  chain_clean <- gsub("\\([^)]*\\)", "", chain_content)
  
  # Split on / or _ (both used as chain separators in HMDB)
  chains <- trimws(unlist(strsplit(chain_clean, "[/_]")))
  
  total_c  <- 0L
  total_db <- 0L
  
  for (chain in chains) {
    chain <- trimws(chain)
    if (chain == "" || chain == "0:0" || chain == "0") next
    
    # Strip sphingoid base prefix (d, t)
    chain <- gsub("^[dt]", "", chain)
    # Strip ether/plasmalogen prefix
    chain <- gsub("^[OP]-", "", chain)
    
    parts <- strsplit(chain, ":")[[1]]
    if (length(parts) == 2) {
      c_val  <- suppressWarnings(floor(as.numeric(trimws(parts[1]))))
      db_val <- suppressWarnings(floor(as.numeric(trimws(parts[2]))))
      if (!is.na(c_val))  total_c  <- total_c  + c_val
      if (!is.na(db_val)) total_db <- total_db + db_val
    }
  }
  
  if (total_c == 0L) return(lipid_name)
  paste0(lipid_class, " ", ether_prefix, total_c, ":", total_db)
}

# ------------------------------------------------------------------
# Helper 4: main annotation parser
# ------------------------------------------------------------------
parse_first_annotation_and_adduct <- function(name_cell, adduct_cell = NA) {
  if (is.na(name_cell) || name_cell %in% c("", "NA")) return(list(name = NA, adduct = NA))
  
  name_cleaned <- as.character(name_cell)
  name_cleaned <- gsub("\\s*\\([^)]*ppm[^)]*\\)", "", name_cleaned, ignore.case = TRUE)
  name_cleaned <- gsub("\\s*\\d+\\.?\\d*\\s*ppm",  "", name_cleaned, ignore.case = TRUE)
  name_cleaned <- trimws(name_cleaned)
  
  parts      <- split_top_level(name_cleaned)
  name_first <- if (length(parts) > 0) trimws(parts[1]) else name_cleaned
  
  adduct_first <- NA
  
  if (is_lipid_name(name_first)) {
    name_first <- convert_to_sum_composition(name_first)
  } else {
    match <- regmatches(name_first, regexec("^(.+?)\\s*\\((.+)\\)$", name_first))[[1]]
    if (length(match) == 3) {
      name_first   <- trimws(match[2])
      adduct_first <- trimws(match[3])
    }
  }
  
  # Resolve adduct: mz_matched_ADDUCT preferred, fall back to Predicted_Adduct
  if (!is.na(adduct_cell) && adduct_cell != "" && adduct_cell != "NA") {
    adduct_first <- trimws(unlist(strsplit(as.character(adduct_cell), "[,;]"))[1])
  }
  
  list(name = name_first, adduct = adduct_first)
}

# ------------------------------------------------------------------
# Build final annotation labels
# ------------------------------------------------------------------
mz_values <- as.numeric(colnames(feature_matrix_df_final_annotated_filtered))
final_annotation_names <- character(length(mz_values))

for (i in seq_along(mz_values)) {
  
  mz_val           <- round(mz_values[i], 4)
  annotation_found <- FALSE
  
  # Priority 1: Formula_Matched_NAME
  #   Adduct: always Predicted_Adduct — mz_matched_ADDUCT belongs to
  #   a different matching method and must NOT be used here.
  if ("Formula_Matched_NAME" %in% rownames(feature_matrix_df_final_annotated_filtered)) {
    
    formula_adduct_source <- NA
    if ("Predicted_Adduct" %in% rownames(feature_matrix_df_final_annotated_filtered))
      formula_adduct_source <- feature_matrix_df_final_annotated_filtered["Predicted_Adduct", i]
    
    parsed <- parse_first_annotation_and_adduct(
      feature_matrix_df_final_annotated_filtered["Formula_Matched_NAME", i],
      formula_adduct_source
    )
    if (!is.na(parsed$name)) {
      final_annotation_names[i] <- paste0(parsed$name, " (", mz_val, ";",
                                          ifelse(is.na(parsed$adduct), "", parsed$adduct), ")")
      annotation_found <- TRUE
    }
  }
  
  # Priority 2: mz_matched_NAME
  #   Adduct: prefer mz_matched_ADDUCT (its own), fall back to Predicted_Adduct.
  if (!annotation_found && "mz_matched_NAME" %in% rownames(feature_matrix_df_final_annotated_filtered)) {
    
    mz_adduct_source <- NA
    if ("mz_matched_ADDUCT" %in% rownames(feature_matrix_df_final_annotated_filtered))
      mz_adduct_source <- feature_matrix_df_final_annotated_filtered["mz_matched_ADDUCT", i]
    if (is.na(mz_adduct_source) || mz_adduct_source == "" || mz_adduct_source == "NA") {
      if ("Predicted_Adduct" %in% rownames(feature_matrix_df_final_annotated_filtered))
        mz_adduct_source <- feature_matrix_df_final_annotated_filtered["Predicted_Adduct", i]
    }
    
    parsed <- parse_first_annotation_and_adduct(
      feature_matrix_df_final_annotated_filtered["mz_matched_NAME", i],
      mz_adduct_source
    )
    if (!is.na(parsed$name)) {
      final_annotation_names[i] <- paste0(parsed$name, " (", mz_val, ";",
                                          ifelse(is.na(parsed$adduct), "", parsed$adduct), ")")
      annotation_found <- TRUE
    }
  }
  
  # Fallback: m/z only
  if (!annotation_found) final_annotation_names[i] <- as.character(mz_val)
}

cat("Annotation parsing complete:",
    sum(!is.na(final_annotation_names) & final_annotation_names != ""),
    "features annotated\n")


#########################################
# MetaboAnalyst CSV with annotated column names
#########################################

df           <- feature_matrix_df_final_annotated_filtered
colnames(df) <- final_annotation_names
df_clean     <- df[-c(1:21), , drop = FALSE]
sample_names <- rownames(df_clean)

label_vector <- dplyr::case_when(
  grepl("processing",  sample_names, ignore.case = TRUE) ~ "Blank_processing",
  grepl("Blank",       sample_names, ignore.case = TRUE) ~ "Blank",
  grepl("QC",          sample_names, ignore.case = TRUE) ~ "QC",
  grepl("NIST",        sample_names, ignore.case = TRUE) ~ "NIST",
  grepl("Conditioning",sample_names, ignore.case = TRUE) ~ "Conditioning",
  TRUE ~ "Sample"
)

df_clean     <- cbind(label = label_vector, df_clean)
df_clean_out <- cbind(name  = rownames(df_clean), df_clean)
rownames(df_clean_out) <- NULL

write.csv(
  df_clean_out,
  file = file.path(data_dir, "feature_matrix_metabo.csv"),
  row.names = FALSE
)

cat("MetaboAnalyst-ready file exported with annotated column names\n")


#########################################
# Transposed CSV with annotated metabolite row names and class, index and batch annotation.
#########################################

df           <- feature_matrix_df_final_annotated_filtered
colnames(df) <- final_annotation_names
df_clean     <- df[-c(1:21), , drop = FALSE]
df_t         <- t(df_clean)
sample_names <- colnames(df_t)

class_row <- dplyr::case_when(
  grepl("processing",  sample_names, ignore.case = TRUE) ~ "Blank_processing",
  grepl("Blank",       sample_names, ignore.case = TRUE) ~ "Blank",
  grepl("QC",          sample_names, ignore.case = TRUE) ~ "QC",
  grepl("NIST",        sample_names, ignore.case = TRUE) ~ "NIST",
  grepl("Conditioning",sample_names, ignore.case = TRUE) ~ "Conditioning",
  TRUE ~ "Sample"
)

index_row <- seq_along(sample_names)
batch_row <- rep("", length(sample_names))

df_final <- rbind(
  class = class_row,
  index = index_row,
  batch = batch_row,
  df_t
)

df_final_out <- cbind(SampleID = rownames(df_final), df_final)
rownames(df_final_out) <- NULL

write.csv(
  df_final_out,
  file = file.path(data_dir, "feature_matrix_transposed_with_metadata.csv"),
  row.names = FALSE
)

cat("Transposed file exported with annotated metabolite names\n")
message("Annotation feature matrices saved. Processing and annotation completed.")