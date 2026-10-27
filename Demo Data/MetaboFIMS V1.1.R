# MetaboFIMS version: 1.2.0
{
  library(MSnbase)
  library(Rdisop)  
  library(R.utils)
  library(webchem)
  library(data.table) 
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

# =============================
# OUTPUT FOLDERS (set by the Shiny app; sensible fallbacks when the script is run on its own)
#   output_dir      : every result file (CSV) of THIS run goes here  -> <data_dir>/<custom name or date_time>
#   mzml_export_dir : exported mzML spectra go here                  -> the processing folder (data_dir)
#   output_subdirs  : optional named list of subfolders inside output_dir
# =============================
{
  if (!exists("mzml_export_dir") || !nzchar(mzml_export_dir)) mzml_export_dir <- data_dir
  if (!exists("output_dir") || !nzchar(output_dir)) {
    output_dir <- file.path(sub("/+$", "", data_dir), format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
  }
  if (!exists("output_subdirs") || is.null(output_subdirs)) output_subdirs <- list()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Path of a result file: inside output_dir, or inside one of the named subfolders of output_subdirs
  out_file <- function(file, sub = NULL) {
    d <- if (!is.null(sub) && !is.null(output_subdirs[[sub]])) output_subdirs[[sub]] else output_dir
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    file.path(d, file)
  }
  message("Output folder (results / CSV): ", output_dir)
  message("mzML export folder: ", mzml_export_dir)
  
  # Scan polarity written into exported mzML files: 1 = positive, 0 = negative (Spectra / mzR convention)
  export_polarity <- if (identical(get0("Ionization_method", inherits = TRUE, ifnotfound = "Positive"), "Negative")) 0L else 1L
  message("Exported mzML polarity: ", if (export_polarity == 0L) "negative (0)" else "positive (1)")
}

# =============================
# REPRODUCIBILITY: run_parameters.yml (written into every output folder)
#   written at the start (status "started"), extended once the samples are known,
#   closed at the very end of the script (status "completed").
# =============================
{
  rp <- new.env()
  rp$data <- list()
  run_started_time <- Sys.time()
  
  update_run_parameters <- function(new) {
    tryCatch({
      rp$data[names(new)] <- new
      yaml::write_yaml(rp$data, file.path(output_dir, "run_parameters.yml"),
                       handlers = list(logical = yaml::verbatim_logical))
    }, error = function(e) warning("Could not write run_parameters.yml: ", conditionMessage(e), call. = FALSE))
    invisible(NULL)
  }
  
  md5_of <- function(f) {
    f <- if (is.null(f) || length(f) != 1L || is.na(f)) "" else as.character(f)
    if (nzchar(f) && file.exists(f)) unname(tools::md5sum(f)) else NA_character_
  }
  
  get_session_info <- function() {
    si  <- sessionInfo()
    ver <- function(l) as.list(vapply(l, function(p) as.character(p$Version), character(1)))
    list(
      r_version         = si$R.version$version.string,
      platform          = si$platform,
      os                = if (is.null(si$running)) NA_character_ else si$running,
      locale            = si$locale,
      attached_packages = ver(si$otherPkgs),
      loaded_packages   = ver(si$loadedOnly),
      base_packages     = as.list(si$basePkgs)
    )
  }
  
  get_var <- function(nm, default = NA_character_) get0(nm, inherits = TRUE, ifnotfound = default)
  list_or_note <- function(x) if (length(x)) as.list(as.character(x)) else "not set - all HMDB rows used"
  
  settings_names <- c("Ionization_method", "rt_min", "rt_max", "mz_min_limit", "mz_max_limit",
                      "EXPORT_SPECTRA", "apply_smoothing", "smoothing_n", "smoothing_p",
                      "intensity_threshold", "ppm_tolerance", "ppm_tolerance_annotation")
  run_settings <- setNames(lapply(settings_names, get_var), settings_names)
  
  script_path_used <- get_var("pipeline_script_path", "")
  
  update_run_parameters(list(
    run_status      = "started",
    run_started     = format(run_started_time, "%Y-%m-%d %H:%M:%S %z"),
    pipeline_script = list(path = script_path_used, md5 = md5_of(script_path_used)),
    folders         = list(data_dir           = get_var("data_dir"),
                           data_dir_raw_files = get_var("data_dir_raw_files"),
                           output_dir         = output_dir,
                           mzml_export_dir    = mzml_export_dir),
    hmdb            = list(csv         = get_var("file_path"),
                           md5         = md5_of(get_var("file_path", "")),
                           biospecimen = list_or_note(get_var("hmdb_biospecimen", NULL)),
                           origin      = list_or_note(get_var("hmdb_origin", NULL))),
    settings        = run_settings,
    session_info    = get_session_info()
  ))
  message("Run parameters: ", file.path(output_dir, "run_parameters.yml"))
}

{
  {
    HMDB_Metabolites <- fread(file_path)
    
    mzxml_files <- list.files(data_dir_raw_files, pattern = "\\.(mzXML|mzML)$", full.names = TRUE)
    
    common_mz_grid <- seq(from = mz_min_limit, to = mz_max_limit, by = 0.0005)
    
    total_files <- length(mzxml_files)
    processed_files <- 0
    
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
    
    resampled_data_list <- list()
    resampled_data_file <- file.path(data_dir_raw_files, "resampled_data.rds")
    
    if (file.exists(resampled_data_file)) {
      message("Resampled data RDS already exists. Loading from file...")
      resampled_data_list <- readRDS(resampled_data_file)
      resampled_data_source <- "rds_cache"
      
    } else {
      message("No RDS found. Processing mzML/mzXML files using RT-window...")
      resampled_data_source <- "raw_files"
      
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
    
    # ---- Reproducibility: record exactly which samples enter this run ----
    # Names are taken from resampled_data_list, i.e. from the .rds when it was loaded.
    {
      sample_names_used <- names(resampled_data_list)
      current_raw_names <- basename(tools::file_path_sans_ext(mzxml_files))
      samples_info <- list(
        source    = if (resampled_data_source == "rds_cache")
          "resampled_data.rds (cached - raw files were NOT re-read)" else
            "raw mzML/mzXML files (RT-window integration run in this session)",
        n_samples = length(sample_names_used),
        names     = as.list(sample_names_used)
      )
      if (resampled_data_source == "rds_cache") {
        samples_info$rds_file     <- resampled_data_file
        samples_info$rds_modified <- format(file.mtime(resampled_data_file), "%Y-%m-%d %H:%M:%S")
        only_rds    <- setdiff(sample_names_used, current_raw_names)
        only_folder <- setdiff(current_raw_names, sample_names_used)
        samples_info$in_rds_but_not_in_raw_folder <- as.list(only_rds)
        samples_info$in_raw_folder_but_not_in_rds <- as.list(only_folder)
        if (length(only_rds) || length(only_folder)) {
          warning("The cached resampled_data.rds does not match the files currently in the raw folder ",
                  "(see 'samples' in run_parameters.yml). Delete the .rds to reprocess the raw files.",
                  call. = FALSE, immediate. = TRUE)
        }
      } else {
        samples_info$raw_files        <- as.list(basename(mzxml_files))
        samples_info$rds_file_written <- resampled_data_file
      }
      update_run_parameters(list(samples = samples_info))
      message("Samples used (", length(sample_names_used), ", source: ", resampled_data_source, "): ",
              paste(sample_names_used, collapse = ", "))
    }
  }
}

if (EXPORT_SPECTRA) {
  
  
  export_dir <- file.path(mzml_export_dir, "accumulated_mzML")
  dir.create(export_dir, showWarnings = FALSE)
  
  for (sample_name in names(resampled_data_list)) {
    spectrum_df <- resampled_data_list[[sample_name]]
    
    spectrum_df <- as.data.table(spectrum_df)
    spectrum_df[, mz := as.numeric(mz)]
    spectrum_df[, intensity := as.numeric(intensity)]
    
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
      scanWindowLowerLimit = mz_min_limit,
      scanWindowUpperLimit = mz_max_limit,
      tic = tic_value,
      scanIndex = 1L,
      polarity = export_polarity
    )
    
    sps <- Spectra(spectra_data, backend = MsBackendMemory())
    
    outfile <- file.path(export_dir, paste0(sample_name, "_accumulated_profile_data.mzML"))
    export(sps, MsBackendMzR(), file = outfile)
    message("Exported: ", outfile, " with TIC = ", tic_value)
  }
} else {
  message("mzML export disabled (EXPORT_SPECTRA == FALSE)")
}


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
    
    resampled_data_list[[sample_name]] <- spectrum
  }
  
} else {
  message("❌ Smoothing DISABLED — centroiding will use raw intensities.")
}

gc()



# =============================
# STEP 2: CENTROIDING  (per-sample apex detection + sub-grid refinement)
# =============================

message("===== STEP 2: Centroiding =====")

## ---------------- parameters ----------------
grid_spacing           <- median(diff(common_mz_grid))   # 0.0005 Da
centroid_halfwidth_ppm <- 20      # max half-window around an apex, in ppm
centroid_frac          <- 0.50    # points >= 50% of apex enter the centroid
min_samples_feature    <- 1       # keep feature only if apex found in >= n samples
centroid_chunk         <- 20000L  # candidates per matrix block (memory control)

get_mz_tolerance <- function(mz, ppm = ppm_tolerance) (mz * ppm) / 1e6

## ---------------- per-sample apex detection + refinement ----------------
centroid_one_sample <- function(mz, int, min_int,
                                halfwidth_ppm = centroid_halfwidth_ppm,
                                frac          = centroid_frac,
                                spacing       = grid_spacing) {
  n <- length(int)
  if (n < 3L) return(NULL)
  
  ## plateau-safe local maxima on the CONTIGUOUS grid (no pre-filtering)
  cand <- which(int[2:(n - 1L)] >  int[1:(n - 2L)] &
                  int[2:(n - 1L)] >= int[3:n]) + 1L
  cand <- cand[int[cand] > min_int]          # intensity cut AFTER detection
  if (!length(cand)) return(NULL)
  
  parts <- split(seq_along(cand), ceiling(seq_along(cand) / centroid_chunk))
  res   <- vector("list", length(parts))
  
  for (ci in seq_along(parts)) {
    cc <- cand[parts[[ci]]]
    
    w    <- pmax(1L, ceiling(mz[cc] * halfwidth_ppm / 1e6 / spacing))
    W    <- max(w)
    offs <- -W:W
    idx  <- outer(cc, offs, "+")
    ok   <- idx >= 1L & idx <= n
    idx[!ok] <- 1L
    
    I    <- matrix(int[idx], nrow = length(cc))
    M    <- matrix(mz[idx],  nrow = length(cc))
    offm <- matrix(offs, nrow = length(cc), ncol = length(offs), byrow = TRUE)
    inwin    <- ok & abs(offm) <= w
    apex_int <- int[cc]
    
    ## suppress ripple maxima: apex must dominate its own window
    win_max <- rep(-Inf, length(cc))
    for (k in seq_len(ncol(I))) win_max <- pmax(win_max, ifelse(inwin[, k], I[, k], -Inf))
    dom <- apex_int >= win_max
    if (!any(dom)) next
    cc <- cc[dom]; I <- I[dom, , drop = FALSE]; M <- M[dom, , drop = FALSE]
    inwin <- inwin[dom, , drop = FALSE]; apex_int <- apex_int[dom]
    
    ## grow outwards while >= frac*apex and still monotonically descending
    keep <- inwin & (I >= frac * apex_int)
    c0   <- W + 1L
    if (ncol(keep) > c0)
      for (k in (c0 + 1L):ncol(keep))
        keep[, k] <- keep[, k] & keep[, k - 1L] & (I[, k] <= I[, k - 1L])
    if (c0 > 1L)
      for (k in (c0 - 1L):1L)
        keep[, k] <- keep[, k] & keep[, k + 1L] & (I[, k] <= I[, k + 1L])
    keep[, c0] <- TRUE
    
    Iw   <- I * keep
    sumI <- rowSums(Iw)
    mzc  <- rowSums(M * Iw) / sumI
    bad  <- !is.finite(mzc) | sumI <= 0
    mzc[bad] <- mz[cc][bad]                   # fall back to apex grid point
    
    res[[ci]] <- data.table(mz = mzc, intensity = apex_int, apex_area = sumI)
  }
  
  rbindlist(res)
}

sample_peaks <- rbindlist(lapply(names(resampled_data_list), function(s) {
  sp  <- resampled_data_list[[s]]
  out <- centroid_one_sample(as.numeric(sp$mz), as.numeric(sp$intensity),
                             intensity_threshold)
  if (is.null(out) || !nrow(out)) return(NULL)
  out[, sample := s][]
}))

message(sprintf("Detected %s refined apexes across %d samples.",
                format(nrow(sample_peaks), big.mark = ","),
                length(resampled_data_list)))

## ---------------- cluster apexes across samples (no chaining) ----------------
setorder(sample_peaks, mz)
mzv   <- sample_peaks$mz
sep   <- c(TRUE, diff(mzv) > (2 * mzv[-length(mzv)] * ppm_tolerance / 1e6))
sample_peaks[, block := cumsum(sep)]

cluster_block <- function(mzv, intv, samv, ppm) {
  ord      <- order(intv, decreasing = TRUE)
  assigned <- rep(FALSE, length(mzv))
  out      <- vector("list", length(mzv))
  nout     <- 0L
  for (i in ord) {
    if (assigned[i]) next
    tol_i <- mzv[i] * ppm / 1e6
    sel   <- which(!assigned & abs(mzv - mzv[i]) <= tol_i)   # window around the SEED
    assigned[sel] <- TRUE
    wt    <- intv[sel]
    nout  <- nout + 1L
    out[[nout]] <- list(
      mz_centroid = sum(mzv[sel] * wt) / sum(wt),
      intensity   = sum(wt),
      n_peaks     = length(sel),
      n_samples   = length(unique(samv[sel])),
      width_ppm   = if (length(sel) > 1L)
        (max(mzv[sel]) - min(mzv[sel])) / mzv[i] * 1e6 else 0
    )
  }
  rbindlist(lapply(out[seq_len(nout)], as.data.table))
}

centroided_peaks <- sample_peaks[, cluster_block(mz, intensity, sample, ppm_tolerance),
                                 by = block]

centroided_peaks <- centroided_peaks[n_samples >= min_samples_feature]
setorder(centroided_peaks, mz_centroid)
centroided_peaks[, group := .I]

message(sprintf("%d features after clustering | median cluster width %.2f ppm | max %.2f ppm",
                nrow(centroided_peaks),
                median(centroided_peaks$width_ppm),
                max(centroided_peaks$width_ppm)))

## ---------------- non-overlapping extraction windows ----------------
centroided_peaks[, ppm_tol := get_mz_tolerance(mz_centroid, ppm_tolerance)]
centroided_peaks[, mz_min  := mz_centroid - ppm_tol]
centroided_peaks[, mz_max  := mz_centroid + ppm_tol]

if (nrow(centroided_peaks) > 1L) {
  nlast <- nrow(centroided_peaks)
  mid   <- (centroided_peaks$mz_centroid[-nlast] + centroided_peaks$mz_centroid[-1]) / 2
  centroided_peaks[seq_len(nlast - 1L), mz_max := pmin(mz_max, mid)]
  centroided_peaks[2:nlast,             mz_min := pmax(mz_min, mid)]
}

stopifnot(!anyDuplicated(centroided_peaks$mz_centroid))
gc()



# =============================
# STEP 3: BUILD FEATURE MATRIX
# =============================

message("===== STEP 3: Feature Matrix Generation =====")

resampled_data_dt <- rbindlist(lapply(names(resampled_data_list), function(sample_name) {
  spectrum <- as.data.table(resampled_data_list[[sample_name]])
  spectrum[, sample := sample_name]
}), fill = TRUE)

## use the non-overlapping windows built in STEP 2 as-is — do NOT recompute them here
centroided_peaks_dt <- as.data.table(centroided_peaks)[
  , .(mz_centroid, mz_min, mz_max, group)
]

## ---------------- integrate every sample over each feature's window (single join) ----------------
all_combinations <- CJ(
  sample      = unique(resampled_data_dt$sample),
  mz_centroid = centroided_peaks_dt$mz_centroid
)

all_combinations <- merge(
  all_combinations,
  centroided_peaks_dt,
  by = "mz_centroid",
  all.x = TRUE
)

integrated <- resampled_data_dt[
  all_combinations,
  on = .(sample, mz >= mz_min, mz <= mz_max),
  allow.cartesian = TRUE
][, .(
  intensity = sum(intensity, na.rm = TRUE),
  apex      = suppressWarnings(max(intensity, na.rm = TRUE)),
  n_points  = sum(!is.na(intensity))
), by = .(sample, mz_centroid)]

integrated[, detected := is.finite(apex) & apex > intensity_threshold]
integrated[!is.finite(apex), apex := 0]

## ---------------- wide feature matrix: SAME full-window sum for every sample, detected or not ----------------
feature_matrix_df <- dcast(
  integrated[, .(sample, mz_centroid, intensity)],
  sample ~ mz_centroid,
  value.var = "intensity",
  fill = 0
)

rownames(feature_matrix_df) <- feature_matrix_df$sample
feature_matrix_df$sample <- NULL

## ---------------- companion detection-flag matrix, for downstream CV/blank filtering ----------------
detection_matrix_df <- dcast(
  integrated[, .(sample, mz_centroid, detected)],
  sample ~ mz_centroid,
  value.var = "detected",
  fill = FALSE
)
rownames(detection_matrix_df) <- detection_matrix_df$sample
detection_matrix_df$sample <- NULL

message(sprintf("✔ Feature matrix complete: %d samples x %d features.",
                nrow(feature_matrix_df), ncol(feature_matrix_df)))

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
  
  export_dir <- file.path(mzml_export_dir, "accumulated_smoothed_mzML")
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
      fromFile = 1L,
      polarity = export_polarity
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
  
  centroided_export_dir <- file.path(mzml_export_dir, "centroided_mzML")
  dir.create(centroided_export_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Export the per-sample refined apexes generated in STEP 2.
  # Do not use combined_dt or centroided_peaks here:
  # - sample_peaks contains the refined apexes and their sample assignments.
  # - centroided_peaks contains cross-sample feature clusters.
  
  if (!exists("sample_peaks") ||
      is.null(sample_peaks) ||
      nrow(sample_peaks) == 0L) {
    stop(
      "Centroided export failed: sample_peaks is missing or empty. ",
      "Check STEP 2: CENTROIDING."
    )
  }
  
  for (sample_name in unique(sample_peaks$sample)) {
    
    # Select the refined apexes belonging to this sample
    sample_spectrum <- sample_peaks[sample == sample_name]
    
    # Keep valid, positive-intensity centroid peaks
    sample_spectrum <- sample_spectrum[
      is.finite(mz) &
        is.finite(intensity) &
        intensity > 0
    ]
    
    if (nrow(sample_spectrum) == 0L) {
      message("Skipping sample without valid centroid peaks: ", sample_name)
      next
    }
    
    # Sort peaks by m/z
    setorder(sample_spectrum, mz)
    
    mzs  <- as.numeric(sample_spectrum$mz)
    ints <- as.numeric(sample_spectrum$intensity)
    tic_value <- sum(ints, na.rm = TRUE)
    
    # Build one centroided MS1 spectrum per sample
    spectra_data <- S4Vectors::DataFrame(
      mz = I(list(mzs)),
      intensity = I(list(ints)),
      rtime = 10,
      msLevel = 1L,
      centroided = TRUE,
      fromFile = 1L,
      scanWindowLowerLimit = min(mzs),
      scanWindowUpperLimit = max(mzs),
      tic = tic_value,
      scanIndex = 1L,
      polarity = export_polarity
    )
    
    sps <- Spectra::Spectra(
      spectra_data,
      backend = Spectra::MsBackendMemory()
    )
    
    outfile <- file.path(
      centroided_export_dir,
      paste0(sample_name, "_centroided.mzML")
    )
    
    # Allow an existing output file to be regenerated
    if (file.exists(outfile)) {
      file.remove(outfile)
    }
    
    Spectra::export(
      sps,
      Spectra::MsBackendMzR(),
      file = outfile
    )
    
    if (!file.exists(outfile) || file.info(outfile)$size <= 0) {
      stop("Centroided mzML export failed for sample: ", sample_name)
    }
    
    message(
      "✔ Exported centroided spectrum for ", sample_name,
      " (", length(mzs), " peaks; TIC = ", signif(tic_value, 6), ")"
    )
  }
  
  message("✔ Centroided export completed: ", centroided_export_dir)
  
} else {
  message("⛔ Skipping centroid export (EXPORT_SPECTRA = FALSE)")
}




gc()

# =============================
# M/z range limiter to reduce computational load. Molecules greater than 850Da require high computational load due to the increasing number of formulas being calculated.
# =============================

{
  if (!exists("mz_min_limit")) mz_min_limit <- 70
  if (!exists("mz_max_limit")) mz_max_limit <- 850
  
  sample_names <- rownames(feature_matrix_df)
  mz_values <- suppressWarnings(as.numeric(colnames(feature_matrix_df)))
  valid_cols <- which(!is.na(mz_values) & mz_values >= mz_min_limit & mz_values <= mz_max_limit)
  
  feature_matrix_df <- as.data.frame(feature_matrix_df)[, valid_cols, drop = FALSE]
  rownames(feature_matrix_df) <- sample_names
  message("m/z range set to ", mz_min_limit, " - ", mz_max_limit)
}

# =============================
# Isotope clustering based on 13C isotope pattern with a set tolerance of 10ppm.
# =============================

{
  mass_diff <- 1.0033548378
  ppm <- 10
  
  detect_isotopic_groups_with_numbers <- function(neutral_masses, intensities, mass_diff = 1.0033548378, ppm = 10) {
    # Sort masses and intensities together for isotope matching, but restore the
    # original feature/column order before returning results to the caller.
    if (length(neutral_masses) != length(intensities)) {
      stop("neutral_masses and intensities must have the same length")
    }
    n <- length(neutral_masses)
    original_order <- order(neutral_masses, na.last = TRUE)
    sorted_masses <- neutral_masses[original_order]
    sorted_intensities <- intensities[original_order]
    
    labels_sorted <- rep(NA_character_, n)
    group_numbers_sorted <- rep(NA_integer_, n)
    group_id <- 1L
    used_indices <- rep(FALSE, n)
    
    for (i in seq_len(n)) {
      if (used_indices[i]) next
      base_mass <- sorted_masses[i]
      base_intensity <- sorted_intensities[i]
      # Missing/non-finite values cannot be used for isotope matching.
      if (!is.finite(base_mass) || !is.finite(base_intensity)) next
      
      tol <- base_mass * ppm / 1e6
      expected_m1 <- base_mass + mass_diff
      m1_index <- which(is.finite(sorted_masses) &
                          abs(sorted_masses - expected_m1) <= tol & !used_indices)
      m1_index <- m1_index[m1_index != i]
      
      if (length(m1_index) > 0 && is.finite(sorted_intensities[m1_index[1]]) &&
          sorted_intensities[m1_index[1]] < base_intensity) {
        m1 <- m1_index[1]
        labels_sorted[i] <- "monoisotopic"
        labels_sorted[m1] <- "M+1"
        group_numbers_sorted[c(i, m1)] <- group_id
        used_indices[c(i, m1)] <- TRUE
        
        expected_m2 <- base_mass + 2 * mass_diff
        m2_index <- which(is.finite(sorted_masses) &
                            abs(sorted_masses - expected_m2) <= tol & !used_indices)
        m2_index <- m2_index[m2_index != i & m2_index != m1]
        if (length(m2_index) > 0 && is.finite(sorted_intensities[m2_index[1]]) &&
            sorted_intensities[m2_index[1]] < base_intensity &&
            sorted_intensities[m2_index[1]] < sorted_intensities[m1]) {
          m2 <- m2_index[1]
          labels_sorted[m2] <- "M+2"
          group_numbers_sorted[m2] <- group_id
          used_indices[m2] <- TRUE
        }
        
        group_id <- group_id + 1L
        next
      }
      
      labels_sorted[i] <- "monoisotopic"
      group_numbers_sorted[i] <- group_id
      used_indices[i] <- TRUE
      group_id <- group_id + 1L
    }
    
    # original_order maps sorted positions -> original positions.
    labels <- rep(NA_character_, n)
    group_numbers <- rep(NA_integer_, n)
    labels[original_order] <- labels_sorted
    group_numbers[original_order] <- group_numbers_sorted
    return(list(labels = labels, group_numbers = group_numbers))
  }
  message("Isotope groups detected with mass difference set to ", mass_diff, " Da and mass tolerance set to ", ppm, " ppm")
}

# =============================
# Adduct-aware neutral mass calculation. Assumes every peak to be a (de)-protonated species and searches for other possible matching adducts.
# =============================


if (Ionization_method == "Positive") {
  
  message("Calculation of Molecular Forumla using RDISOP (Ionization_method = Positive)")
  
  proton_mass <- 1.007276  # [M+H]+
  na_mass     <- 22.989221 # [M+Na]+
  k_mass      <- 38.963158 # [M+K]+
  ppm <- 10
  
  detect_neutral_masses <- function(mz_values, ppm = 10) {
    adduct_masses <- list(H = proton_mass, Na = na_mass, K = k_mass)
    
    neutral_masses <- rep(NA, length(mz_values))
    adducts <- rep(NA_character_, length(mz_values))
    used <- rep(FALSE, length(mz_values))
    
    for (i in seq_along(mz_values)) {
      if (used[i]) next
      
      candidate_neutral <- mz_values[i] - proton_mass
      
      tol <- mz_values[i] * ppm / 1e6
      
      na_idx <- which(abs(mz_values - (candidate_neutral + na_mass)) <= tol & !used)
      k_idx  <- which(abs(mz_values - (candidate_neutral + k_mass)) <= tol & !used)
      
      neutral_masses[i] <- candidate_neutral
      adducts[i] <- "M+H"
      used[i] <- TRUE
      
      if (length(na_idx) > 0) { neutral_masses[na_idx[1]] <- candidate_neutral; adducts[na_idx[1]] <- "M+Na"; used[na_idx[1]] <- TRUE }
      if (length(k_idx)  > 0) { neutral_masses[k_idx[1]]  <- candidate_neutral; adducts[k_idx[1]]  <- "M+K";  used[k_idx[1]]  <- TRUE }
    }
    
    attr(neutral_masses, "adduct") <- adducts
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
    adducts <- rep(NA_character_, length(mz_values))
    used <- rep(FALSE, length(mz_values))
    
    for (i in seq_along(mz_values)) {
      if (used[i]) next
      
      candidate_neutral <- mz_values[i] + proton_mass
      
      tol <- mz_values[i] * ppm / 1e6
      
      expected_cl <- candidate_neutral + cl_mass
      cl_idx <- which(abs(mz_values - expected_cl) <= tol & !used)
      
      expected_fa <- candidate_neutral + fa_mass
      fa_idx <- which(abs(mz_values - expected_fa) <= tol & !used)
      
      neutral_masses[i] <- candidate_neutral
      adducts[i] <- "M-H"
      used[i] <- TRUE
      
      if (length(cl_idx) > 0) {
        neutral_masses[cl_idx[1]] <- candidate_neutral
        adducts[cl_idx[1]] <- "M+Cl"
        used[cl_idx[1]] <- TRUE
      }
      if (length(fa_idx) > 0) {
        neutral_masses[fa_idx[1]] <- candidate_neutral
        adducts[fa_idx[1]] <- "M+HCOO"
        used[fa_idx[1]] <- TRUE
      }
    }
    
    attr(neutral_masses, "adduct") <- adducts
    return(neutral_masses)
  }
}


{
  mz_values <- as.numeric(colnames(feature_matrix_df))
  intensities <- as.numeric(colMeans(feature_matrix_df, na.rm = TRUE))  # mean over all samples (was: sample 1 only)
  
  neutral_mass_values <- detect_neutral_masses(mz_values, ppm = ppm)
  # Adduct labels from this (single) calculation, kept per feature for the final output matrix
  initial_adducts <- setNames(attr(neutral_mass_values, "adduct"), colnames(feature_matrix_df))
  attr(neutral_mass_values, "adduct") <- NULL
  
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
  
  for (i in 1:nrow(feature_matrix_df_isotope)) {
    if (!(rownames(feature_matrix_df_isotope)[i] %in% c("Isotopic Labels"))) {
      feature_matrix_df_isotope[i, ] <- as.numeric(feature_matrix_df_isotope[i, ])
    }
  }
  
  write.csv(feature_matrix_df_isotope,
            file = out_file("feature_matrix_df_isotope.csv"),
            row.names = TRUE)
}




# =============================
# RDISOP formula calculation with explicit application of the implemented Golden-Rule-inspired filters.
# =============================


{
  
  getEssentialElements <- function(neutral_mass) {
    Da <- abs(neutral_mass)
    
    if (Da <= 140) {
      return(initializeElements(c("C", "H", "N", "O", "P", "S")))  
    } else if (Da <= 180) { 
      return(initializeElements(c("C", "H", "N", "O", "P", "S")))  
    } else if (Da <= 250) { 
      return(initializeElements(c("C", "H", "N", "O", "P", "S")))  
    } else if (Da <= 1000) {
      return(initializeElements(c("C", "H", "N", "O", "P", "S")))  
    }
  }
  
  
  getmassaccuracy <- function (neutral_mass) {
    Da <- abs(neutral_mass)
    
    if (Da <= 300) {
      return(10)  # CHNO
    } else if (Da <= 1000) { 
      return(5)
    }}
  message("RDISOP formula prediction: element-count limits, mass accuracy, isotope-pattern scoring, H/C, heteroatom ratios, element-probability, DBE, valence, and Rdisop validity filters are evaluated; final formula selection requires all hard filters to pass. Mass tolerance: 10 ppm through 300 Da, 5 ppm above 300 Da.")
  
  
  #Golden rule #1: Restriction for element numbers
  
  getMaxElementLimits <- function(neutral_mass) {
    Da <- abs(neutral_mass)  
    
    if (Da <= 100) {
      return("C7H14N5O4P2S2") 
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
    } else if (Da <= 900) { 
      return("C59H106N7O23P3S1")  
    } else if (Da <= 950) { 
      return("C59H102N12O25P5S1")  
    } else if (Da <= 1000) { 
      return("C63H112N10O27P6")  
    } else { 
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
    } else if (Da <= 1000) { 
      return("C35H70")  
    } else { 
      return(NULL) 
    }
  }
  
  
  parseFormulaElements <- function(formula) {
    matches <- regmatches(formula, gregexpr("([A-Z][a-z]?)([0-9]*)", formula))[[1]]
    elements <- gsub("[0-9]", "", matches)
    counts <- as.numeric(gsub("[^0-9]", "", matches))
    counts[is.na(counts)] <- 1
    structure(tapply(counts, elements, sum), names = unique(elements))
  }
  
  
  detectLipidClass <- function(formula) {
    el <- parseFormulaElements(formula)
    
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
  
  neutral_mass_values <- as.numeric(feature_matrix_df_isotope[1, ])
  # Use the isotope groups calculated once above (feature_matrix_df_isotope), no second calculation
  peak_groups <- as.integer(unlist(feature_matrix_df_isotope["Peak Group Numbers", ], use.names = FALSE))
  
  unique_groups <- unique(na.omit(peak_groups))
  
  
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
  mono_index <- indices[which.min(neutral_mass_values[indices])]
  mono_mass <- neutral_mass_values[mono_index]
  
  max_elements <- getMaxElementLimits(mono_mass)
  min_elements <- getMinElementLimits(mono_mass)
  essential_elements <- getEssentialElements(mono_mass)  
  mz_dependent_mass_accuracy <- getmassaccuracy(mono_mass)
  
  
  if (!is.null(max_elements)) {
    masses <- neutral_mass_values[indices]
    intensities <- colMeans(feature_matrix_df[, indices, drop = FALSE], na.rm = TRUE)
    
    # Rules #2-3: pass the mass-accuracy tolerance and observed isotope
    # mass/intensity pattern to Rdisop for formula decomposition.
    molecules <- tryCatch(
      decomposeIsotopes(masses, intensities, ppm = mz_dependent_mass_accuracy, 
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
        required_elements <- c("N", "O", "P", "S")
        for (el in required_elements) {
          if (!(el %in% names(elements))) {
            elements[el] <- 0  
          }
        }
        
        elements <- elements[names(elements) %in% required_elements]
        
        r1 <- na.omit(elements)
        r2 <- r1[r1 > 1]
        
        if (length(r2) <= 2) return(TRUE)
        
        if (length(r2) == 4 && 
            (r2["N"] >= 10 || r2["O"] >= 20 || r2["P"] >= 4 || r2["S"] >= 3)) {
          return(FALSE)
        }
        
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
      
      
      element_probabilities <- sapply(formulas, function(f) {
        el <- parseFormulaElements(f)
        checkElementProbability(el)
      })
      
      
      calculateDBE <- function(formula) {
        el <- parseFormulaElements(formula)
        C <- ifelse("C" %in% names(el), el["C"], 0)
        H <- ifelse("H" %in% names(el), el["H"], 0)
        N <- ifelse("N" %in% names(el), el["N"], 0)
        X <- ifelse("Halogen" %in% names(el), el["Halogen"], 0) 
        P <- ifelse("P" %in% names(el), el["P"], 0)
        S <- ifelse("S" %in% names(el), el["S"], 0)
        
        dbe <- C - H / 2 + N / 2 + P / 2 + 1
        return(dbe)
      }
      dbe_values <- sapply(formulas, calculateDBE)
      
      valid_dbe <- !is.na(dbe_values) & dbe_values >= 0
      
      
      
      # Golden Rule 7: Valence check
      checkValenceCompliance <- function(formula) {
        el <- parseFormulaElements(formula)
        
        valences <- c(
          H = 1, C = 4, N = 3, O = 2, P = 3, S = 2, Na = 1
        )
        
        total_bonds <- 0
        
        for (element in names(el)) {
          if (!element %in% names(valences)) {
            return(FALSE) 
          }
          total_bonds <- total_bonds + valences[element] * el[element]
        }
        
        return(total_bonds %% 2 == 0)
      }
      
      valid_valence <- sapply(formulas, checkValenceCompliance)
      
      lipid_info <- lapply(formulas, detectLipidClass)
      lipid_scores <- sapply(lipid_info, function(x) x$score)
      lipid_classes <- sapply(lipid_info, function(x) x$class)
      
      adjusted_scores <- scores + lipid_scores
      
      
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
        Valid_Rdisop = !is.na(validities) & tolower(trimws(as.character(validities))) == "valid",
        Valid_HC = !is.na(hc_ratios) & is.finite(hc_ratios) &
          hc_ratios >= 0.2 & hc_ratios <= 3.1,
        Valid_Heteroatom_Ratios =
          is.finite(hetero_ratios[, "NtoC"]) &
          is.finite(hetero_ratios[, "OtoC"]) &
          is.finite(hetero_ratios[, "PtoC"]) &
          is.finite(hetero_ratios[, "StoC"]) &
          hetero_ratios[, "NtoC"] < 1.3 &
          hetero_ratios[, "OtoC"] < 1.2 &
          hetero_ratios[, "PtoC"] < 0.3 &
          hetero_ratios[, "StoC"] < 0.8,
        Valid_Element_Probability = !is.na(element_probabilities) & element_probabilities,
        Valid_DBE = !is.na(dbe_values) & is.finite(dbe_values) & dbe_values >= 0,
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
      
      # Apply the formula filters together. These are neutral metabolite
      # formulas, so require carbon before applying the organic H/C and
      # heteroatom rules; an undefined ratio must not pass as NA.
      has_carbon <- vapply(formulas, function(f) {
        el <- parseFormulaElements(f)
        "C" %in% names(el) && is.finite(as.numeric(el["C"])) && as.numeric(el["C"]) > 0
      }, logical(1))
      
      # Golden Rule #4: hydrogen-to-carbon ratio.
      valid_hc <- has_carbon & is.finite(hc_ratios) &
        hc_ratios >= 0.2 & hc_ratios <= 3.1
      
      # Golden Rule #5: heteroatom-to-carbon ratios.
      # Pipeline thresholds: N/C < 1.3; O/C < 1.2; P/C < 0.3; S/C < 0.8.
      valid_hetero <- has_carbon &
        is.finite(hetero_ratios[, "NtoC"]) &
        is.finite(hetero_ratios[, "OtoC"]) &
        is.finite(hetero_ratios[, "PtoC"]) &
        is.finite(hetero_ratios[, "StoC"]) &
        hetero_ratios[, "NtoC"] < 1.3 &
        hetero_ratios[, "OtoC"] < 1.2 &
        hetero_ratios[, "PtoC"] < 0.3 &
        hetero_ratios[, "StoC"] < 0.8
      
      # Golden Rule #6: element-probability heuristic.
      valid_element_probability <- !is.na(element_probabilities) & element_probabilities
      
      # Require Rdisop's own validity flag as an additional hard filter.
      valid_rdisop <- !is.na(validities) & tolower(trimws(as.character(validities))) == "valid"
      
      # DBE and valence are hard filters; missing/unknown results fail closed.
      valid_dbe <- !is.na(dbe_values) & is.finite(dbe_values) & dbe_values >= 0
      valid_valence <- !is.na(valid_valence) & valid_valence
      
      valid_indices <- which(
        valid_rdisop &
          valid_hc &
          valid_hetero &
          valid_element_probability &
          valid_dbe &
          valid_valence
      )
      
      if (length(valid_indices) > 0) {
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

# =============================
# Convert peak group RDISOP information to df for csv export
# =============================


{
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
    molecular_info_row_heteroprob,  
    feature_matrix_df_isotope[-(1:3), ]
  )
  
  
  rownames(feature_matrix_df_final)[4:14] <- c(
    "Molecular_Formula", "Score", "Validity", "DBE", "Valence_Valid", "H/C ratio",
    "N/C ratio", "O/C ratio", "P/C ratio", "S/C ratio", "Heteroatom Probability"  
  )
  
  
  write.csv(feature_matrix_df_final, file = out_file("feature_matrix_with_molecular_info.csv"), row.names = T)
  write.csv(molecular_data, file = out_file("molecular_data.csv"), row.names = F)
}

message("RDISOP calculations complete. Export of feature_matrix_with_molecular_info.csv and molecular_data.csv successful")

message("Performing annotation of features by matching mz values and predicted formulas using the HMDB database and a mass tolerance of ", ppm_tolerance_annotation, " ppm" )


# =============================
# HMDB SUBSET SELECTION (Biospecimen / Origin) based on the SHINY APP settings
# =============================

# hmdb_biospecimen / hmdb_origin are injected by the Shiny app. If they are missing or empty
# (e.g. script run standalone, or an older app version), ALL rows of the CSV would be used,
# so this is reported loudly instead of silently.
if (!exists("hmdb_biospecimen")) hmdb_biospecimen <- NULL
if (!exists("hmdb_origin"))      hmdb_origin <- NULL

message("HMDB selection received from the app -> Biospecimen: ",
        if (length(hmdb_biospecimen)) paste(hmdb_biospecimen, collapse = ", ") else "<NONE>",
        " | Origin: ",
        if (length(hmdb_origin)) paste(hmdb_origin, collapse = ", ") else "<NONE>")
if (length(hmdb_biospecimen) == 0 || length(hmdb_origin) == 0) {
  warning("No HMDB Biospecimen/Origin selection received - using ALL rows of the HMDB CSV. ",
          "Check that the Shiny app passes 'hmdb_biospecimen' and 'hmdb_origin' in the config.",
          call. = FALSE, immediate. = TRUE)
}

select_hmdb_subset <- function(hmdb, biospecimen, origin) {
  hmdb <- as.data.frame(hmdb, stringsAsFactors = FALSE)
  
  missing_cols <- setdiff(c("Biospecimen", "Origin", "HMDB_ID"), names(hmdb))
  if (length(missing_cols)) {
    stop("HMDB CSV is missing required column(s): ", paste(missing_cols, collapse = ", "),
         "\nColumns found: ", paste(names(hmdb), collapse = ", "), call. = FALSE)
  }
  
  # Cells may hold several values (e.g. "Endogenous; Food") -> split into tokens.
  # The selections come from the values found in this CSV (Shiny app), so exact
  # (case-insensitive) token matching is used. NULL/empty selection = use everything.
  split_tokens <- function(x) {
    lapply(strsplit(tolower(trimws(as.character(x))), "[;,|/]"), trimws)
  }
  row_has_any <- function(tokens, wanted) {
    vapply(tokens, function(tk) any(tk %in% tolower(wanted)), logical(1))
  }
  distinct_values <- function(x) {
    v <- trimws(unlist(strsplit(as.character(x), "[;,|/]")))
    v <- v[!is.na(v) & nzchar(v)]
    sort(v[!duplicated(tolower(v))])
  }
  
  if (is.null(biospecimen) || length(biospecimen) == 0) biospecimen <- distinct_values(hmdb$Biospecimen)
  if (is.null(origin) || length(origin) == 0)           origin      <- distinct_values(hmdb$Origin)
  
  keep_bio <- row_has_any(split_tokens(hmdb$Biospecimen), biospecimen)
  
  origin_tokens <- split_tokens(hmdb$Origin)
  origin_hits <- lapply(origin, function(o) row_has_any(origin_tokens, o))
  names(origin_hits) <- origin
  keep_origin <- Reduce(`|`, origin_hits)
  
  keep <- keep_bio & keep_origin
  
  # ---- Report ----
  message("HMDB subset selection for annotation:")
  message("  Biospecimen values in CSV: ",
          paste0(names(table(hmdb$Biospecimen)), " (", table(hmdb$Biospecimen), ")", collapse = ", "))
  message("  Origin values in CSV: ",
          paste0(names(table(hmdb$Origin)), " (", table(hmdb$Origin), ")", collapse = ", "))
  message("  Selected biospecimen: ", paste(biospecimen, collapse = ", "),
          " -> ", sum(keep_bio), " rows")
  for (o in origin) message("  Selected origin '", o, "': ", sum(origin_hits[[o]]), " rows")
  
  hmdb_sub <- hmdb[keep, , drop = FALSE]
  
  # The same metabolite can occur in several rows; keep each HMDB_ID once
  n_before <- nrow(hmdb_sub)
  hmdb_sub <- hmdb_sub[!duplicated(hmdb_sub$HMDB_ID), , drop = FALSE]
  message("  -> ", nrow(hmdb_sub), " unique HMDB entries retained for annotation (of ",
          nrow(hmdb), " rows in CSV; ", n_before - nrow(hmdb_sub), " duplicate HMDB_IDs removed)")
  
  if (nrow(hmdb_sub) == 0) {
    stop("No HMDB entries left after Biospecimen/Origin filtering. ",
         "Check the selection and the spelling of the values in the CSV.", call. = FALSE)
  }
  
  data.table::as.data.table(hmdb_sub)
}

HMDB_Metabolites_full <- HMDB_Metabolites
HMDB_Metabolites <- select_hmdb_subset(HMDB_Metabolites_full, hmdb_biospecimen, hmdb_origin)


# =============================
# ANNOTATION BASED ON MOLECULAR FORMULA using the provided HMDB CSV defined in the SHINY APP
# =============================


feature_matrix_df_final_annotated <- feature_matrix_df_final

library(dplyr)


feature_matrix_df_final_annotated["Molecular_Formula", ] <- tolower(trimws(as.character(unlist(feature_matrix_df_final_annotated["Molecular_Formula", ]))))
HMDB_Metabolites$CHEMICAL_FORMULA <- trimws(tolower(HMDB_Metabolites$CHEMICAL_FORMULA))


feature_matrix_df_final_annotated["Formula_Matched_NAME", ] <- sapply(feature_matrix_df_final_annotated["Molecular_Formula", ], function(MF) {
  if (is.na(MF) || MF == "") return(NA)
  match_idx <- which(HMDB_Metabolites$CHEMICAL_FORMULA == MF)
  if (length(match_idx) > 0) paste(HMDB_Metabolites$NAME[match_idx], collapse = ", ") else NA
})


feature_matrix_df_final_annotated["Formula_Matched_HMDB", ] <- sapply(feature_matrix_df_final_annotated["Molecular_Formula", ], function(MF) {
  if (is.na(MF) || MF == "") return(NA)
  match_idx <- which(HMDB_Metabolites$CHEMICAL_FORMULA == MF)
  if (length(match_idx) > 0) paste(HMDB_Metabolites$HMDB_ID[match_idx], collapse = ", ") else NA
})


feature_matrix_df_final_annotated["Formula_Matched_SMILES", ] <- sapply(feature_matrix_df_final_annotated["Molecular_Formula", ], function(MF) {
  if (is.na(MF) || MF == "") return(NA)
  match_idx <- which(HMDB_Metabolites$CHEMICAL_FORMULA == MF)
  if (length(match_idx) > 0) paste(HMDB_Metabolites$SMILES[match_idx], collapse = ", ") else NA
})


list_columns <- sapply(feature_matrix_df_final_annotated, is.list)
if (any(list_columns)) {
  feature_matrix_df_final_annotated[list_columns] <- lapply(feature_matrix_df_final_annotated[list_columns], function(x) {
    sapply(x, function(y) if (is.list(y)) paste(unlist(y), collapse = ", ") else y)
  })
}

write.csv(feature_matrix_df_final_annotated,
          file = out_file("feature_matrix_with_formula_annotations.csv"),
          row.names = TRUE)

# =============================
# ANNOTATION BASED ON m/z based on used ionization mode using the provided HMDB CSV defined in the SHINY APP
# =============================


mz_features <- as.numeric(colnames(feature_matrix_df_final_annotated))

calculate_ppm_error <- function(mz_observed, mz_theoretical) {
  abs((mz_observed - mz_theoretical) / mz_theoretical) * 1e6
}

# =============================
# POSITIVE ION MODE
# =============================

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
  
  # =============================
  # NEGATIVE ION MODE
  # =============================
  
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

# =============================
# RUN MATCHING ON ALL FEATURES
# =============================


mz_names        <- vector("character", length(mz_features))
mz_hmdb_ids     <- vector("character", length(mz_features))
mz_adduct_types <- vector("character", length(mz_features))

for (i in seq_along(mz_features)) {
  
  mz_value <- mz_features[i]
  matches  <- match_mz(mz_value, HMDB_Metabolites, ppm_tolerance_annotation)
  
  if (nrow(matches) > 0) {
    
    matches <- matches %>%
      rowwise() %>%
      mutate(
        min_ppm = min(c_across(starts_with("ppm_error_")), na.rm = TRUE)
      ) %>%
      ungroup() %>%
      arrange(min_ppm)
    
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

# ---- Check: every m/z match must come from the filtered HMDB subset ----
{
  matched_ids <- unique(trimws(unlist(strsplit(na.omit(mz_hmdb_ids), ","))))
  matched_ids <- matched_ids[nzchar(matched_ids)]
  outside <- setdiff(matched_ids, HMDB_Metabolites$HMDB_ID)
  if (length(outside)) {
    warning(length(outside), " m/z-matched HMDB_IDs are not in the filtered subset!", call. = FALSE, immediate. = TRUE)
  }
  origin_of_matches <- HMDB_Metabolites$Origin[HMDB_Metabolites$HMDB_ID %in% matched_ids]
  message("m/z annotation used ", nrow(HMDB_Metabolites), " HMDB entries; ",
          length(matched_ids), " distinct HMDB_IDs matched. Origin of matches: ",
          if (length(origin_of_matches)) paste0(names(table(origin_of_matches)), " (", table(origin_of_matches), ")", collapse = ", ") else "none")
}

feature_matrix_df_final_annotated["mz_matched_NAME", ]   <- mz_names
feature_matrix_df_final_annotated["mz_matched_HMDB", ]   <- mz_hmdb_ids
feature_matrix_df_final_annotated["mz_matched_ADDUCT", ] <- mz_adduct_types

# =============================
# PREDICTED ADDUCT LABELING to return the adduct prediction into the dataframe
# =============================


# Re-use the adduct labels from the first calculation (detect_neutral_masses), looked up by feature
predicted_adducts <- unname(initial_adducts[colnames(feature_matrix_df_final_annotated)])

adduct_labels_row <- as.data.frame(t(predicted_adducts), stringsAsFactors = FALSE)
colnames(adduct_labels_row) <- colnames(feature_matrix_df_final_annotated)
rownames(adduct_labels_row) <- "Predicted_Adduct"

# =============================
# ASSEMBLE FINAL ANNOTATED MATRIX
# =============================


annotations_df <- feature_matrix_df_final_annotated[
  c("Formula_Matched_NAME",
    "Formula_Matched_HMDB",
    "Formula_Matched_SMILES",
    "mz_matched_NAME",
    "mz_matched_HMDB",
    "mz_matched_ADDUCT"), , drop = FALSE]

feature_matrix_df_final_annotated <- feature_matrix_df_final_annotated[
  !rownames(feature_matrix_df_final_annotated) %in% rownames(annotations_df), , drop = FALSE]

feature_matrix_df_final_annotated["Molecular_Formula", ] <-
  toupper(trimws(as.character(unlist(feature_matrix_df_final_annotated["Molecular_Formula", ]))))

feature_matrix_df_final_annotated <- rbind(
  feature_matrix_df_final_annotated[1:2, , drop = FALSE],
  adduct_labels_row,
  feature_matrix_df_final_annotated[-(1:2), , drop = FALSE]
)

feature_matrix_df_final_annotated <- rbind(
  feature_matrix_df_final_annotated[1:15, ],  
  annotations_df,
  feature_matrix_df_final_annotated[-(1:15), ]
)

list_columns <- sapply(feature_matrix_df_final_annotated, is.list)
if (any(list_columns)) {
  feature_matrix_df_final_annotated[list_columns] <-
    lapply(feature_matrix_df_final_annotated[list_columns], function(x) {
      sapply(x, function(y) if (is.list(y)) paste(unlist(y), collapse = ", ") else y)
    })
}

cat("Predicted_Adduct row added at row 3.\n")


# =============================
# EXPORT FULL MATRIX of ANNOTATED AND NON-ANNOTATED FEATURES
# =============================


write.csv(feature_matrix_df_final_annotated,
          file = out_file("feature_matrix_with_formula_and_mz_annotations.csv"),
          row.names = TRUE)

# =============================
# EXPORT FILTERED MATRIX (ANNOTATED FEATURES ONLY)
# =============================


na_cols <- which(
  is.na(feature_matrix_df_final_annotated["Formula_Matched_NAME", ]) &
    is.na(feature_matrix_df_final_annotated["mz_matched_NAME", ])
)

# Guard: -integer(0) would drop ALL columns if every feature is annotated
feature_matrix_df_final_annotated_filtered <- if (length(na_cols) > 0) {
  feature_matrix_df_final_annotated[, -na_cols, drop = FALSE]
} else {
  feature_matrix_df_final_annotated
}

write.csv(feature_matrix_df_final_annotated_filtered,
          file = out_file("feature_matrix_with_formula_and_mz_annotations_filtered.csv"),
          row.names = TRUE)



# =============================
# Parse annotations to replace m/z column names for METABOANALYST Styled output.
# =============================



library(stringr)


split_top_level <- function(x, require_space = TRUE) {
  chars  <- strsplit(x, "")[[1]]
  n      <- length(chars)
  depth  <- 0L
  result <- character(0)
  current <- ""
  
  for (k in seq_len(n)) {
    ch <- chars[k]
    if (ch == "(") {
      depth <- depth + 1L; current <- paste0(current, ch)
    } else if (ch == ")") {
      depth <- depth - 1L; current <- paste0(current, ch)
    } else if (ch %in% c(",", ";") && depth == 0L &&
               (!require_space || (k < n && grepl("\\s", chars[k + 1L])))) {
      result <- c(result, trimws(current)); current <- ""
    } else {
      current <- paste0(current, ch)
    }
  }
  if (nchar(trimws(current)) > 0) result <- c(result, trimws(current))
  result
}


LIPID_CLASS_REGEX <- paste0(
  "^(TG|DG|MG|",
  "PC|PE|PS|PI|PG|PA|",
  "LysoPC|LysoPE|LysoPI|LysoPG|LysoPS|LysoPA|",
  "SM|Cer|GlcCer|GalCer|LacCer|HexCer|",
  "Trihexosylceramide|Tetrahexosylceramide|",
  "Ganglioside(\\s+\\S+)?|",
  "CE|FAHFA|WE",
  ")\\s*\\("
)

is_lipid_name <- function(name) {
  grepl(LIPID_CLASS_REGEX, name) && grepl("\\d+:\\d+", name)
}

SPHINGO_CLASSES <- c("SM", "Cer", "GlcCer", "GalCer", "LacCer", "HexCer",
                     "Trihexosylceramide", "Tetrahexosylceramide", "Ganglioside")


convert_to_sum_composition <- function(lipid_name) {
  if (is.na(lipid_name) || lipid_name == "") return(lipid_name)
  
  name_norm <- gsub("\\s+\\(", "(", trimws(lipid_name))
  
  name_norm <- gsub("^(Ganglioside)\\s*\\S+\\s*\\(", "\\1(", name_norm)
  
  m <- regmatches(name_norm, regexec("^([A-Za-z]+)\\((.+)\\)$", name_norm))[[1]]
  if (length(m) != 3) return(lipid_name)
  
  lipid_class   <- m[2]
  chain_content <- m[3]
  
  chain_clean <- gsub("\\([^)]*\\)", "", chain_content)
  
  chains <- trimws(unlist(strsplit(chain_clean, "[/_]")))
  
  total_c      <- 0L
  total_db     <- 0L
  ether_prefix <- ""  
  
  sphingo_prefix <- ""
  if (lipid_class %in% SPHINGO_CLASSES) {
    if (grepl("^d", chain_clean)) {
      sphingo_prefix <- "d"
    } else if (grepl("^t", chain_clean)) {
      sphingo_prefix <- "t"
    }
  }
  
  for (chain in chains) {
    chain <- trimws(chain)
    if (chain == "" || chain == "0:0" || chain == "0") next
    
    chain <- gsub("^[dt]", "", chain)
    
    if (grepl("^O-", chain)) {
      ether_prefix <- "-O"
      chain        <- sub("^O-", "", chain)
    } else if (grepl("^P-", chain)) {
      ether_prefix <- "-P"
      chain        <- sub("^P-", "", chain)
    }
    
    parts <- strsplit(chain, ":")[[1]]
    if (length(parts) == 2) {
      c_val  <- suppressWarnings(floor(as.numeric(trimws(parts[1]))))
      db_val <- suppressWarnings(floor(as.numeric(trimws(parts[2]))))
      if (!is.na(c_val))  total_c  <- total_c  + c_val
      if (!is.na(db_val)) total_db <- total_db + db_val
    }
  }
  
  if (total_c == 0L) return(lipid_name)
  
  paste0(lipid_class, ether_prefix, " ", sphingo_prefix, total_c, ":", total_db)
}


split_annotation_entries <- function(x) {
  x <- trimws(as.character(x))
  
  # 1) unambiguous delimiter (see note below)
  if (grepl(" | ", x, fixed = TRUE))
    return(trimws(strsplit(x, " | ", fixed = TRUE)[[1]]))
  
  # 2) ppm-annotated cells: anchor on the "(... ppm ...)" blocks
  m <- gregexpr("\\([^()]*ppm[^()]*\\)", x, perl = TRUE)[[1]]
  if (m[1] != -1) {
    ends    <- as.integer(m) + attr(m, "match.length") - 1L
    starts  <- c(1L, head(ends, -1L) + 1L)
    entries <- substring(x, starts, ends)
    return(trimws(sub("^[,;]\\s*", "", trimws(entries))))
  }
  
  # 3) last resort
  split_top_level(x)
}

parse_first_annotation_and_adduct <- function(name_cell, adduct_cell = NA) {
  if (is.na(name_cell) || name_cell %in% c("", "NA")) return(list(name = NA, adduct = NA))
  
  entries     <- split_annotation_entries(name_cell)
  first_entry <- if (length(entries)) entries[1] else as.character(name_cell)
  
  name_first <- gsub("\\s*\\([^)]*ppm[^)]*\\)", "", first_entry, ignore.case = TRUE)
  name_first <- trimws(gsub("\\s*\\d+\\.?\\d*\\s*ppm", "", name_first, ignore.case = TRUE))
  
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
  
  if (!is.na(adduct_cell) && adduct_cell != "" && adduct_cell != "NA")
    adduct_first <- trimws(unlist(strsplit(as.character(adduct_cell), "[,;]"))[1])
  
  list(name = name_first, adduct = adduct_first)
}


mz_values <- as.numeric(colnames(feature_matrix_df_final_annotated_filtered))
final_annotation_names <- character(length(mz_values))

for (i in seq_along(mz_values)) {
  
  mz_val           <- round(mz_values[i], 4)
  annotation_found <- FALSE
  
  # Priority 1: Formula_Matched_NAME
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

# =============================
# MetaboAnalyst CSV with annotated column names
# =============================


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
  file = out_file("feature_matrix_metabo.csv"),
  row.names = FALSE
)

cat("MetaboAnalyst-ready file exported with annotated column names\n")

# =============================
# Transposed CSV with annotated metabolite row names and class, index and batch annotation.
# =============================



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
  file = out_file("feature_matrix_transposed_with_metadata.csv"),
  row.names = FALSE
)

cat("Transposed file exported with annotated metabolite names\n")
message("Annotation feature matrices saved. Processing and annotation completed.")

# =============================
# REPRODUCIBILITY: close run_parameters.yml
# =============================
update_run_parameters(list(
  run_status       = "completed",
  run_finished     = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
  run_duration_min = round(as.numeric(difftime(Sys.time(), run_started_time, units = "mins")), 1),
  output_files     = as.list(sort(list.files(output_dir, recursive = TRUE))),
  session_info     = get_session_info()
))
message("Run parameters saved: ", file.path(output_dir, "run_parameters.yml"))