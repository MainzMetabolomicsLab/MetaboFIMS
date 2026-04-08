# app.R (enhanced with reset functionality)
packages <- list(
  bioc = c("xcms", "MSnbase", "CAMERA", "Rdisop", "BiocParallel", "Spectra", "S4Vectors", "mzR"),
  cran = c("R.utils", "webchem", "data.table", "dplyr", "tidyverse", "ggplot2",
           "stringdist", "shiny", "shinyFiles", "DT", "callr")
)

install_if_missing <- function(pkg, type = c("cran", "bioc")) {
  type <- match.arg(type)
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, " from ", type)
    if (type == "cran") {
      install.packages(pkg)
    } else if (type == "bioc") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }
  library(pkg, character.only = TRUE)
}

invisible(lapply(packages$bioc, install_if_missing, type = "bioc"))
invisible(lapply(packages$cran, install_if_missing, type = "cran"))

message("✅ All packages installed and loaded.")

sys_info <- list(
  R_version = R.version.string,
  OS_type = .Platform$OS.type,
  OS_name = Sys.info()[["sysname"]],
  release = Sys.info()[["release"]],
  version = Sys.info()[["version"]],
  machine = Sys.info()[["machine"]],
  user = Sys.info()[["user"]],
  n_cores = parallel::detectCores(logical = TRUE),
  total_mem_GB = round(ps::ps_system_memory()$total / 1024^3, 2)
)

cat("===== System Info =====\n")
for (nm in names(sys_info)) cat(sprintf("%-20s: %s\n", nm, sys_info[[nm]]))
cat("=======================\n\n")
library(shiny)
library(shinyFiles)
library(DT)
library(callr)
# ensure YAML for config save/load
ensure_packages <- function(pkgs) {
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
}
ensure_packages(c("yaml"))
library(yaml)

# ========== OS detection helpers ==========
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

# make logs folder
if (!dir.exists("logs")) dir.create("logs", recursive = TRUE)

# ========== UI ==========
ui <- fluidPage(
  titlePanel("FIMS Data processing of positive or negative ion mode data"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Ionization mode"),
      radioButtons("Ionization_method", "Select ionization mode",
                   choices = c("Positive", "Negative"), selected = "Positive", inline = TRUE),
      hr(),
      h4("Script & folders"),
      fluidRow(
        column(9, textInput("script_path", "R Script",
                            value = default_path(
                              win = file.path("C:/Users/fabia/Seafile/Meine Bibliothek/FIA",
                                              "FIA_XCMS_isotope_neutral_golden_rulez_incl_smoothing_only_RT_pos_neg_combined.R"),
                              mac = file.path("/Users/Fabian Schmitt/Seafile/Meine Bibliothek/FIA",
                                              "FIA_XCMS_isotope_neutral_golden_rulez_incl_smoothing_only_RT_pos_neg_combined.R"),
                              linux = file.path("/home", "User", "etc.",
                                                "FIA_XCMS_isotope_neutral_golden_rulez_incl_smoothing_only_RT_pos_neg_combined.R")
                            ))),
        column(3, shinyFilesButton("browse_script", "Browse", "Choose R script", multiple = FALSE))
      ),
      fluidRow(
        column(9, textInput("data_dir", "Processing data directory (data_dir)",
                            value = default_path(
                              win = file.path("D:", "Arbeit", "FIA", "processing"),
                              mac = file.path("/Volumes/T7/Arbeit/FIA/processing"),
                              linux = file.path("~/processing")
                            ))),
        column(3, shinyDirButton("browse_data_dir", "Browse", "Choose folder"))
      ),
      fluidRow(
        column(9, textInput("raw_dir", "Raw files directory (data_dir_raw_files)",
                            value = default_path(
                              win = file.path("D:", "Arbeit", "FIA", "processing", "raw_mzML"),
                              mac = file.path("/Volumes/T7/Arbeit/FIA/processing/raw_mzML"),
                              linux = file.path("~/processing/raw_mzML")
                            ))),
        column(3, shinyDirButton("browse_raw_dir", "Browse", "Choose folder"))
      ),
      fluidRow(
        column(9, textInput("file_path", "HMDB annotation CSV (file_path)",
                            value = default_path(
                              win = file.path("D:", "Arbeit", "FIA", "HMDB_BLOOD_ENDO_FOOD_DRUG_DRUGMET.csv"),
                              mac = file.path("/Volumes/T7/Arbeit/FIA/HMDB_BLOOD_ENDO_FOOD_DRUG_DRUGMET.csv"),
                              linux = file.path("~/HMDB_BLOOD.csv")
                            ))),
        column(3, shinyFilesButton("browse_file_path", "Browse", "Choose CSV", multiple = FALSE))
      ),
      hr(),
      h4("RT integration window (seconds) for 1D spectra accumulation"),
      numericInput("rt_min", "rt_min (s)", value = 5, min = 0, step = 1),
      numericInput("rt_max", "rt_max (s)", value = 25, min = 0, step = 1),
      hr(),
      h4("Export of Spectra"),
      checkboxInput("EXPORT_SPECTRA", "Export 1D spectra as .mzML", value = FALSE),
      hr(),
      h4("Smoothing and Centroiding parameters"),
      checkboxInput("apply_smoothing", "Apply smoothing", value = TRUE),
      numericInput("smoothing_n", "Smoothing window (n)", value = 7, min = 1, step = 2),
      numericInput("smoothing_p", "Smoothing polynomial (p)", value = 2, min = 0, step = 1),
      numericInput("intensity_threshold", "Intensity threshold", value = 5000, min = 0, step = 1),
      numericInput("ppm_tolerance", "PPM tolerance", value = 20, min = 0, step = 1),
      hr(),
      h4("Annotation parameters"),
      numericInput("ppm_tolerance_annotation", "PPM tolerance (annotation)", value = 20, min = 0, step = 1),
      hr(),
      h4("Run Controls & Config"),
      fluidRow(
        column(6, actionButton("run_btn", "Run pipeline", class = "btn-primary")),
        column(6, actionButton("stop_btn", "Stop (graceful)", class = "btn-danger"))
      ),
      br(),
      fluidRow(
        column(8, textInput("profile_name", "Profile name - Save to working directory (getwd())", value = "")),
        column(4, actionButton("save_config", "Save config"))
      ),
      fileInput("load_config_file", "Load config (.yml)", accept = c(".yml", ".yaml")),
      br(),
      downloadButton("download_log", "Download log")
    ),
    mainPanel(
      width = 8,
      h4("System Info"),
      tags$div(
        style = "height:180px; overflow-y:scroll; border:1px solid #ccc; padding:5px; font-family:monospace; background:#f0f0f0;",
        verbatimTextOutput("system_info_ui")
      ),
      hr(),
      h4("Live Console / Log"),
      tags$div(
        style = "height:300px; overflow-y:scroll; border:1px solid #ccc; padding:5px; font-family:monospace; background:#f8f8f8;",
        verbatimTextOutput("console_out")
      ),
      hr(),
      h4("Environment / Run Info"),
      verbatimTextOutput("env_info"),
      hr(),
      h4("Data Directory Explorer"),
      fluidRow(
        column(9, textInput("csv_explorer_dir", "Directory to explore", value = "")),
        column(3, shinyDirButton("browse_csv_explorer", "Browse", "Choose folder"))
      ),
      fluidRow(
        column(6, actionButton("refresh_files", "Refresh file list")),
        column(6, actionButton("open_output_folder", "Open output folder"))
      ),
      br(),
      DTOutput("file_browser"),
      hr(),
      h4("Preview selected CSV"),
      DTOutput("file_preview"),
      hr(),
      h4("Progress"),
      uiOutput("progress_ui"),
      hr(),
      h5("Notes:"),
      p("App sources the provided R script in a clean environment. Ensure the script respects the variables passed in the config.")
    )
  )
)

# ========== SERVER ==========
server <- function(input, output, session) {
  # ---------- roots for shinyFiles ----------
  if (OS_TYPE == "windows") {
    roots <- c(home = normalizePath("~"), "C:" = "C:/", "D:" = "D:/", "E:" = "E:/", "F:" = "F:/")
  } else if (OS_TYPE == "macos") {
    roots <- c(home = normalizePath("~"), root = "/", Volumes = "/Volumes")
  } else {
    roots <- c(home = normalizePath("~"), root = "/")
  }
  
  shinyFileChoose(input, "browse_script", roots = roots, session = session)
  shinyFileChoose(input, "browse_file_path", roots = roots, session = session)
  shinyFileChoose(input, "browse_DDA_path", roots = roots, session = session)
  shinyDirChoose(input, "browse_data_dir", roots = roots, session = session)
  shinyDirChoose(input, "browse_raw_dir", roots = roots, session = session)
  shinyDirChoose(input, "browse_csv_explorer", roots = roots, session = session)
  
  output$system_info_ui <- renderPrint({
    sys_info <- list(
      R_version = R.version.string,
      OS_type = .Platform$OS.type,
      OS_name = Sys.info()[["sysname"]],
      release = Sys.info()[["release"]],
      version = Sys.info()[["version"]],
      machine = Sys.info()[["machine"]],
      user = Sys.info()[["user"]],
      n_cores = parallel::detectCores(logical = TRUE),
      total_mem_GB = round(ps::ps_system_memory()$total / 1024^3, 2)
    )
    cat("===== SYSTEM INFO =====\n")
    for (nm in names(sys_info)) {
      cat(sprintf("%-15s: %s\n", nm, sys_info[[nm]]))
    }
    cat("=======================\n")
  })
  
  # ---------- reactive values ----------
  rv <- reactiveValues(
    log = character(),
    running = FALSE,
    r_process = NULL,
    file_list = NULL,
    preview_data = NULL,
    progress = 0,
    progress_total = NA,
    est_seconds_per_mzml = 15  # crude default estimate
  )
  
  append_log <- function(msg) {
    line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", msg)
    rv$log <- c(rv$log, line)
    if (length(rv$log) > 10000) rv$log <- tail(rv$log, 10000)
    # write to rotating daily log
    logfn <- file.path("logs", paste0(Sys.Date(), "_fia_run_log.txt"))
    cat(line, file = logfn, sep = "\n", append = TRUE)
  }
  
  output$console_out <- renderText({ paste(tail(rv$log, 200), collapse = "\n") })
  
  # ---------- parse shinyFiles selections ----------
  observeEvent(input$browse_script, {
    fileinfo <- parseFilePaths(roots, input$browse_script)
    if (nrow(fileinfo) > 0) updateTextInput(session, "script_path", value = as.character(fileinfo$datapath))
  })
  observeEvent(input$browse_file_path, {
    fileinfo <- parseFilePaths(roots, input$browse_file_path)
    if (nrow(fileinfo) > 0) updateTextInput(session, "file_path", value = as.character(fileinfo$datapath))
  })
  observeEvent(input$browse_DDA_path, {
    fileinfo <- parseFilePaths(roots, input$browse_DDA_path)
    if (nrow(fileinfo) > 0) updateTextInput(session, "data_dir_DDA", value = as.character(fileinfo$datapath))
  })
  observeEvent(input$browse_data_dir, {
    dirpath <- parseDirPath(roots, input$browse_data_dir)
    if (length(dirpath) > 0) updateTextInput(session, "data_dir", value = dirpath)
  })
  observeEvent(input$browse_raw_dir, {
    dirpath <- parseDirPath(roots, input$browse_raw_dir)
    if (length(dirpath) > 0) updateTextInput(session, "raw_dir", value = dirpath)
  })
  observeEvent(input$browse_csv_explorer, {
    dirpath <- parseDirPath(roots, input$browse_csv_explorer)
    if (length(dirpath) > 0) {
      updateTextInput(session, "csv_explorer_dir", value = dirpath)
      refresh_file_list()  # immediate refresh
    }
  })
  
  # ---------- auto-detect polarity (simple heuristic) ----------
  detect_polarity <- function(raw_dir) {
    if (is.null(raw_dir) || !nzchar(raw_dir) || !dir.exists(raw_dir)) return(NULL)
    files <- list.files(raw_dir, recursive = TRUE, full.names = FALSE)
    files_lc <- tolower(files)
    pos_hits <- sum(grepl("pos", files_lc) | grepl("_p_", files_lc) | grepl("_pos", files_lc))
    neg_hits <- sum(grepl("neg", files_lc) | grepl("_n_", files_lc) | grepl("_neg", files_lc))
    if (pos_hits > neg_hits && pos_hits > 0) return("Positive")
    if (neg_hits > pos_hits && neg_hits > 0) return("Negative")
    return(NULL)
  }
  
  observeEvent(input$raw_dir, {
    pol <- detect_polarity(input$raw_dir)
    if (!is.null(pol)) {
      # only update if user didn't manually override recently
      # do not overwrite user's explicit selection unless it's different from detected
      isolate({
        if (is.null(input$.last_manual_ion) || input$.last_manual_ion != TRUE) {
          updateRadioButtons(session, "Ionization_method", selected = pol)
        }
      })
    }
  })
  
  # track manual ionization change so auto-detect doesn't overwrite user choice
  observeEvent(input$Ionization_method, {
    session$sendCustomMessage("ion_manual", list(val = TRUE))
    # we'll set an input flag via JS below (but we can store in reactive)
    # store a small flag in session$userData
    session$userData$ion_manual <- TRUE
  }, ignoreInit = TRUE)
  
  # ---------- Save / Load config ----------
  make_config_list <- function() {
    list(
      Ionization_method = input$Ionization_method,
      data_dir = input$data_dir,
      data_dir_raw_files = input$raw_dir,
      file_path = input$file_path,
      data_dir_DDA = input$data_dir_DDA,
      EXPORT_SPECTRA = input$EXPORT_SPECTRA,
      apply_smoothing = input$apply_smoothing,
      smoothing_n = input$smoothing_n,
      smoothing_p = input$smoothing_p,
      intensity_threshold = input$intensity_threshold,
      ppm_tolerance = input$ppm_tolerance,
      ppm_tolerance_annotation = input$ppm_tolerance_annotation,
      rt_min = input$rt_min,
      rt_max = input$rt_max
    )
  }
  
  observeEvent(input$save_config, {
    cfg <- make_config_list()
    fname <- input$profile_name
    if (!nzchar(fname)) fname <- paste0("config_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    fname <- paste0(fname, ".yml")
    yaml::write_yaml(cfg, fname)
    showNotification(paste("Saved config to", fname), type = "message")
    append_log(paste("Saved config to", fname))
  })
  
  observeEvent(input$load_config_file, {
    f <- input$load_config_file
    req(f)
    tryCatch({
      cfg <- yaml::read_yaml(f$datapath)
      # populate UI
      updateRadioButtons(session, "Ionization_method", selected = cfg$Ionization_method %||% input$Ionization_method)
      updateTextInput(session, "data_dir", value = cfg$data_dir %||% input$data_dir)
      updateTextInput(session, "raw_dir", value = cfg$data_dir_raw_files %||% input$raw_dir)
      updateTextInput(session, "file_path", value = cfg$file_path %||% input$file_path)
      updateTextInput(session, "data_dir_DDA", value = cfg$data_dir_DDA %||% input$data_dir_DDA)
      updateCheckboxInput(session, "EXPORT_SPECTRA", value = cfg$EXPORT_SPECTRA %||% input$EXPORT_SPECTRA)
      updateCheckboxInput(session, "apply_smoothing", value = cfg$apply_smoothing %||% input$apply_smoothing)
      updateNumericInput(session, "smoothing_n", value = cfg$smoothing_n %||% input$smoothing_n)
      updateNumericInput(session, "smoothing_p", value = cfg$smoothing_p %||% input$smoothing_p)
      updateNumericInput(session, "intensity_threshold", value = cfg$intensity_threshold %||% input$intensity_threshold)
      updateNumericInput(session, "ppm_tolerance", value = cfg$ppm_tolerance %||% input$ppm_tolerance)
      updateNumericInput(session, "ppm_tolerance_annotation", value = cfg$ppm_tolerance_annotation %||% input$ppm_tolerance_annotation)
      updateNumericInput(session, "rt_min", value = cfg$rt_min %||% input$rt_min)
      updateNumericInput(session, "rt_max", value = cfg$rt_max %||% input$rt_max)
      showNotification("Config loaded", type = "message")
      append_log(paste("Loaded config from", f$name))
    }, error = function(e) {
      showNotification(paste("Failed to load config:", e$message), type = "error")
    })
  })
  
  # ---------- Environment / Info ----------
  output$env_info <- renderPrint({
    script_exists <- file.exists(input$script_path)
    mzml_count <- if (dir.exists(input$raw_dir)) length(list.files(input$raw_dir, pattern = "\\.(mzXML|mzML)$", ignore.case = TRUE)) else 0
    script_md5 <- tryCatch(digest::digest(file = input$script_path, algo = "md5"), error = function(e) NA)
    list(
      R_version = R.version.string,
      script_exists = script_exists,
      script_md5 = script_md5,
      mzml_count = mzml_count,
      Ionization_method = input$Ionization_method,
      pipeline_status = if(rv$running) "RUNNING" else "IDLE"
    )
  })
  
  # Switch directories based on ionization mode
  observeEvent(input$Ionization_method, {
    
    if (input$Ionization_method == "Positive") {
      
      updateTextInput(session, "data_dir",
                      value = default_path(
                        win = file.path("D:", "Arbeit", "FIA", "processing"),
                        mac = file.path("/Volumes/T7/Arbeit/FIA/processing"),
                        linux = file.path("~/processing")
                      ))
      
      updateTextInput(session, "raw_dir",
                      value = default_path(
                        win = file.path("D:", "Arbeit", "FIA", "processing", "raw_mzML"),
                        mac = file.path("/Volumes/T7/Arbeit/FIA/processing/raw_mzML"),
                        linux = file.path("~/processing/raw_mzML")
                      ))
      
      updateTextInput(session, "data_dir_DDA",
                      value = default_path(
                        win = file.path("D:", "Arbeit", "FIA", "DDA", "processing", "filtered_feature_df_annotated.csv"),
                        mac = file.path("/Volumes/T7/Arbeit/FIA/DDA/processing/filtered_feature_df_annotated.csv"),
                        linux = file.path("~/filtered_feature_df_annotated.csv")
                      ))
      
    } else {
      
      updateTextInput(session, "data_dir",
                      value = default_path(
                        win = file.path("D:", "Arbeit", "FIA", "processing_neg"),
                        mac = file.path("/Volumes/T7/Arbeit/FIA/processing_neg"),
                        linux = file.path("~/processing_neg")
                      ))
      
      updateTextInput(session, "raw_dir",
                      value = default_path(
                        win = file.path("D:", "Arbeit", "FIA", "processing_neg", "raw_mzML"),
                        mac = file.path("/Volumes/T7/Arbeit/FIA/processing_neg/raw_mzML"),
                        linux = file.path("~/processing_neg/raw_mzML")
                      ))
      
      updateTextInput(session, "data_dir_DDA",
                      value = default_path(
                        win = file.path("D:", "Arbeit", "FIA", "DDA", "processing_neg", "filtered_feature_df_annotated_neg.csv"),
                        mac = file.path("/Volumes/T7/Arbeit/FIA/DDA/processing_neg/filtered_feature_df_annotated_neg.csv"),
                        linux = file.path("~/filtered_feature_df_annotated_neg.csv")
                      ))
    }
  }, ignoreNULL = TRUE, ignoreInit = FALSE)
  
  
  # ---------- File explorer ----------
  refresh_file_list <- function() {
    dirpath <- if (nzchar(input$csv_explorer_dir)) input$csv_explorer_dir else input$data_dir
    if (is.null(dirpath) || !nzchar(dirpath)) {
      rv$file_list <- NULL; return()
    }
    if (dir.exists(dirpath)) {
      files <- list.files(dirpath, full.names = TRUE)
      if (length(files) == 0) {
        rv$file_list <- data.frame(name=character(0), path=character(0), type=character(0), size_kb=numeric(0))
      } else {
        rv$file_list <- data.frame(
          name = basename(files),
          path = files,
          type = tools::file_ext(files),
          size_kb = round(file.info(files)$size / 1024, 2),
          stringsAsFactors = FALSE
        )
      }
    } else {
      rv$file_list <- NULL
    }
  }
  
  observeEvent(input$refresh_files, { refresh_file_list(); append_log("File list refreshed") })
  observeEvent(input$csv_explorer_dir, { refresh_file_list() })
  output$file_browser <- renderDT({
    req(rv$file_list)
    datatable(rv$file_list, selection = "single", options = list(pageLength = 15))
  })
  
  observeEvent(input$file_browser_rows_selected, {
    row <- input$file_browser_rows_selected
    req(row)
    file <- rv$file_list$path[row]
    ext <- tools::file_ext(file)
    if (tolower(ext) == "csv") {
      tryCatch({
        rv$preview_data <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
      }, error = function(e) {
        rv$preview_data <- data.frame(error = paste("Failed to read CSV:", e$message))
      })
    } else rv$preview_data <- NULL
  })
  
  output$file_preview <- renderDT({
    req(rv$preview_data)
    datatable(rv$preview_data, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # ---------- Open output folder ----------
  observeEvent(input$open_output_folder, {
    outdir <- input$data_dir
    if (!nzchar(outdir) || !dir.exists(outdir)) {
      showNotification("Output directory does not exist", type = "error"); return()
    }
    if (OS_TYPE == "windows") shell.exec(normalizePath(outdir))
    else if (OS_TYPE == "macos") system2("open", args = c(shQuote(outdir)))
    else system2("xdg-open", args = c(shQuote(outdir)))
  })
  
  # ---------- Validate inputs before run ----------
  validate_inputs <- function() {
    errs <- character()
    if (!nzchar(input$script_path) || !file.exists(input$script_path)) errs <- c(errs, "Pipeline script not found.")
    if (!nzchar(input$data_dir) || !dir.exists(input$data_dir)) errs <- c(errs, "Processing data_dir not found.")
    if (!nzchar(input$raw_dir) || !dir.exists(input$raw_dir)) errs <- c(errs, "Raw files directory not found.")
    if (!nzchar(input$file_path) || !file.exists(input$file_path)) errs <- c(errs, "HMDB CSV not found.")
    mzmls <- list.files(input$raw_dir, pattern = "\\.(mzXML|mzML)$", ignore.case = TRUE)
    if (length(mzmls) == 0) errs <- c(errs, "No mzML files found in raw_dir.")
    if (length(errs)) return(paste(errs, collapse = "\n"))
    return(NULL)
  }
  
  # ---------- Progress UI ----------
  output$progress_ui <- renderUI({
    total <- rv$progress_total
    if (is.na(total) || is.null(total)) {
      tagList(
        tags$p("No progress information available."),
        tags$p(paste("Estimated progress:", round(rv$progress, 1), "%"))
      )
    } else {
      pct <- if (total > 0) round(100 * rv$progress / total) else 0
      tagList(
        tags$div(
          style = "width:100%; background:#eee; border-radius:4px; height:18px; overflow:hidden;",
          tags$div(style = sprintf("width:%s%%; height:100%%; background:steelblue;", pct))
        ),
        tags$p(paste("Processing:", rv$progress, "/", total, "files (", pct, "%)"))
      )
    }
  })
  
  # ---------- Run pipeline (with confirmation modal) ----------
  observeEvent(input$run_btn, {
    if (rv$running) { showNotification("Pipeline already running", type = "warning"); return() }
    # preflight
    err <- validate_inputs()
    if (!is.null(err)) {
      showModal(modalDialog(title = "Validation failed", paste(err, collapse = "\n"), easyClose = TRUE, footer = modalButton("Close")))
      append_log("Validation failed: \n" %+% err)
      return()
    }
    # estimate time
    mzmls <- list.files(input$raw_dir, pattern = "\\.(mzXML|mzML)$", ignore.case = TRUE, full.names = TRUE)
    est_sec <- length(mzmls) * rv$est_seconds_per_mzml
    est_text <- paste0("Estimated time: ~", round(est_sec/60,1), " minutes for ", length(mzmls), " files.")
    showModal(modalDialog(
      title = "Confirm Run",
      paste("You're about to start the pipeline.\n", est_text),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_run", "Start", class = "btn-primary")
      ), easyClose = TRUE
    ))
  })
  
  # confirm_run starts the background process
  observeEvent(input$confirm_run, {
    removeModal()
    # build args, write a small config file next to script
    cfg <- make_config_list()
    append_log("Starting pipeline with direct config injection")
    rv$progress <- 0
    rv$progress_total <- NA
    rv$running <- TRUE
    
    script <- input$script_path
    if (!file.exists(script)) {
      append_log("Script not found at run time. Aborting.")
      showNotification("Script not found", type = "error")
      rv$running <- FALSE
      return()
    }
    
    # spawn background R process that sources the script with config
    # We call Rscript and pass config path as an arg; the user's script should accept it, or we source it with local env below
    # Safer: run in background and capture stdout/stderr
    rv$r_process <- r_bg(
      function(script, config) {
        # create a clean environment for the script
        env <- new.env()
        # assign all config elements as variables
        list2env(config, envir = env)
        # optionally, print a log message
        cat("===== Pipeline started =====\n")
        # source the user script in that environment
        source(script, local = env, echo = TRUE)
      },
      args = list(script = script, config = cfg)
    )
    
    # monitor process in an observer
    observe({
      req(rv$running)  # only run while rv$running is TRUE
      invalidateLater(100, session)
      req(rv$r_process)
      
      if (rv$r_process$is_alive()) {
        out <- rv$r_process$read_output_lines()
        err <- rv$r_process$read_error_lines()
        if (length(out) && any(nzchar(out))) {
          for (ln in out[nzchar(out)]) {
            append_log(ln)
            # parse protocol lines: if user script emits "PROGRESS: i/N"
            if (grepl("^PROGRESS:", ln)) {
              nums <- as.numeric(unlist(regmatches(ln, gregexpr("[0-9]+", ln))))
              if (length(nums) >= 2) {
                rv$progress <- nums[1]
                rv$progress_total <- nums[2]
              }
            }
          }
        }
        if (length(err) && any(nzchar(err))) {
          for (ln in err[nzchar(err)]) append_log(paste0(">: ", ln))
        } else {
          # attempt to update estimated progress if no PROGRESS messages found
          # crude: use number of mzML files processed if the script prints processed filenames; detect filename lines
          if (is.na(rv$progress_total)) {
            mzmls <- list.files(input$raw_dir, pattern = "\\.(mzXML|mzML)$", ignore.case = TRUE, full.names = TRUE)
            if (length(mzmls) > 0) {
              rv$progress_total <- length(mzmls)
              # try to count how many mzML filenames appeared in log
              processed <- sum(sapply(mzmls, function(p) any(grepl(basename(p), tail(rv$log, 200), fixed = TRUE))))
              rv$progress <- processed
            }
          }
        }
      } else {
        # process finished - read final output
        out <- rv$r_process$read_output_lines()
        err <- rv$r_process$read_error_lines()
        if (length(out) && any(nzchar(out))) for (ln in out[nzchar(out)]) append_log(ln)
        if (length(err) && any(nzchar(err))) for (ln in err[nzchar(err)]) append_log(paste0("ERR: ", ln))
        
        # Reset state for next run
        rv$running <- FALSE
        rv$r_process <- NULL
        rv$progress <- 0
        rv$progress_total <- NA
        
        append_log("Pipeline finished. Ready for next run.")
        showNotification("Pipeline completed successfully!", type = "message", duration = 5)
        
        # Refresh file list to show new outputs
        refresh_file_list()
      }
    })
  }, ignoreInit = TRUE)
  
  # ---------- Stop pipeline ----------
  observeEvent(input$stop_btn, {
    req(rv$r_process)
    if (rv$r_process$is_alive()) {
      rv$r_process$kill()
      append_log("Pipeline stopped by user.")
      rv$running <- FALSE
      rv$r_process <- NULL
      rv$progress <- 0
      rv$progress_total <- NA
      showNotification("Pipeline stopped", type = "warning")
    } else {
      append_log("No running pipeline to stop.")
    }
  })
  
  # ---------- Download log ----------
  output$download_log <- downloadHandler(
    filename = function() paste0("fia_run_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"),
    content = function(file) writeLines(rv$log, con = file)
  )
}

shinyApp(ui, server)
