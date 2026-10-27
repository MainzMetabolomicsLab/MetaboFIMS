# MetaboFIMS SHINY version: 1.2.0
# MetaboFIMS SHINY application with first-start setup wizard.
# Start this app from the folder that should contain your data (the "working folder").
# On first start (or whenever paths.yml is missing / invalid / from another OS) a wizard asks to either
#   (a) set the four paths manually with the built-in file browser, or
#   (b) download the Demo Data from GitHub into the working folder and point the app to it.
# The result is stored in paths.yml next to the app, so the wizard only appears again when needed.

# ---- Repos & parallel install ----
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = max(1L, parallel::detectCores(logical = TRUE) - 1L))

packages <- list(
  bioc = c("MSnbase", "Rdisop", "BiocParallel", "Spectra", "S4Vectors", "mzR"),
  cran = c(
    # app
    "shiny", "shinyFiles", "DT", "callr", "yaml", "ps", "digest", "rlang", "jsonlite",
    # pipeline
    "data.table", "dplyr", "tidyverse", "ggplot2", "stringr", "signal",
    "R.utils", "webchem", "stringdist"
  )
)

is_installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# BiocManager first (also checks R/Bioc compatibility)
if (!is_installed("BiocManager")) install.packages("BiocManager")

missing_cran <- Filter(Negate(is_installed), packages$cran)
if (length(missing_cran)) {
  message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

missing_bioc <- Filter(Negate(is_installed), packages$bioc)
if (length(missing_bioc)) {
  message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# ---- Verify: stop early with a clear message ----
still_missing <- Filter(Negate(is_installed), unlist(packages))
if (length(still_missing)) {
  stop("Could not install: ", paste(still_missing, collapse = ", "),
       "\nCheck Rtools (Windows) / Xcode CLT (macOS) / system libraries (Linux).",
       call. = FALSE)
}

# ---- Attach only what the Shiny app itself needs ----
# (the pipeline runs in a separate r_bg process and loads its own packages)
library(shiny)
library(shinyFiles)
library(DT)
library(callr)
library(yaml)

# `%||%` is only in base R >= 4.4.0
`%||%` <- rlang::`%||%`

message("✅ All required packages are installed.")

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

get_os <- function() {
  sys <- Sys.info()[["sysname"]]
  if (.Platform$OS.type == "windows") return("windows")
  if (sys == "Darwin") return("macos")
  if (sys == "Linux") return("linux")
  return("unknown")
}
OS_TYPE <- get_os()

if (!dir.exists("logs")) dir.create("logs", recursive = TRUE)

# ======================================================================
# First-start setup: paths.yml, auto-detection, GitHub demo-data download
# ======================================================================

# Folder the app was started from. Demo data and hmdb.csv are placed here.
APP_DIR    <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
PATHS_FILE <- file.path(APP_DIR, "paths.yml")

# App settings that persist across starts (kept separate from paths.yml)
SETTINGS_FILE <- file.path(APP_DIR, "app_settings.yml")
read_app_settings <- function() {
  if (!file.exists(SETTINGS_FILE)) return(list())
  st <- tryCatch(yaml::read_yaml(SETTINGS_FILE), error = function(e) NULL)
  if (is.list(st)) st else list()
}
# Automatic update check at app start: on unless switched off (missing file = on)
read_auto_update <- function() !isFALSE(read_app_settings()$auto_update_check)
save_auto_update <- function(on) {
  st <- read_app_settings()
  st$auto_update_check <- isTRUE(on)
  yaml::write_yaml(st, SETTINGS_FILE)
}

GH_OWNER   <- "MainzMetabolomicsLab"
GH_REPO    <- "MetaboFIMS"
GH_BRANCH  <- "main"
GH_DEMO    <- "Demo Data"

PATH_KEYS   <- c("script_path", "data_dir", "data_dir_raw_files", "file_path")
PATH_LABELS <- c(script_path = "R script",
                 data_dir = "Processing data directory",
                 data_dir_raw_files = "Raw files directory",
                 file_path = "HMDB annotation CSV")

.s      <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) "" else as.character(x)
is_file <- function(x) { x <- .s(x); nzchar(x) && file.exists(x) && !dir.exists(x) }
is_dir  <- function(x) { x <- .s(x); nzchar(x) && dir.exists(x) }
norm_path <- function(x) sub("/+$", "", normalizePath(x, winslash = "/", mustWork = FALSE))

# TRUE/FALSE per required path
paths_status <- function(p) c(
  script_path        = is_file(p$script_path),
  data_dir           = is_dir(p$data_dir),
  data_dir_raw_files = is_dir(p$data_dir_raw_files),
  file_path          = is_file(p$file_path)
)

# Reads paths.yml. Returns an empty list if missing, unreadable, or written on a different OS.
read_saved_paths <- function() {
  if (!file.exists(PATHS_FILE)) return(list())
  p <- tryCatch(yaml::read_yaml(PATHS_FILE), error = function(e) NULL)
  if (is.null(p) || !identical(p$os, OS_TYPE)) return(list())
  p[intersect(PATH_KEYS, names(p))]
}

save_saved_paths <- function(p) {
  out <- c(lapply(p[PATH_KEYS], .s), list(os = OS_TYPE, saved = format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  yaml::write_yaml(out, PATHS_FILE)
}

# Pipeline script sitting next to the app: "MetaboFIMS.R" first; otherwise (older installs, e.g.
# "MetaboFIMS V1.1.R") the highest MetaboFIMS*.R name wins. App files themselves are ignored.
PIPELINE_SCRIPT_NAME <- "MetaboFIMS.R"
find_pipeline_script <- function() {
  exact <- file.path(APP_DIR, PIPELINE_SCRIPT_NAME)
  if (file.exists(exact)) return(exact)
  f <- list.files(APP_DIR, pattern = "^MetaboFIMS.*\\.[Rr]$", full.names = TRUE)
  f <- f[!grepl("shiny|app", basename(f), ignore.case = TRUE)]
  if (length(f)) sort(f, decreasing = TRUE)[1] else ""
}

# The Shiny app itself ("MetaboFIMS SHINY.R"); any MetaboFIMS*.R with "shiny"/"app" in its name counts as app file.
APP_SCRIPT_NAME <- "MetaboFIMS SHINY.R"
is_app_file <- function(name) grepl("^MetaboFIMS.*\\.[Rr]$", name) & grepl("shiny|app", name, ignore.case = TRUE)
find_local_app <- function() {
  exact <- file.path(APP_DIR, APP_SCRIPT_NAME)
  if (file.exists(exact)) return(exact)
  f <- list.files(APP_DIR, full.names = TRUE)
  f <- f[is_app_file(basename(f))]
  if (length(f)) sort(f, decreasing = TRUE)[1] else ""
}

# HMDB csv sitting directly in the working folder
find_hmdb_csv <- function() {
  f <- list.files(APP_DIR, pattern = "^hmdb.*\\.csv$", full.names = TRUE, ignore.case = TRUE)
  if (length(f)) f[1] else ""
}

# Downloads everything below "Demo Data/" of the GitHub repo into `dest` (folder structure kept,
# the "Demo Data" prefix is dropped). Uses ONE API call for the file list, then raw downloads.
# Already complete files are skipped, so an interrupted download can simply be restarted.
download_demo_data <- function(dest, progress = function(i, n, f) NULL) {
  old <- options(timeout = max(3600, getOption("timeout"))); on.exit(options(old), add = TRUE)
  api <- sprintf("https://api.github.com/repos/%s/%s/git/trees/%s?recursive=1", GH_OWNER, GH_REPO, GH_BRANCH)
  res <- tryCatch(jsonlite::fromJSON(api),
                  error = function(e) stop("Could not read the GitHub file list (offline, or API rate limit reached?): ",
                                           conditionMessage(e), call. = FALSE))
  if (isTRUE(res$truncated)) warning("GitHub returned a truncated file list - some demo files may be missing.", call. = FALSE)
  tree <- res$tree
  prefix <- paste0(GH_DEMO, "/")
  blobs  <- tree[tree$type == "blob" & startsWith(tree$path, prefix), , drop = FALSE]
  # The Shiny app is not part of the demo download (it is updated via "Check for updates" instead)
  blobs  <- blobs[!is_app_file(basename(blobs$path)), , drop = FALSE]
  if (!nrow(blobs)) stop("No files found below '", GH_DEMO, "' in the GitHub repository.", call. = FALSE)
  
  enc <- function(path) paste(vapply(strsplit(path, "/", fixed = TRUE)[[1]],
                                     utils::URLencode, "", reserved = TRUE), collapse = "/")
  raw_url   <- function(path) sprintf("https://raw.githubusercontent.com/%s/%s/%s/%s", GH_OWNER, GH_REPO, GH_BRANCH, enc(path))
  media_url <- function(path) sprintf("https://media.githubusercontent.com/media/%s/%s/%s/%s", GH_OWNER, GH_REPO, GH_BRANCH, enc(path))
  is_lfs_pointer <- function(f) {
    file.size(f) < 1024 && isTRUE(startsWith(tryCatch(readLines(f, n = 1, warn = FALSE), error = function(e) ""),
                                             "version https://git-lfs"))
  }
  
  rel     <- substring(blobs$path, nchar(prefix) + 1)
  targets <- file.path(dest, rel)
  n <- length(rel)
  for (i in seq_len(n)) {
    progress(i, n, rel[i])
    dir.create(dirname(targets[i]), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(targets[i]) && isTRUE(file.size(targets[i]) == blobs$size[i])) next
    utils::download.file(raw_url(blobs$path[i]), targets[i], mode = "wb", quiet = TRUE)
    if (is_lfs_pointer(targets[i])) {   # repo stores this file with Git LFS
      utils::download.file(media_url(blobs$path[i]), targets[i], mode = "wb", quiet = TRUE)
    }
  }
  rel
}

# Works out the four app paths from the files that were just downloaded.

# Works out the four app paths from the files that were just downloaded.
resolve_demo_paths <- function(dest, rel_files) {
  # HMDB CSV: must end up directly in the working folder
  hm <- rel_files[grepl("^hmdb.*\\.csv$", basename(rel_files), ignore.case = TRUE)]
  file_path <- ""
  if (length(hm)) {
    top <- file.path(dest, basename(hm[1]))
    if (!file.exists(top)) file.copy(file.path(dest, hm[1]), top)
    file_path <- top
  }
  
  # Raw files folder: folder holding the most mzML/mzXML files.
  # Prefer non-negative-mode folders where possible.
  ms <- rel_files[grepl("\\.(mzml|mzxml)$", rel_files, ignore.case = TRUE)]
  raw_dir <- ""
  data_dir <- ""
  
  if (length(ms)) {
    cand <- names(sort(table(dirname(ms)), decreasing = TRUE))
    nonneg <- cand[!grepl("neg", cand, ignore.case = TRUE)]
    rel_dir <- if (length(nonneg)) nonneg[1] else cand[1]
    
    raw_dir <- if (rel_dir == ".") dest else file.path(dest, rel_dir)
    
    # If raw files are in a dedicated raw_mzML folder,
    # use its parent as the processing directory.
    if (tolower(basename(raw_dir)) == "raw_mzml") {
      data_dir <- dirname(raw_dir)
    } else {
      # Backward compatibility: mzML files directly in processing
      data_dir <- raw_dir
    }
  }
  
  list(
    data_dir = if (nzchar(data_dir)) norm_path(data_dir) else "",
    data_dir_raw_files = if (nzchar(raw_dir)) norm_path(raw_dir) else "",
    file_path = if (nzchar(file_path)) norm_path(file_path) else ""
  )
}


# ---- Update check: compare the local pipeline script with the one on GitHub (git blob SHA-1) ----
# A different hash is taken as "a newer version is available".
git_blob_sha <- function(bytes) {
  digest::digest(c(charToRaw(sprintf("blob %d", length(bytes))), as.raw(0L), bytes),
                 algo = "sha1", serialize = FALSE)
}
# Hashes of a local file: as is, and with CRLF -> LF (Windows checkouts / editors may convert line endings)
local_script_shas <- function(path) {
  b  <- readBin(path, "raw", file.size(path))
  lf <- b[!(b == as.raw(13L) & c(b[-1L], as.raw(0L)) == as.raw(10L))]
  unique(c(git_blob_sha(b), git_blob_sha(lf)))
}
# Version lines in the file headers:
#   pipeline: "# MetaboFIMS version: 1.2.0"        (old style "#V1.1" is still recognised)
#   app:      "## MetaboFIMS SHINY version: 1.2.0"
VERSION_RX <- list(
  pipeline = c("^#+\\s*MetaboFIMS version:\\s*([0-9]+(\\.[0-9]+)*)", "^#\\s*V([0-9]+(\\.[0-9]+)*)"),
  app      = c("^#+\\s*MetaboFIMS SHINY version:\\s*([0-9]+(\\.[0-9]+)*)")
)
COMPONENT_LABELS <- c(pipeline = "Pipeline script", app = "Shiny app")

# Version from the file header. Returns NA if no version line is found in the first 30 lines.
read_script_version <- function(path, rx = VERSION_RX$pipeline) {
  if (!is_file(path)) return(NA_character_)
  l <- tryCatch(readLines(path, n = 30, warn = FALSE, encoding = "UTF-8"), error = function(e) character(0))
  l <- sub("^\ufeff", "", l)   # byte order mark written by some Windows editors
  for (r in rx) {
    m <- regmatches(l, regexec(r, l, ignore.case = TRUE))
    hit <- vapply(m, function(x) if (length(x) >= 2) x[2] else NA_character_, "")
    hit <- hit[!is.na(hit)]
    if (length(hit)) return(hit[1])
  }
  NA_character_
}
fmt_version <- function(v) if (is.null(v) || length(v) != 1L || is.na(v)) "unknown" else v

# File list of the GitHub repo (one API call, shared by all components of one check)
gh_tree <- function() {
  api <- sprintf("https://api.github.com/repos/%s/%s/git/trees/%s?recursive=1", GH_OWNER, GH_REPO, GH_BRANCH)
  jsonlite::fromJSON(api)$tree
}

# Exact content of a file version on GitHub, fetched by its git blob SHA through the API.
# Unlike raw.githubusercontent.com this is never served from a stale cache right after an upload.
gh_blob_bytes <- function(sha) {
  api <- sprintf("https://api.github.com/repos/%s/%s/git/blobs/%s", GH_OWNER, GH_REPO, sha)
  b <- jsonlite::fromJSON(api)
  if (!identical(b$encoding, "base64")) stop("unexpected blob encoding from GitHub", call. = FALSE)
  jsonlite::base64_dec(gsub("[\r\n]", "", b$content))
}
# Writes the GitHub file with this SHA to `dest`: via the blob API, falling back to the raw URL.
# Returns TRUE only if the written file really has that SHA.
fetch_online_file <- function(sha, url, dest) {
  ok <- tryCatch({ writeBin(gh_blob_bytes(sha), dest); TRUE }, error = function(e) FALSE)
  if (!ok || !(sha %in% local_script_shas(dest)))
    ok <- tryCatch({ utils::download.file(url, dest, mode = "wb", quiet = TRUE); TRUE }, error = function(e) FALSE)
  ok && file.exists(dest) && sha %in% local_script_shas(dest)
}

# Pipeline script ("Demo Data/MetaboFIMS.R") or Shiny app ("Demo Data/MetaboFIMS SHINY.R") on GitHub;
# without the exact name the highest matching name wins. Also reads the version from its header.
find_online_script <- function(kind = "pipeline", tree = gh_tree()) {
  nm   <- basename(tree$path)
  is_kind <- if (kind == "app") is_app_file(nm) else grepl("^MetaboFIMS.*\\.[Rr]$", nm) & !is_app_file(nm)
  hit  <- tree$type == "blob" & startsWith(tree$path, paste0(GH_DEMO, "/")) & is_kind
  b    <- tree[hit, , drop = FALSE]
  if (!nrow(b)) stop("No MetaboFIMS ", if (kind == "app") "Shiny app" else "pipeline script",
                     " found in '", GH_DEMO, "' of the GitHub repository.", call. = FALSE)
  exact <- basename(b$path) == (if (kind == "app") APP_SCRIPT_NAME else PIPELINE_SCRIPT_NAME)
  b    <- if (any(exact)) b[which(exact)[1], , drop = FALSE] else b[order(basename(b$path), decreasing = TRUE)[1], , drop = FALSE]
  enc  <- paste(vapply(strsplit(b$path, "/", fixed = TRUE)[[1]], utils::URLencode, "", reserved = TRUE), collapse = "/")
  url  <- sprintf("https://raw.githubusercontent.com/%s/%s/%s/%s", GH_OWNER, GH_REPO, GH_BRANCH, enc)
  # Online version: read from exactly the file version with this SHA
  version <- NA_character_
  tmp <- tempfile(fileext = ".R"); on.exit(unlink(tmp), add = TRUE)
  if (fetch_online_file(b$sha, url, tmp)) version <- read_script_version(tmp, VERSION_RX[[kind]])
  list(path = b$path, name = basename(b$path), sha = b$sha, url = url, version = version)
}

# Compares one local file with its GitHub counterpart. state: "up to date", "update available",
# "local file modified" (same version, other content), "local version is newer", "no local file", "check failed"
assess_update <- function(kind, local_path, tree) {
  res <- list(kind = kind, local = local_path,
              local_version = read_script_version(local_path, VERSION_RX[[kind]]), version = NA_character_)
  if (!is_file(local_path)) return(c(res, list(state = "no local file")))
  on <- tryCatch(find_online_script(kind, tree), error = function(e) e)
  if (inherits(on, "error")) return(c(res, list(state = "check failed", error = conditionMessage(on))))
  res[names(on)] <- on
  if (on$sha %in% local_script_shas(local_path)) return(c(res, list(state = "up to date")))
  cmp <- if (!is.na(res$local_version) && !is.na(on$version)) utils::compareVersion(on$version, res$local_version) else NA
  state <- if (isTRUE(cmp == -1)) "local version is newer" else if (isTRUE(cmp == 0)) "local file modified" else "update available"
  c(res, list(state = state))
}

# Switching Positive <-> Negative: use the sibling "<data_dir>_neg" folder if it exists
swap_polarity_dirs <- function(dd, rd, to) {
  dd <- sub("/+$", "", dd); rd <- sub("/+$", "", rd)
  base   <- sub("_neg$", "", dd)
  new_dd <- if (to == "Negative") paste0(base, "_neg") else base
  if (!dir.exists(new_dd) || new_dd == dd) return(NULL)
  new_rd <- if (startsWith(rd, dd)) paste0(new_dd, substring(rd, nchar(dd) + 1)) else rd
  if (!dir.exists(new_rd)) new_rd <- rd
  list(data_dir = new_dd, raw_dir = new_rd)
}

# ---- Output folder: <data_dir>/<custom name or YYYY-MM-DD_HH-MM-SS> ----
clean_folder_name <- function(x) {
  x <- trimws(.s(x))
  x <- gsub('[\\\\/:*?"<>|]', "_", x)   # characters not allowed in folder names
  x <- gsub("^\\.+|\\.+$", "", x)       # no leading / trailing dots
  trimws(x)
}
make_output_dir <- function(data_dir, name = "") {
  nm <- clean_folder_name(name)
  if (!nzchar(nm)) nm <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  file.path(sub("/+$", "", data_dir), nm)
}
# Subfolders created inside every output folder, e.g. c("annotation", "peaklists").
# They are passed to the pipeline as the named list `output_subdirs`.
OUTPUT_SUBFOLDERS <- character(0)

# ---- What the UI starts with: saved paths if still valid, otherwise whatever can be auto-detected ----
saved_paths <- read_saved_paths()
initial_paths <- list(
  script_path        = if (is_file(saved_paths$script_path)) .s(saved_paths$script_path) else find_pipeline_script(),
  data_dir           = if (is_dir(saved_paths$data_dir)) .s(saved_paths$data_dir) else "",
  data_dir_raw_files = if (is_dir(saved_paths$data_dir_raw_files)) .s(saved_paths$data_dir_raw_files) else "",
  file_path          = if (is_file(saved_paths$file_path)) .s(saved_paths$file_path) else find_hmdb_csv()
)
NEEDS_SETUP <- !all(paths_status(initial_paths))

setup_modal <- function(msg = NULL) {
  modalDialog(
    title = "First-start setup",
    p("Detected operating system: ", strong(switch(OS_TYPE, windows = "Windows", macos = "macOS",
                                                   linux = "Linux", "unknown"))),
    p("Working folder (the app was started from here): ", tags$code(APP_DIR)),
    if (!is.null(msg)) tags$p(style = "color:#b71c1c;", msg),
    p("The R script, the processing folder, the raw files folder and the HMDB CSV are not set yet."),
    radioButtons(
      "setup_mode", "How do you want to set them up?",
      choices = c("Set the paths manually (built-in file browser)" = "manual",
                  "Download the Demo Data from GitHub into the working folder" = "download"),
      selected = "manual"
    ),
    helpText("Download: all files of 'Demo Data' are copied into the working folder, hmdb.csv is placed directly in it, ",
             "and the folders are filled in automatically. The pipeline R script is searched next to the app; ",
             "if it is not found you will be asked to select it."),
    footer = tagList(
      actionButton("setup_skip", "Later"),
      actionButton("setup_go", "Continue", class = "btn-primary")
    ),
    easyClose = FALSE
  )
}

# Reads the distinct Biospecimen / Origin values from the HMDB CSV (only these two columns are loaded).
# Cells with several values (separated by ; , | /) are split. Returns NULL if the CSV or columns are missing.
read_hmdb_levels <- function(csv_path) {
  if (is.null(csv_path) || !nzchar(csv_path) || !file.exists(csv_path)) return(NULL)
  hdr <- tryCatch(names(data.table::fread(csv_path, nrows = 0)), error = function(e) NULL)
  if (is.null(hdr)) return(NULL)
  cols <- intersect(c("Biospecimen", "Origin"), hdr)
  if (length(cols) < 2) return(NULL)
  dt <- tryCatch(
    data.table::fread(csv_path, select = cols, colClasses = "character", data.table = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt)) return(NULL)
  
  distinct_values <- function(x) {
    v <- trimws(unlist(strsplit(as.character(x), "[;,|/]")))
    v <- v[!is.na(v) & nzchar(v)]
    v <- v[!duplicated(tolower(v))]   # "blood" and "Blood" count once
    sort(v)
  }
  list(biospecimen = distinct_values(dt$Biospecimen),
       origin      = distinct_values(dt$Origin))
}

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
      uiOutput("setup_status"),
      fluidRow(
        column(6, actionButton("open_setup", "Run setup wizard")),
        column(6, actionButton("save_paths", "Save paths as default"))
      ),
      br(),
      actionButton("check_updates", "Check for updates"),
      checkboxInput("auto_update", "Check for updates automatically at app start", value = read_auto_update()),
      uiOutput("version_status"),
      br(),
      fluidRow(
        column(9, textInput("script_path", "R Script", value = initial_paths$script_path)),
        column(3, shinyFilesButton("browse_script", "Browse", "Choose R script", multiple = FALSE))
      ),
      fluidRow(
        column(9, textInput("data_dir", "Processing data directory (data_dir)", value = initial_paths$data_dir)),
        column(3, shinyDirButton("browse_data_dir", "Browse", "Choose folder"))
      ),
      fluidRow(
        column(9, textInput("raw_dir", "Raw files directory (data_dir_raw_files)", value = initial_paths$data_dir_raw_files)),
        column(3, shinyDirButton("browse_raw_dir", "Browse", "Choose folder"))
      ),
      fluidRow(
        column(9, textInput("file_path", "HMDB annotation CSV (file_path)", value = initial_paths$file_path)),
        column(3, shinyFilesButton("browse_file_path", "Browse", "Choose CSV", multiple = FALSE))
      ),
      textInput("output_name", "Output folder name (optional)", value = "",
                placeholder = "leave empty = date_time folder"),
      uiOutput("output_dir_preview"),
      helpText("Biospecimen and origin options are read from the selected HMDB CSV."),
      checkboxGroupInput(
        "hmdb_biospecimen", "HMDB biospecimen",
        choices  = character(0),
        selected = NULL,
        inline   = TRUE
      ),
      checkboxGroupInput(
        "hmdb_origin", "HMDB origin",
        choices  = character(0),
        selected = NULL,
        inline   = TRUE
      ),
      hr(),
      h4("RT integration window (seconds) for 1D spectra accumulation"),
      numericInput("rt_min", "rt_min (s)", value = 5, min = 0, step = 1),
      numericInput("rt_max", "rt_max (s)", value = 25, min = 0, step = 1),
      hr(),
      h4("Mass range (m/z) limits"),
      numericInput("mz_min_limit", "mz_min_limit", value = 70, min = 0, step = 1),
      numericInput("mz_max_limit", "mz_max_limit", value = 1000, min = 0, step = 1),
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
      numericInput("ppm_tolerance_annotation", "PPM tolerance (annotation)", value = 10, min = 0, step = 1),
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

server <- function(input, output, session) {
  # ---------- roots for shinyFiles ----------
  # "App folder" comes first so the working folder is one click away in every file browser
  if (OS_TYPE == "windows") {
    roots <- c("App folder" = APP_DIR, home = normalizePath("~"), "C:" = "C:/", "D:" = "D:/", "E:" = "E:/", "F:" = "F:/")
  } else if (OS_TYPE == "macos") {
    roots <- c("App folder" = APP_DIR, home = normalizePath("~"), root = "/", Volumes = "/Volumes")
  } else {
    roots <- c("App folder" = APP_DIR, home = normalizePath("~"), root = "/")
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
  
  rv <- reactiveValues(
    log = character(),
    running = FALSE,
    r_process = NULL,
    file_list = NULL,
    preview_data = NULL,
    progress = 0,
    progress_total = NA,
    est_seconds_per_mzml = 15,
    last_output_dir = NULL
  )
  
  append_log <- function(msg) {
    line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", msg)
    # isolate(): append_log() may also be called outside a reactive consumer (e.g. at session start)
    cur <- c(isolate(rv$log), line)
    if (length(cur) > 10000) cur <- tail(cur, 10000)
    rv$log <- cur
    logfn <- file.path("logs", paste0(Sys.Date(), "_fia_run_log.txt"))
    cat(line, file = logfn, sep = "\n", append = TRUE)
  }
  
  output$console_out <- renderText({ paste(tail(rv$log, 200), collapse = "\n") })
  
  # ======================= First-start setup =======================
  setup_pending <- reactiveVal(NEEDS_SETUP)
  
  current_paths <- reactive(list(
    script_path        = input$script_path,
    data_dir           = input$data_dir,
    data_dir_raw_files = input$raw_dir,
    file_path          = input$file_path
  ))
  current_paths_d <- debounce(current_paths, 800)
  
  output$setup_status <- renderUI({
    st <- paths_status(current_paths())
    if (all(st)) {
      tags$p(style = "color:#2e7d32;", "\u2714 All paths are set.")
    } else {
      tags$p(style = "color:#b71c1c;", "\u2716 Missing or invalid: ", paste(PATH_LABELS[!st], collapse = ", "))
    }
  })
  
  output$output_dir_preview <- renderUI({
    dd <- .s(input$data_dir)
    if (!nzchar(dd)) return(NULL)
    nm <- clean_folder_name(input$output_name)
    shown <- file.path(sub("/+$", "", dd), if (nzchar(nm)) nm else "<date_time of run start>")
    helpText("Output folder of the next run: ", tags$code(shown))
  })
  
  if (NEEDS_SETUP) {
    append_log(paste0("First start / incomplete paths (OS: ", OS_TYPE, ") - opening setup wizard."))
    session$onFlushed(function() showModal(setup_modal(), session = session), once = TRUE)
  } else {
    append_log(paste0("Paths loaded from ", PATHS_FILE, " (OS: ", OS_TYPE, ")."))
  }
  
  observeEvent(input$open_setup, { showModal(setup_modal()) })
  observeEvent(input$setup_skip, { removeModal() })
  
  observeEvent(input$setup_go, {
    removeModal()
    if (identical(input$setup_mode, "manual")) {
      setup_pending(TRUE)
      showNotification(
        "Use the Browse buttons to set the R script, processing folder, raw files folder and HMDB CSV. They are saved automatically as soon as all four are valid.",
        type = "message", duration = 20)
      append_log("Setup: manual mode selected.")
      return()
    }
    
    append_log(paste0("Setup: downloading Demo Data from GitHub into ", APP_DIR))
    res <- tryCatch({
      rel <- withProgress(message = "Downloading Demo Data from GitHub", value = 0, {
        download_demo_data(APP_DIR, progress = function(i, n, f) {
          setProgress(value = (i - 1) / n, detail = sprintf("%d/%d  %s", i, n, f))
        })
      })
      resolve_demo_paths(APP_DIR, rel)
    }, error = function(e) e)
    
    if (inherits(res, "error")) {
      append_log(paste("Setup download failed:", conditionMessage(res)))
      showModal(setup_modal(msg = paste("Download failed:", conditionMessage(res))))
      return()
    }
    
    script <- if (is_file(input$script_path)) input$script_path else find_pipeline_script()
    updateTextInput(session, "script_path", value = script)
    if (nzchar(res$data_dir))           updateTextInput(session, "data_dir", value = res$data_dir)
    if (nzchar(res$data_dir_raw_files)) updateTextInput(session, "raw_dir", value = res$data_dir_raw_files)
    if (nzchar(res$file_path))          updateTextInput(session, "file_path", value = res$file_path)
    append_log(paste0("Demo data ready. data_dir=", res$data_dir, " | raw=", res$data_dir_raw_files,
                      " | HMDB=", res$file_path))
    
    p <- list(script_path = script, data_dir = res$data_dir,
              data_dir_raw_files = res$data_dir_raw_files, file_path = res$file_path)
    if (all(paths_status(p))) {
      save_saved_paths(p)
      setup_pending(FALSE)
      showNotification("Demo data downloaded and paths saved.", type = "message")
    } else {
      setup_pending(TRUE)   # the auto-save observer below finishes the job once the rest is set
      showNotification(
        paste0("Demo data downloaded. Still to set: ", paste(PATH_LABELS[!paths_status(p)], collapse = ", "),
               ". Paths are saved automatically once all four are valid."),
        type = "warning", duration = 20)
    }
  })
  
  # While setup is open: save paths.yml as soon as all four are valid
  observe({
    p <- current_paths_d()
    req(setup_pending())
    if (all(paths_status(p))) {
      save_saved_paths(p)
      setup_pending(FALSE)
      showNotification("All paths valid - saved to paths.yml.", type = "message")
      append_log(paste("Setup finished, paths saved to", PATHS_FILE))
    }
  })
  
  # Manual save of the current paths as new default
  observeEvent(input$save_paths, {
    p <- current_paths()
    st <- paths_status(p)
    if (!all(st)) {
      showNotification(paste("Cannot save, missing or invalid:", paste(PATH_LABELS[!st], collapse = ", ")), type = "error")
      return()
    }
    if (identical(input$Ionization_method, "Negative")) {   # always store the positive-mode folders
      p$data_dir <- sub("_neg(?=/|$)", "", p$data_dir, perl = TRUE)
      p$data_dir_raw_files <- sub("_neg(?=/|$)", "", p$data_dir_raw_files, perl = TRUE)
    }
    save_saved_paths(p)
    showNotification("Paths saved as default (paths.yml).", type = "message")
    append_log(paste("Paths saved to", PATHS_FILE))
  })
  
  # ======================= Update check (pipeline script vs. GitHub) =======================
  update_info  <- reactiveVal(NULL)
  online_state <- reactiveVal(list())   # per component: list(version, status) of the last check
  dialog_auto  <- reactiveVal(FALSE)    # TRUE while the dialog of the automatic start check is open
  
  output$version_status <- renderUI({
    st <- online_state()
    line <- function(kind, local_path) {
      s <- st[[kind]]
      tags$div(COMPONENT_LABELS[[kind]], " - local: ", tags$code(fmt_version(read_script_version(local_path, VERSION_RX[[kind]]))),
               " | GitHub: ", tags$code(if (is.null(s)) "not checked" else fmt_version(s$version)),
               if (!is.null(s)) paste0(" (", s$status, ")"))
    }
    helpText(line("pipeline", input$script_path), line("app", find_local_app()))
  })
  
  check_for_updates <- function(manual = FALSE) {
    tree <- tryCatch(gh_tree(), error = function(e) e)
    if (inherits(tree, "error")) {
      online_state(list(pipeline = list(version = NA_character_, status = "check failed"),
                        app      = list(version = NA_character_, status = "check failed")))
      append_log(paste("Update check skipped (offline or GitHub not reachable):", conditionMessage(tree)))
      if (manual) showNotification(paste("Update check failed:", conditionMessage(tree)), type = "error", session = session)
      return(invisible(NULL))
    }
    res <- list(pipeline = assess_update("pipeline", isolate(.s(input$script_path)), tree),
                app      = assess_update("app", find_local_app(), tree))
    online_state(lapply(res, function(r) list(version = r$version, status = r$state)))
    for (r in res) append_log(paste0("Update check - ", COMPONENT_LABELS[[r$kind]], ": ", r$state,
                                     " (local ", fmt_version(r$local_version), ", GitHub ", fmt_version(r$version), ")",
                                     if (!is.null(r$error)) paste0(" - ", r$error) else ""))
    
    todo <- Filter(function(r) r$state %in% c("update available", "local file modified"), res)
    if (!length(todo)) {
      if (manual) showNotification(paste(vapply(res, function(r) paste0(COMPONENT_LABELS[[r$kind]], ": ", r$state), ""),
                                         collapse = " | "), type = "message", session = session)
      return(invisible(NULL))
    }
    update_info(todo)
    dialog_auto(!manual)
    
    choice_names <- lapply(todo, function(r) {
      tagList(strong(COMPONENT_LABELS[[r$kind]]), ": ", tags$code(fmt_version(r$local_version)), " \u2192 ",
              tags$code(fmt_version(r$version)), " (", r$name, ")",
              if (r$state == "local file modified") tags$span(style = "color:#b71c1c;",
                                                              " - same version number, but your local file was modified; updating overwrites your changes"))
    })
    showModal(modalDialog(
      title = "Updates available on GitHub",
      p("Select what should be updated. The selected files are downloaded and replace your local copies."),
      checkboxGroupInput("update_items", NULL,
                         choiceNames  = unname(choice_names),
                         choiceValues = names(todo),
                         selected     = names(todo)[vapply(todo, function(r) r$state == "update available", TRUE)]),
      if ("app" %in% names(todo)) helpText("A Shiny app update takes effect after restarting the app."),
      if (!manual) tagList(
        hr(),
        checkboxInput("disable_auto", "Don't check for updates automatically at start anymore", value = FALSE),
        helpText("You can still check with 'Check for updates' and switch the automatic check on again there.")
      ),
      footer = tagList(
        actionButton("update_keep", "Keep my versions"),
        actionButton("update_go", "Update selected", class = "btn-primary")
      ),
      easyClose = FALSE
    ), session = session)
  }
  
  # Automatic check at start (not while the setup wizard is needed, and only if not switched off)
  if (!NEEDS_SETUP && read_auto_update()) {
    session$onFlushed(function() check_for_updates(manual = FALSE), once = TRUE)
  } else if (!NEEDS_SETUP) {
    append_log("Automatic update check at start is switched off - use 'Check for updates'.")
  }
  
  observeEvent(input$auto_update, {
    save_auto_update(input$auto_update)
    append_log(paste0("Automatic update check at start ", if (isTRUE(input$auto_update)) "switched on." else "switched off."))
    showNotification(if (isTRUE(input$auto_update)) "Automatic update check at start: on" else
      "Automatic update check at start: off - use 'Check for updates' to check manually.", type = "message")
  }, ignoreInit = TRUE)
  observeEvent(input$check_updates, { check_for_updates(manual = TRUE) })
  
  # Opt-out ticked in the automatic dialog -> switch the setting off (the checkbox observer saves it)
  apply_disable_auto <- function() {
    if (isTRUE(dialog_auto()) && isTRUE(input$disable_auto))
      updateCheckboxInput(session, "auto_update", value = FALSE)
    dialog_auto(FALSE)
  }
  
  observeEvent(input$update_keep, {
    removeModal()
    apply_disable_auto()
    append_log("Update declined - continuing with the local versions.")
  })
  
  # Downloads one component, checks its hash and replaces the local file. Returns the new path.
  replace_with_online <- function(r) {
    tmp <- tempfile(fileext = ".R"); on.exit(unlink(tmp), add = TRUE)
    if (!fetch_online_file(r$sha, r$url, tmp)) stop("could not download the GitHub version (or its hash did not match)", call. = FALSE)
    dest <- file.path(dirname(r$local), r$name)
    if (!file.copy(tmp, dest, overwrite = TRUE)) stop("could not write ", dest, call. = FALSE)
    if (norm_path(dest) != norm_path(r$local)) unlink(r$local)   # old file replaced by the new one
    norm_path(dest)
  }
  
  observeEvent(input$update_go, {
    sel  <- input$update_items
    info <- update_info()
    removeModal()
    apply_disable_auto()
    req(info)
    if (!length(sel)) { append_log("Update: nothing selected."); return() }
    if (rv$running) {
      showNotification("A pipeline run is in progress - update after it has finished.", type = "warning")
      return()
    }
    st <- online_state()
    new_app <- NULL
    for (k in sel) {
      r <- info[[k]]
      dest <- tryCatch(replace_with_online(r), error = function(e) e)
      if (inherits(dest, "error")) {
        append_log(paste0("Update failed (", COMPONENT_LABELS[[k]], "): ", conditionMessage(dest)))
        showNotification(paste0("Update failed (", COMPONENT_LABELS[[k]], "): ", conditionMessage(dest)), type = "error")
        next
      }
      if (k == "pipeline") {
        updateTextInput(session, "script_path", value = dest)
        saved <- read_saved_paths()
        if (length(saved)) { saved$script_path <- dest; save_saved_paths(saved) }
      }
      if (k == "app") new_app <- dest
      st[[k]] <- list(version = r$version, status = "up to date")
      append_log(paste0(COMPONENT_LABELS[[k]], " updated from GitHub: ", basename(r$local), " (", fmt_version(r$local_version),
                        ") -> ", r$name, " (", fmt_version(r$version), ")"))
      if (k != "app") showNotification(paste0(COMPONENT_LABELS[[k]], " updated to ", r$name, " (", fmt_version(r$version), ")"),
                                       type = "message")
    }
    online_state(st)
    update_info(NULL)
    if (!is.null(new_app)) {
      new_app_path(new_app)
      showModal(modalDialog(
        title = "Shiny app updated",
        p("The new app version (", tags$code(fmt_version(st$app$version)), ") has been saved. ",
          "It takes effect after a restart."),
        p("'Restart now' restarts the R session in RStudio and starts the new app automatically. ",
          "This window is then no longer connected - the new app opens like with 'Run App' in RStudio."),
        footer = tagList(
          actionButton("restart_later", "Later"),
          actionButton("restart_app", "Restart now", class = "btn-primary")
        ),
        easyClose = FALSE
      ))
    }
  })
  
  # ---- Restart after an app update (RStudio: restart the R session, then run the new app) ----
  new_app_path <- reactiveVal(NULL)
  
  observeEvent(input$restart_later, {
    removeModal()
    append_log("App restart postponed - the new app version is used after the next start.")
    showNotification("The new app version is used after the next start.", type = "message")
  })
  
  observeEvent(input$restart_app, {
    removeModal()
    app_file <- new_app_path()
    req(app_file)
    if (rv$running) {
      showNotification("A pipeline run is in progress - restart after it has finished.", type = "warning")
      return()
    }
    can_restart <- requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
      "command" %in% names(formals(rstudioapi::restartSession))
    if (!can_restart) {
      append_log("Automatic restart not available (needs RStudio and a recent 'rstudioapi') - please restart manually.")
      showNotification("Automatic restart needs RStudio - please stop the app and start it again.",
                       type = "warning", duration = NULL)
      return()
    }
    # After the restart only setwd() runs as restart command, so RStudio can finish "Restoring R session".
    # runApp() is then sent to the console like typed code and opens where RStudio's "Run App" setting says
    # (no launch.browser here, so RStudio's own viewer setting - e.g. "Run in Window" - is used).
    run_cmd <- sprintf("shiny::runApp(%s)", deparse(norm_path(app_file)))
    cmd <- sprintf("setwd(%s); rstudioapi::sendToConsole(%s, execute = TRUE, focus = FALSE)",
                   deparse(APP_DIR), deparse(run_cmd))
    append_log(paste("Restarting R session with the new app:", basename(app_file)))
    rstudioapi::restartSession(command = cmd)
    stopApp()
  })
  
  # ---- HMDB biospecimen / origin options, read from the selected CSV ----
  hmdb_levels_available <- debounce(reactive(read_hmdb_levels(input$file_path)), 500)
  
  # Selection waiting to be applied (e.g. from a loaded config) once the choices exist
  hmdb_pending <- reactiveValues(biospecimen = NULL, origin = NULL)
  
  set_hmdb_choices <- function(lv) {
    pick <- function(choices, wanted) {
      if (is.null(wanted)) return(choices)               # default: everything ticked
      sel <- choices[tolower(choices) %in% tolower(unlist(wanted))]
      if (length(sel)) sel else choices
    }
    updateCheckboxGroupInput(session, "hmdb_biospecimen", choices = lv$biospecimen,
                             selected = pick(lv$biospecimen, hmdb_pending$biospecimen), inline = TRUE)
    updateCheckboxGroupInput(session, "hmdb_origin", choices = lv$origin,
                             selected = pick(lv$origin, hmdb_pending$origin), inline = TRUE)
    hmdb_pending$biospecimen <- NULL
    hmdb_pending$origin <- NULL
  }
  
  observeEvent(hmdb_levels_available(), {
    lv <- hmdb_levels_available()
    if (is.null(lv)) {
      updateCheckboxGroupInput(session, "hmdb_biospecimen", choices = character(0), selected = character(0), inline = TRUE)
      updateCheckboxGroupInput(session, "hmdb_origin", choices = character(0), selected = character(0), inline = TRUE)
      if (nzchar(input$file_path) && file.exists(input$file_path)) {
        showNotification("HMDB CSV needs 'Biospecimen' and 'Origin' columns.", type = "warning")
        append_log("HMDB CSV has no 'Biospecimen'/'Origin' columns - no selection options available.")
      }
      return()
    }
    set_hmdb_choices(lv)
    append_log(paste0("HMDB CSV options - Biospecimen: ", paste(lv$biospecimen, collapse = ", "),
                      " | Origin: ", paste(lv$origin, collapse = ", ")))
  }, ignoreNULL = FALSE)
  
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
      refresh_file_list()  
    }
  })
  
  detect_polarity <- function(raw_dir) {
    if (is.null(raw_dir) || !nzchar(raw_dir) || !dir.exists(raw_dir)) return(NULL)
    files <- list.files(raw_dir, pattern = "\\.(mzXML|mzML)$", ignore.case = TRUE, full.names = FALSE)
    files_lc <- tolower(files)
    pos_hits <- sum(grepl("pos", files_lc) | grepl("_p_", files_lc) | grepl("_pos", files_lc))
    neg_hits <- sum(grepl("neg", files_lc) | grepl("_neg", files_lc) | grepl("_n_", files_lc))
    if (pos_hits > neg_hits && pos_hits > 0) return("Positive")
    if (neg_hits > pos_hits && neg_hits > 0) return("Negative")
    return(NULL)
  }
  
  observeEvent(input$raw_dir, {
    pol <- detect_polarity(input$raw_dir)
    if (!is.null(pol)) {
      
      isolate({
        if (is.null(input$.last_manual_ion) || input$.last_manual_ion != TRUE) {
          updateRadioButtons(session, "Ionization_method", selected = pol)
        }
      })
    }
  })
  
  observeEvent(input$Ionization_method, {
    session$sendCustomMessage("ion_manual", list(val = TRUE))
    
    session$userData$ion_manual <- TRUE
  }, ignoreInit = TRUE)
  
  make_config_list <- function() {
    list(
      Ionization_method = input$Ionization_method,
      data_dir = input$data_dir,
      data_dir_raw_files = input$raw_dir,
      file_path = input$file_path,
      output_folder_name = input$output_name,
      hmdb_biospecimen = input$hmdb_biospecimen,
      hmdb_origin = input$hmdb_origin,
      data_dir_DDA = input$data_dir_DDA,
      EXPORT_SPECTRA = input$EXPORT_SPECTRA,
      apply_smoothing = input$apply_smoothing,
      smoothing_n = input$smoothing_n,
      smoothing_p = input$smoothing_p,
      intensity_threshold = input$intensity_threshold,
      ppm_tolerance = input$ppm_tolerance,
      ppm_tolerance_annotation = input$ppm_tolerance_annotation,
      rt_min = input$rt_min,
      rt_max = input$rt_max,
      mz_min_limit = input$mz_min_limit,
      mz_max_limit = input$mz_max_limit
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
      updateRadioButtons(session, "Ionization_method", selected = cfg$Ionization_method %||% input$Ionization_method)
      updateTextInput(session, "data_dir", value = cfg$data_dir %||% input$data_dir)
      updateTextInput(session, "raw_dir", value = cfg$data_dir_raw_files %||% input$raw_dir)
      updateTextInput(session, "file_path", value = cfg$file_path %||% input$file_path)
      updateTextInput(session, "output_name", value = cfg$output_folder_name %||% "")
      # Selection is applied once the options for the (possibly new) CSV have been read
      hmdb_pending$biospecimen <- cfg$hmdb_biospecimen
      hmdb_pending$origin <- cfg$hmdb_origin
      if (identical(cfg$file_path %||% input$file_path, input$file_path)) {
        lv <- hmdb_levels_available()   # same CSV: options won't be re-read, so apply now
        if (!is.null(lv)) set_hmdb_choices(lv)
      }
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
      updateNumericInput(session, "mz_min_limit", value = cfg$mz_min_limit %||% input$mz_min_limit)
      updateNumericInput(session, "mz_max_limit", value = cfg$mz_max_limit %||% input$mz_max_limit)
      showNotification("Config loaded", type = "message")
      append_log(paste("Loaded config from", f$name))
    }, error = function(e) {
      showNotification(paste("Failed to load config:", e$message), type = "error")
    })
  })
  
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
      hmdb_biospecimen = paste(input$hmdb_biospecimen, collapse = ", "),
      hmdb_origin = paste(input$hmdb_origin, collapse = ", "),
      pipeline_status = if(rv$running) "RUNNING" else "IDLE",
      last_output_dir = rv$last_output_dir %||% "-"
    )
  })
  
  # Positive <-> Negative: switch to the sibling "<data_dir>_neg" folders if they exist.
  # (Runs only on a user change, never at startup, so configured paths are not overwritten.)
  observeEvent(input$Ionization_method, {
    sw <- swap_polarity_dirs(input$data_dir, input$raw_dir, input$Ionization_method)
    if (is.null(sw)) {
      append_log(paste0("Ionization mode set to ", input$Ionization_method,
                        " - no matching '_neg' / positive folder found, paths unchanged."))
      return()
    }
    updateTextInput(session, "data_dir", value = sw$data_dir)
    updateTextInput(session, "raw_dir", value = sw$raw_dir)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  
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
  
  observeEvent(input$open_output_folder, {
    outdir <- rv$last_output_dir
    if (is.null(outdir) || !dir.exists(outdir)) outdir <- input$data_dir
    if (!nzchar(outdir) || !dir.exists(outdir)) {
      showNotification("Output directory does not exist", type = "error"); return()
    }
    if (OS_TYPE == "windows") shell.exec(normalizePath(outdir))
    else if (OS_TYPE == "macos") system2("open", args = c(shQuote(outdir)))
    else system2("xdg-open", args = c(shQuote(outdir)))
  })
  
  validate_inputs <- function() {
    errs <- character()
    if (!nzchar(input$script_path) || !file.exists(input$script_path)) errs <- c(errs, "Pipeline script not found.")
    if (!nzchar(input$data_dir) || !dir.exists(input$data_dir)) errs <- c(errs, "Processing data_dir not found.")
    if (!nzchar(input$raw_dir) || !dir.exists(input$raw_dir)) errs <- c(errs, "Raw files directory not found.")
    if (!nzchar(input$file_path) || !file.exists(input$file_path)) errs <- c(errs, "HMDB CSV not found.")
    mzmls <- list.files(input$raw_dir, pattern = "\\.(mzXML|mzML)$", ignore.case = TRUE)
    if (length(mzmls) == 0) errs <- c(errs, "No mzML files found in raw_dir.")
    if (is.na(input$mz_min_limit) || is.na(input$mz_max_limit)) {
      errs <- c(errs, "mz_min_limit and mz_max_limit must be set.")
    } else if (input$mz_min_limit >= input$mz_max_limit) {
      errs <- c(errs, "mz_min_limit must be smaller than mz_max_limit.")
    }
    if (nzchar(input$file_path) && file.exists(input$file_path) && is.null(read_hmdb_levels(input$file_path))) {
      errs <- c(errs, "HMDB CSV needs 'Biospecimen' and 'Origin' columns.")
    } else {
      if (length(input$hmdb_biospecimen) == 0) errs <- c(errs, "Select at least one HMDB biospecimen.")
      if (length(input$hmdb_origin) == 0) errs <- c(errs, "Select at least one HMDB origin.")
    }
    if (length(errs)) return(paste(errs, collapse = "\n"))
    return(NULL)
  }
  
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
  
  observeEvent(input$run_btn, {
    if (rv$running) { showNotification("Pipeline already running", type = "warning"); return() }
    err <- validate_inputs()
    if (!is.null(err)) {
      showModal(modalDialog(title = "Validation failed", paste(err, collapse = "\n"), easyClose = TRUE, footer = modalButton("Close")))
      append_log(paste0("Validation failed: \n", err))
      return()
    }
    mzmls <- list.files(input$raw_dir, pattern = "\\.(mzXML|mzML)$", ignore.case = TRUE, full.names = TRUE)
    est_sec <- length(mzmls) * rv$est_seconds_per_mzml
    est_text <- paste0("Estimated time: ~", round(est_sec/60,1), " minutes for ", length(mzmls), " files.")
    nm <- clean_folder_name(input$output_name)
    out_text <- if (nzchar(nm)) {
      od <- make_output_dir(input$data_dir, nm)
      paste0("Output folder: ", od, if (dir.exists(od)) "  (already exists - existing files may be overwritten)" else "")
    } else {
      paste0("Output folder: ", file.path(sub("/+$", "", input$data_dir), "<date_time>"), "  (created at start)")
    }
    showModal(modalDialog(
      title = "Confirm Run",
      tagList(p("You're about to start the pipeline."), p(est_text), p(out_text)),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_run", "Start", class = "btn-primary")
      ), easyClose = TRUE
    ))
  })
  
  observeEvent(input$confirm_run, {
    removeModal()
    cfg <- make_config_list()
    
    out_dir <- make_output_dir(input$data_dir, input$output_name)
    ok <- dir.exists(out_dir) || dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok) {
      append_log(paste("Could not create output folder:", out_dir))
      showNotification(paste("Could not create output folder:", out_dir), type = "error")
      return()
    }
    subdirs <- stats::setNames(file.path(out_dir, OUTPUT_SUBFOLDERS), OUTPUT_SUBFOLDERS)
    for (d in subdirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    rv$last_output_dir  <- out_dir
    cfg$output_dir      <- out_dir                  # CSV results go here
    cfg$output_subdirs  <- as.list(subdirs)         # named list of subfolders inside output_dir
    cfg$mzml_export_dir <- input$data_dir           # exported mzML spectra go directly into the processing folder
    cfg$pipeline_script_path <- input$script_path   # recorded (with md5) in run_parameters.yml
    append_log(paste0("Output folder: ", out_dir))
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
    
    rv$r_process <- r_bg(
      function(script, config) {
        env <- new.env()
        list2env(config, envir = env)
        cat("===== Pipeline started =====\n")
        source(script, local = env, echo = TRUE)
      },
      args = list(script = script, config = cfg)
    )
    
    observe({
      req(rv$running)  
      invalidateLater(100, session)
      req(rv$r_process)
      
      if (rv$r_process$is_alive()) {
        out <- rv$r_process$read_output_lines()
        err <- rv$r_process$read_error_lines()
        if (length(out) && any(nzchar(out))) {
          for (ln in out[nzchar(out)]) {
            append_log(ln)
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
          
          if (is.na(rv$progress_total)) {
            mzmls <- list.files(input$raw_dir, pattern = "\\.(mzXML|mzML)$", ignore.case = TRUE, full.names = TRUE)
            if (length(mzmls) > 0) {
              rv$progress_total <- length(mzmls)
              processed <- sum(sapply(mzmls, function(p) any(grepl(basename(p), tail(rv$log, 200), fixed = TRUE))))
              rv$progress <- processed
            }
          }
        }
      } else {
        out <- rv$r_process$read_output_lines()
        err <- rv$r_process$read_error_lines()
        if (length(out) && any(nzchar(out))) for (ln in out[nzchar(out)]) append_log(ln)
        if (length(err) && any(nzchar(err))) for (ln in err[nzchar(err)]) append_log(paste0("ERR: ", ln))
        
        rv$running <- FALSE
        rv$r_process <- NULL
        rv$progress <- 0
        rv$progress_total <- NA
        
        append_log("Pipeline finished. Ready for next run.")
        showNotification("Pipeline completed successfully!", type = "message", duration = 5)
        if (!is.null(rv$last_output_dir)) updateTextInput(session, "csv_explorer_dir", value = rv$last_output_dir)
        
        refresh_file_list()
      }
    })
  }, ignoreInit = TRUE)
  
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
  
  output$download_log <- downloadHandler(
    filename = function() paste0("fia_run_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"),
    content = function(file) writeLines(rv$log, con = file)
  )
}

shinyApp(ui, server)