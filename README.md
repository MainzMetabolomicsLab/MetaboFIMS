# MetaboFIMS — Shiny App

A cross-platform R Shiny application for automated processing of **Flow Injection Mass Spectrometry (FIMS/FIA)** data in positive and negative ion mode. 
---

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Configuration](#configuration)
- [Pipeline Script Interface](#pipeline-script-interface)
- [Output & File Explorer](#output--file-explorer)
- [Troubleshooting](#troubleshooting)
- [Project Structure](#project-structure)

---

## Features

- **Positive / Negative ion mode** switching with automatic directory update
- **Polarity auto-detection** from raw file names (`pos`, `neg`, etc.)
- **Background process execution** — the pipeline runs in a separate R process, keeping the UI responsive
- **Live console log** with timestamps, scrollable output, and downloadable log file
- **Progress bar** driven by `PROGRESS: i/N` messages emitted by the pipeline script
- **File explorer** with CSV preview inside the app
- **Save / load YAML config profiles** for reproducible runs
- **Input validation** with a confirmation dialog and estimated runtime before each run
- **Graceful stop** — kill the background process at any time

---

## Requirements

### R packages

The app installs missing packages automatically on first launch. You will need an internet connection for the initial setup. The following are required:

**Bioconductor:**
`xcms`, `MSnbase`, `CAMERA`, `Rdisop`, `BiocParallel`, `Spectra`, `S4Vectors`, `mzR`

**CRAN:**
`R.utils`, `webchem`, `data.table`, `dplyr`, `tidyverse`, `ggplot2`, `stringdist`, `shiny`, `shinyFiles`, `DT`, `callr`, `yaml`

### System

| Item | Minimum |
|---|---|
| R | ≥ 4.3 |
| OS | Windows 10+, macOS 12+, Linux (Ubuntu 20.04+) |
| RAM | 8 GB (16 GB+ recommended for large datasets) |
| Input files | `.mzML` or `.mzXML` centroided mass spectrometry data |

---

## Installation

1. **Clone or download** this repository:
   ```bash
   git clone https://github.com/MainzMetabolomicsLab/MetaboFIMS.git
   cd fims-pipeline
   ```

2. **Open `app.R`** in RStudio.

3. **Run the app:**
   ```r
   shiny::runApp("app.R")
   ```
   On first launch, missing packages are installed automatically. This may take several minutes, depending on how many packages need to be installed first.

4. Alternatively, launch from the terminal:
   ```bash
   Rscript -e "shiny::runApp('app.R')"
   ```

---

## Usage

### 1. Select ionization mode

Use the **Positive / Negative** radio buttons at the top of the sidebar. Switching mode automatically updates the default directory paths. The app will also attempt to auto-detect polarity from your raw file names when you set the raw files directory, if files are only present in either of the directories.

### 2. Set file paths

| Field | Description |
|---|---|
| **R Script** | Path to processing pipeline `.R` script |
| **Processing data directory** | Output/working directory for the pipeline (`data_dir`) |
| **Raw files directory** | Folder containing your `.mzML` / `.mzXML` files |
| **HMDB annotation CSV** | Path to the HMDB metabolite database CSV used for annotation |

Use the **Browse** buttons next to each field to navigate your filesystem.

### 3. Set processing parameters

| Parameter | Description |
|---|---|
| `rt_min` / `rt_max` | Retention time window (seconds) for 1D spectra accumulation |
| Export Spectra | Export 1D spectra as `.mzML` files |
| Apply smoothing | Enable Savitzky–Golay smoothing before centroiding |
| Smoothing window (n) | Window size (odd integer, e.g. 7) |
| Smoothing polynomial (p) | Polynomial degree (e.g. 2) |
| Intensity threshold | Minimum peak intensity to retain |
| PPM tolerance | Mass accuracy window for peak picking |
| PPM tolerance (annotation) | Mass accuracy window for HMDB annotation |

### 4. Run the pipeline

Click **Run pipeline**. A confirmation dialog will appear showing the number of detected `.mzML` files and an estimated runtime. Click **Start** to launch the background process.

Monitor progress in the **Live Console / Log** panel. If your pipeline script emits lines like:

```
PROGRESS: 3/12
```

the progress bar will update automatically.

### 5. Stop the pipeline

Click **Stop (graceful)** to kill the background process at any time. The UI resets and is ready for the next run.

---

## Configuration

### Saving a profile

1. Enter a name in the **Profile name** field.
2. Click **Save config**.

A `.yml` file is written to the current working directory (`getwd()`).

### Loading a profile

Click **Load config (.yml)** and select a previously saved `.yml` file. All fields are populated automatically.

Config files are plain YAML and can be edited by hand:

```yaml
Ionization_method: Positive
data_dir: /Volumes/T7/Arbeit/FIA/processing
data_dir_raw_files: /Volumes/T7/Arbeit/FIA/processing/raw_mzML
file_path: /Volumes/T7/Arbeit/FIA/HMDB_BLOOD_ENDO_FOOD_DRUG_DRUGMET.csv
rt_min: 5
rt_max: 25
apply_smoothing: true
smoothing_n: 7
smoothing_p: 2
intensity_threshold: 5000
ppm_tolerance: 20
ppm_tolerance_annotation: 20
EXPORT_SPECTRA: false
```

---

## Pipeline Script Interface

The app sources your pipeline R script in a **clean environment** and injects all config values as variables. Your script can use them directly by name, for example:

```r
# These variables are available inside your pipeline script:
Ionization_method        # "Positive" or "Negative"
data_dir                 # processing output directory
data_dir_raw_files       # raw mzML directory
file_path                # HMDB CSV path
rt_min                   # numeric
rt_max                   # numeric
apply_smoothing          # logical
smoothing_n              # integer
smoothing_p              # integer
intensity_threshold      # numeric
ppm_tolerance            # numeric
ppm_tolerance_annotation # numeric
EXPORT_SPECTRA           # logical
```


## Output & File Explorer

Use the **Data Directory Explorer** panel to browse the output folder:

- Set the directory path manually or use the **Browse** button.
- Click **Refresh file list** to update after a run.
- Click **Open output folder** to open the directory in your system file manager.
- Click any `.csv` row in the file table to preview its contents directly in the app.

Log files are written to a `logs/` subfolder in the app's working directory, named by date (e.g. `logs/2025-06-01_fia_run_log.txt`). Use **Download log** to save the current session log.

---

## Troubleshooting

**App won't start / package errors**
Run the package installation block manually in an R console and check for errors, particularly for Bioconductor packages which require a matching `BiocManager` version.

**"Script not found" error**
Make sure the path in the **R Script** field points to an existing `.R` file. Use the Browse button to avoid typos.

**No mzML files detected**
Confirm that your raw files directory contains files with a `.mzML` or `.mzXML` extension (case-insensitive). The validation check runs before the pipeline starts.

**Pipeline appears to hang**
Check the Live Console for the last output line. If the process is stuck, use **Stop** and check whether the pipeline script itself has an issue outside the app.

**Progress bar stays at 0%**
The progress bar only advances if your pipeline script emits `PROGRESS: i/N` lines. Without these, the app falls back to counting how many raw filenames have appeared in the log.

---

## Project Structure

```
fims-pipeline/
├── app.R                         # Main Shiny application
├── logs/                         # Auto-created; daily run logs
├── *.yml                         # Saved config profiles (created at runtime)
└── README.md                     # This file
```

The pipeline script itself (e.g. `FIA_XCMS_isotope_neutral_golden_rulez_incl_smoothing_only_RT_pos_neg_combined.R`) is referenced by path and is not bundled with the app.

---

## License

MIT — see `LICENSE` for details.
