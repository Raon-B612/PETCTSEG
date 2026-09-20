# libraries
suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(readxl)
  library(tidyverse)
  })

# paths and config
home_dir    <- path(Sys.getenv("USERPROFILE"), "Documents", "PETCTSEG")
script_dir  <- path(home_dir, "scripts")
prep_dir    <- path(script_dir, "prep")
measure_dir <- path(home_dir, "measurements/current")
summary_dir <- path(home_dir, "summaries")

# Load VOI labels
voi    <- read_xlsx(path(script_dir, "voi_labels_v3.xlsx"))
voi_dt <- as.data.table(voi)[, .(module, intensities, label)]
setkey(voi_dt, module, intensities)

cols2select <- c("suv_max", "suv_mean", "hu_mean", "volume", "area")

hc_modules  <- unique(voi$module)
modules     <- c("tseg_total", "tseg_tissue_types", "tseg_vertebrae_body", "tseg_vertebrae_pp_refined", "moose_body_composition")
tmj_modules <- c("tseg_head_glands_cavities", "tseg_head_muscles", "tseg_headneck_bones_vessels", "tseg_headneck_muscles", "tseg_craniofacial_structures")
ent_modules <- c(tmj_modules, modules)

# Function to load all utility and core scripts in order [Restored]
load_utils <- function() {
  utils <- c("core_logic.R", "utils_file_system.R", 'utils_database.R')
  for (u in utils) {
    u_path <- path(script_dir, 'utils', u)
    if (file.exists(u_path)) source(u_path) else stop('[ERROR] Loading libaries incomplete..')
  }
}

# Execute loader [Restored]
load_utils()

# Initialize directories (Must be after load_utils to use get_filtered_dirs)
gz_dir  <- "C:/Temp/PETCTSRC"
gz_dir1 <- get_filtered_dirs("D:")
gz_dir2 <- get_filtered_dirs("E:")
gz_all  <- c(gz_dir1, gz_dir2) |> sort()

hc_dir  <- gz_all |> keep(\(x) str_detect(x, "HC[BMPVX]$"))
tmj_dir <- gz_all |> keep(\(x) str_detect(x, "TMJ[BMPVX]$"))

# doi configuration (Preserved commented-out lines)
# doi <- c("CA", "C1", "C2", "CM", "EN", "EY", "G1", "G4", "G5", "G6", "G7", "G8", "G9", "HC", "M1", "M2", "NU", "OB", "PSMA2", "PSMA3", "PSMA", "TMJ", "TR")