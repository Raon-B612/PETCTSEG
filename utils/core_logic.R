# Build DICOM to NIfTI mapping
read_dicom_dir <- function(src_dir = 'C:/Temp/PETCTSRC/', dcm_dir = dcm_path) {
  message('[INFO] ', path(dcm_dir), ' --> ', path(src_dir))
  
  dcm_df <- dir_ls(dcm_dir, type = "directory",  recurse = TRUE, regexp = "(CT|PT)$") |>
    tibble(dir = _) |>
    mutate(
      subjid_long = str_extract(dir, patterns$subjid_long),
      mod = case_when(
        str_detect(dir, "CT$")    ~ 'ctdcm',
        str_detect(dir, "PT$")    ~ 'petdcm',
        str_detect(dir, "SPECT$") ~ 'spectdcm'
      )
    ) |>
    pivot_wider(
      names_from  = mod,
      values_from = dir
    ) |>
    left_join(seg_db |> select(subjid_long, qdate), by = "subjid_long") |>
    filter(is.na(qdate)) |>
    mutate(
      subjid      = str_extract(subjid_long, patterns$subjid),
      pdate       = str_extract(subjid_long, '\\d{8}$') |> ymd(),
      prefix      = str_extract(subjid_long, patterns$prefix),
      output_dir  = path(src_dir, prefix, subjid_long),
      ct          = path(output_dir, paste0('CT_', subjid, '.nii.gz')),
      suv         = path(output_dir, paste0('SUV_', subjid, '.nii.gz')),
    )
  
  write_csv(dcm_df, path(prep_dir, "petctdcm.csv"))
  message('[INFO] petctdcm.csv: saved.\n')
  dcm_df
}

# Python script execution for CT NIfTI conversion
dcm2nii_ct <- function(dicom_df=dcm_df, skip = TRUE) {
  nii_df <- if (isTRUE(skip)) filter(dicom_df, !file.exists(ct)) else dicom_df
  write_csv(nii_df, path(prep_dir, "dcm2nii_prep.csv"))
  
  py_script <- path(script_dir, "conversion", "convert_ct_to_nifti.py")
  cmd_str   <- sprintf('start "DICOM2NII" cmd /k python -u "%s"', py_script)
  shell(cmd_str, wait = FALSE)
}

# LIFEx execution for PET SUV conversion
dcm2nii_pet <- function(dicom_df = dcm_df) {
  lx_df <- filter(dicom_df, !file.exists(suv))
  
  sx <- c("LIFEx.Script = Main", "LIFEx.Script.Version = 25.06.1", "")
  sx <- c(sx, unlist(map2(seq_len(nrow(lx_df)) - 1, lx_df$output_dir, \(i, outdir) {
    pfx <- paste0("LIFEx.Patient", i, ".Series0")
    c(paste0(pfx, "=", lx_df$petdcm[i + 1]), paste0(pfx, ".Operation0=SUVbw"), 
      paste0(pfx, ".Operation0.Output.Directory=", str_replace_all(outdir, "\\\\", "/"), "/"), 
      paste0(pfx, ".Operation0=Save nii float32"), "")
  })))
  
  write_lines(sx, path(prep_dir, "convert_suv_to_nifti_lifex.txt"))
  
  suppressMessages(
    shell('"C:\\Users\\Heartwork101\\AppData\\Local\\LIFEx-25.06.1\\LIFEx-25.06.1.exe"', wait = FALSE)
    )
}

# SUV image resampling to CT space
resample_suv_to_ct <- function(source_dir, overwrite = 'FALSE') {
  rename_lx(source_dir)
  move_deprecated(source_dir)
  
  gz_df <- dir_parsing(source_dir)
  
  incomplete <- gz_df |> filter(!file.exists(ct)|!file.exists(suv))
  if (nrow(incomplete) != 0) stop(cat('[ERROR] missing NifTI files in', unlist(select(incomplete, dir)), '..'))
  
  if (overwrite) {
    rs_df <- gz_df
  } else {
    rs_df <- gz_df |> filter(file.exists(ct), file.exists(suv), !file.exists(pet))
  }
  write_csv(rs_df, path(prep_dir, "resample_suv_prep.csv"))
  
  if (nrow(rs_df) > 0) {
    is_hdd <- str_detect(as.character(rs_df$ct[1]), '^[DEde]:')
    hdd_flag <- if (is_hdd) "--hdd" else ""
    
    py_script <- path(script_dir, "preprocessing", "resample_suv_to_ct.py")
    cmd_str   <- sprintf('start "RESAMPLE_PET_TO_CT" cmd /k python -u "%s" %s', py_script, hdd_flag)
    
    shell(cmd_str, wait = FALSE)
  }
}

# Run TotalSegmentator model
run_totalsegmentator <- function() {
  py_script <- path(script_dir, "segmentation", "segment_total.py")
  cmd_str   <- sprintf(
    'start /wait "TOTALSEGMENTATOR" cmd /k ""C:\\Users\\Heartwork101\\miniconda3\\Scripts\\activate.bat" "C:\\Users\\Heartwork101\\miniconda3" && conda activate totalseg_env && python -u "%s""',
    py_script
  )
  shell(cmd_str, wait = TRUE)
}

# Run MOOSE model
run_moose <- function() {
  py_script <- path(script_dir, "segmentation", "segment_moose.py")
  cmd_str   <- sprintf(
    'start /wait "MOOSE_SEG" cmd /k ""C:/Users/Heartwork101/moose/Scripts/activate.bat" && python -u "%s""',
    py_script
  )
  shell(cmd_str, wait = TRUE)
}

# Prepare CSV files for segmentation
prepare_segmentation <- function(source_dir = gz_dir, task = 'all') {
  # source_dir <- hc_dir; task <- 'all'
  rename_moose(source_dir)
  
  gz_df <- dir_parsing(source_dir) |> generate_expanded_df()
  
  incomplete <- gz_df |> dplyr::filter(!file.exists(ct))
  if (nrow(incomplete) != 0) stop(sprintf("[!!] CT image missing in %s ..", paste(incomplete$dir, collapse = ", ")))
  
  seg_df <- filter(gz_df, seg_missing) |> mutate(output_dir = dir) |> arrange(ct)
  
  ts_df <- seg_df |> 
    filter(str_detect(module, "tseg")) |> 
    select(ct, seg, module) |> 
    arrange(ct)
  write_csv(ts_df, path(prep_dir, 'tseg_prep.csv')) # all segmentation queued
  
  
  ts_df <- switch(
    task,
    total = ts_df[str_detect(ts_df$module, "total"), ],
    other = ts_df[!str_detect(ts_df$module, "total"), ],
    ts_df
  )
  write_csv(ts_df |> arrange(desc(ct)), path(prep_dir, 'tseg_prep_emma.csv')) # queded for emma
  
  # convert path for wsl env
  ts_df |>
    mutate(
      ct  = ct  |> str_replace('^C\\:', '/mnt/c') |> str_replace('^D\\:', '/mnt/d') |> str_replace('^E\\:', '/mnt/e'),
      seg = seg |> str_replace('^C\\:', '/mnt/c') |> str_replace('^D\\:', '/mnt/d') |> str_replace('^E\\:', '/mnt/e')
    ) |>
    arrange(desc(ct)) |>
    write_csv(path(prep_dir, 'tseg_prep_wsl.csv'))
  
  # create moose segmentation list
  ms_df <- seg_df |> 
    filter(str_detect(module, "moose")) |> 
    select(ct, module, output_dir) |> 
    mutate(outfile = path(output_dir, paste0('clin_CT_', str_remove(module, 'moose_'), '_segmentation_', basename(ct)))) |>
    arrange(ct) |>
    filter(!file.exists(outfile))
  write_csv(ms_df, path(prep_dir, "moose_prep.csv"))
  write_csv(ms_df |> arrange(desc(ct)), path(prep_dir, 'moose_prep_emma.csv')) # queded for emma
  
  # convert path for wsl env
  ms_df |>
    mutate(
      ct  = ct  |> str_replace('^C\\:', '/mnt/c') |> str_replace('^D\\:', '/mnt/d') |> str_replace('^E\\:', '/mnt/e'),
      output_dir = output_dir |> str_replace('^C\\:', '/mnt/c') |> str_replace('^D\\:', '/mnt/d') |> str_replace('^E\\:', '/mnt/e'),
      outfile = outfile |> str_replace('^C\\:', '/mnt/c') |> str_replace('^D\\:', '/mnt/d') |> str_replace('^E\\:', '/mnt/e')
    ) |>
    arrange(desc(ct)) |>
    write_csv(path(prep_dir, 'moose_prep_wsl.csv'))
  
  cat('\n[+] Saving segmentation lists ... done.')
}

# High-level segmentation wrapper
segment_ct <- function(source_dir) {
  run_totalsegmentator()
  run_moose()
  rename_moose(source_dir)
}

# Run ROI measurement
measure_roi <- function(source_dir = gz_dir, overwrite = FALSE) {
  # source_dir <- 'D:/PSMAV/'
  df <- dir_parsing(source_dir) |> generate_expanded_df()
  
  incomplete <- df |> filter(img_missing|seg_missing)
  if (nrow(incomplete) != 0) stop(cat('>> missing image or segmentation files in', unlist(select(incomplete, dir)), '..'))
  
  if (overwrite == FALSE) df <- filter(df, csv_missing)
    
  measure_df <- df |> 
    filter(!img_missing&!seg_missing) |> 
    mutate(output_dir = dirname(csv)) |>
    select(subjid, ct, pet, seg, module, output_dir, csv)
  write_csv(measure_df, path(prep_dir, 'measure_prep.csv'))
  
  is_hdd <- str_detect(as.character(measure_df$ct[1]), '^[DEde]:')
  hdd_flag <- if (is_hdd) "--hdd" else ""
  
  py_script <- path(script_dir, "measurement", "measure_roi.py")
  cmd_str <- sprintf('start "MEASURE_ROI" cmd /k python -u "%s" %s', py_script, hdd_flag)
  shell(cmd_str, wait = FALSE)
  # message("\n[INFO] ROI measurement completed: ", unlist(source_dir), '.')
}

# Helper: Empty QC data table
empty_qc_dt <- function() { data.table(dir = character(), subjid_long = character(), subjid = character(), prefix = character()) }

# Robust CSV reader
# [Modified] Enhanced intensities mapping while preserving original filtering
safe_read_measure <- function(file, cols2select) {
  info <- file.info(file)
  if (!file.exists(file)) return(list(ok = FALSE, reason = "file_not_found", dt = NULL))
  if (is.na(info$size) || info$size == 0) return(list(ok = FALSE, reason = "empty_file", dt = NULL))
  tryCatch({
    probe <- fread(file, nrows = 0, showProgress = FALSE)
    if (!length(names(probe))) return(list(ok = FALSE, reason = "no_columns", dt = NULL))
    dt <- fread(file, showProgress = FALSE, na.strings = c("", "NA"))
    if (!nrow(dt)) return(list(ok = FALSE, reason = "no_rows", dt = NULL))
    
    # Mapping logic for intensities
    if (!"intensities" %in% names(dt)) { 
      if ("label" %in% names(dt)) {
        setnames(dt, "label", "intensities")
      } else {
        if (ncol(dt) < 1L) return(list(ok = FALSE, reason = "cannot_assign_intensities", dt = NULL))
        setnames(dt, 1L, "intensities") 
      }
    }
    
    dt <- dt[!is.na(intensities) & intensities != 0]; if (!nrow(dt)) return(list(ok = FALSE, reason = "filtered_empty", dt = NULL))
    suppressWarnings(set(dt, j = "intensities", value = as.integer(dt[["intensities"]]))); miss <- setdiff(cols2select, names(dt))
    if (length(miss)) { for (m in miss) set(dt, j = m, value = NA_real_) }
    for (nm in cols2select) { suppressWarnings(set(dt, j = nm, value = as.numeric(dt[[nm]]))) }
    list(ok = TRUE, reason = NA_character_, dt = dt)
  }, error = function(e) { list(ok = FALSE, reason = conditionMessage(e), dt = NULL) })
}

# Collect subject directories
get_subjects_dt <- function(dept_dir) {
  dirs <- fs::dir_ls(dept_dir, type = "directory", recurse = FALSE, regexp = "_\\d{8}$")
  if (!length(dirs)) return(data.table())
  dt <- data.table(dir = as.character(dirs))
  dt[, subjid_long := basename(dir)][, subjid := str_remove(subjid_long, "_\\d{8}")][, prefix := str_remove(subjid, "\\d{10,13}")][]
}

# Map subjects to analysis modules
build_subject_module_map <- function(df, modules, tmj_modules, ent_modules, hc_modules) {
  if (!nrow(df)) return(data.table()); out <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    subjdir <- df$dir[i]
    modset <- if (str_detect(subjdir, "HC[BMPVX]")) hc_modules else if (str_detect(subjdir, "TMJ[BMPVX]")) tmj_modules else if (str_detect(subjdir, "EN[BMPVX]")) ent_modules else modules
    tmp <- data.table(dir = df$dir[i], subjid_long = df$subjid_long[i], subjid = df$subjid[i], prefix = df$prefix[i], module = modset)
    tmp[, csv := path(measure_dir, prefix, subjid_long, paste0(module, "_", subjid, ".csv"))][, csv_missing := !file.exists(csv)]; out[[i]] <- tmp
  }
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

# Aggregate results and perform QC
summarize_subject_measures <- function(overwrite = FALSE) {
  depts <- fs::dir_ls(measure_dir, type = "directory", recurse = FALSE, regexp = "[A-Z0-9]+[BMPVX]$")
  qc_fail_file <- path(script_dir, "logs", "qc_fail.csv"); qc_fail_cols <- c("dir", "prefix", "subjid_long", "subjid")
  qc_fail_dt <- if (file.exists(qc_fail_file)) tryCatch(fread(qc_fail_file), error = function(e) empty_qc_dt()) else empty_qc_dt()
  
  for (dept_to_summarize in depts) {
    df <- get_subjects_dt(dept_to_summarize); if (!nrow(df)) next
    csv_df <- build_subject_module_map(df, modules, tmj_modules, ent_modules, hc_modules)
    if (any(csv_df$csv_missing)) { 
      cat("\n[ERROR] Missing CSVs in:", basename(dept_to_summarize))
      filter(csv_df, csv_missing==TRUE) |> pull(subjid) |> unique() |> print()
      next 
      }
    
    csv_split <- split(csv_df, by = "dir"); cat("Summarizing", basename(dept_to_summarize), "...")
    dept_qc_accum <- vector("list", 0L)
    
    for (subjdir in names(csv_split)) {
      subj_csv_df <- csv_split[[subjdir]]; subjid <- str_extract(subjdir, "[A-Z0-9]+\\d{10,13}"); dept <- str_extract(dirname(subjdir), "[A-Z0-9]+[BMPVX]$")
      out <- path(summary_dir, dept, paste0("summary_", subjid, ".csv"))
      if (file.exists(out) && !isTRUE(overwrite)) next
      
      out_list <- vector("list", nrow(subj_csv_df)); ok_subject <- TRUE
      for (i in seq_len(nrow(subj_csv_df))) {
        res <- safe_read_measure(subj_csv_df$csv[i], cols2select); if (!isTRUE(res$ok)) { ok_subject <- FALSE; break }
        dt <- res$dt; set(dt, j = "module", value = subj_csv_df$module[i]); out_list[[i]] <- dt[, c("module", "intensities", intersect(cols2select, names(dt))), with = FALSE]
      }
      
      if (!ok_subject || !length(Filter(Negate(is.null), out_list))) { cat("\n[SKIP] Failed to read or empty:", subjid); next }
      
      summary_df <- voi_dt[rbindlist(out_list, fill = TRUE), on = .(module, intensities)]
      
      # [Modified] Final QC with TMJ Bypass logic
      is_tmj <- str_detect(subjid, "TMJ[BMPVX]")
      t1_missing <- is.na(summary_df[tolower(label) == "vertebrae_t1", suv_mean][1])
      
      if (!is_tmj && t1_missing) { 
        cat("\n[ERROR] QC FAIL - T1 missing for:", subjid)
        dept_qc_accum[[length(dept_qc_accum) + 1L]] <- unique(subj_csv_df[, ..qc_fail_cols], by = "dir"); next 
      }
      
      if (nrow(summary_df) > 0) {
        write_atomic_csv(summary_df, out)
      } else {
        cat("\n[WARN] Empty result after join for:", subjid)
      }
    }
    
    # [Added] Save accumulated QC failures to file
    if (length(dept_qc_accum) > 0) {
      new_fails <- rbindlist(dept_qc_accum)
      qc_fail_dt <- unique(rbind(qc_fail_dt, new_fails), by = "dir")
      if (!dir.exists(dirname(qc_fail_file))) dir.create(dirname(qc_fail_file), recursive = TRUE)
      fwrite(qc_fail_dt, qc_fail_file)
    }
  }
  cat('\n[+] Summaries are successfully saved.')
}
