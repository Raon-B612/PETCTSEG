# Create parent directory for a file
mkdir_file <- function(file) {
  if (!dir.exists(dirname(file))) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  }
}

# Create a directory
mkdir <- function(files_or_folders) {
  dl <- path(files_or_folders)
  newdir <- c(dl[is_dir(dl)], dirname(dl[!is_dir(dl)])) |> unique()
  dir_create(newdir)
  cat(paste0("Successfully created ", newdir, "\n"), sep = "")
}

mkdir_dir <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
}


# Ensure parent directory exists
ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

# Retrieve directories ending with department suffixes
get_filtered_dirs <- function(drive) {
  list.dirs(drive, full.names = TRUE, recursive = FALSE) %>%
    (\(x) x[str_detect(x, "[BMPVX]$")])() %>%
    normalizePath(winslash = "/", mustWork = FALSE)
}

# Retrieve directories ending with department suffixes
get_filtered_dirs <- function(drive) {
  list.dirs(drive, full.names = TRUE, recursive = FALSE) %>%
    (\(x) x[str_detect(x, "[BMPVX]$")])() %>%
    normalizePath(winslash = "/", mustWork = FALSE)
}

# Import PACS Excel list
import_pacs_db <- function(src) {
  null_file <- list.files("D:/LIST/", full.names = TRUE, recursive = FALSE, pattern = src) %>%
    sort(decreasing = TRUE) %>% .[1]
  # Strictly preserved Korean column names
  choose.files(default = null_file) %>%
    read_xls() %>%
    rename(
      accn  = `ACCESSION NUM`, pname = `ENG NAME`, pid = 등록번호,
      pdate = 검사일시, exam = 검사명, dept = 처방과, room = 촬영실_상세
    ) %>%
    mutate(pdate = str_extract(pdate, "\\d{4}-\\d{2}-\\d{2}") %>% ymd()) %>%
    distinct() %>% filter(!is.na(pid))
}

# Write CSV atomically using temp file
write_atomic_csv <- function(dt, out) {
  ensure_parent_dir(out)
  tmp <- paste0(out, ".tmp")
  if (file.exists(tmp)) file.remove(tmp)
  fwrite(dt, tmp)
  ok <- file.rename(tmp, out)
  if (!ok) {
    if (file.exists(tmp)) file.remove(tmp)
    stop(sprintf("atomic write failed: %s", out))
  }
  invisible(out)
}

# Global environment to cache file list for memory-based lookup
.seg_cache <- new.env()

# Recursive directory listing optimized for HDD speed
dir_ls_recursive <- function(path) {
  # Optimization: Replaced manual R recursion with optimized system calls using fs::dir_ls(recurse = TRUE)
  items <- try(dir_ls(path, recurse = TRUE, type = "any", fail = FALSE), silent = TRUE)
  if (inherits(items, "try-error")) return(character(0))
  
  # Exclude system-protected directories
  items <- items[!grepl("\\$RECYCLE.BIN|System Volume Information", items)]
  
  # Filter directories that match the subject date pattern (8 digits)
  dirs <- items[is_dir(items) & grepl("\\d{8}$", items)]
  as.character(dirs)
}

# Parse subject directory metadata with bulk I/O optimization
dir_parsing <- function(source_dir) {
  # Optimization: Scan all files in the source directory once and store in memory cache
  .seg_cache$all_files <- as.character(dir_ls(source_dir, recurse = TRUE, type = "file", fail = FALSE))
  
  dir_ls_recursive(source_dir) |> tibble(dir = _) |>
    mutate(
      subjid_long = basename(dir), 
      subjid = str_remove(subjid_long, "_\\d{8}"), 
      prefix = str_remove(subjid, "\\d{10,13}"), 
      pdate_flat = str_extract(subjid_long, "_\\d{8}") |> str_remove("_"), 
      pdate = ymd(pdate_flat), 
      accn = str_extract(subjid, "\\d{10,13}")
    ) |>
    mutate(
      ct = file.path(dir, paste0("CT_", subjid, ".nii.gz")), 
      suv = file.path(dir, paste0("SUV_", subjid, ".nii.gz")), 
      pet = case_when(
        str_detect(prefix, "TMJ") ~ file.path(dir, paste0("SPECT_", subjid, ".nii.gz")), 
        TRUE ~ file.path(dir, paste0("PET_", subjid, ".nii.gz"))
      )
    )
}

# Generate structured analysis dataframe with memory-based existence checks
generate_expanded_df <- function(parsed_df) {
  # # prevent using generate_expanded_df() without dir_parsing()
  # if (is.null(.seg_cache$all_files)) {
  #   stop("[!!] .seg_cache$all_files가 비어있습니다 — dir_parsing()을 먼저 호출하세요.")
  # }
  
  df_hc    <- parsed_df |> filter( str_detect(subjid, '^HC[BMPVX]'))
  df_tmj   <- parsed_df |> filter( str_detect(subjid, '^TMJ[BMPVX]'))
  df_ent   <- parsed_df |> filter( str_detect(subjid, '^EN[BMPVX]'))
  df_nonhc <- parsed_df |> filter(!str_detect(subjid, '^(HC[BMPVX]|TMJ[BMPVX]|EN[BMPVX])'))
  
  bind_rows(
    crossing(df_tmj, tmj_modules) |> rename(module = tmj_modules), 
    crossing(df_ent, ent_modules) |> rename(module = ent_modules), 
    crossing(df_hc, hc_modules) |> rename(module = hc_modules),
    crossing(df_nonhc, modules) |> rename(module = modules)
  ) |>
    mutate(
      seg = file.path(dir, paste0(module, "_", subjid, ".nii.gz")), 
      csv = file.path(measure_dir, prefix, subjid_long, paste0(module, "_", subjid, ".csv")),
      # Optimization: Replaced file_exists() with in-memory set matching to avoid HDD latency
      ct_exists = ct %in% .seg_cache$all_files, 
      pet_exists = pet %in% .seg_cache$all_files, 
      suv_exists = suv %in% .seg_cache$all_files, 
      seg_missing = !(seg %in% .seg_cache$all_files), 
      # Keep file_exists for CSV if they are on a separate drive (SSD/OneDrive)
      csv_missing = !file_exists(csv), 
      img_missing = !(ct_exists & pet_exists & suv_exists)
    )
}

# Rename MOOSE model outputs
rename_moose <- function(subjdir) {
  mkdir_dir("C:/Temp/todel/moosez")
  c(subjdir) |>
    map(~ dir_ls(.x, regexp = "clin_CT_.*_segmentation_CT.*\\.nii.gz", recurse = TRUE)) |>
    flatten_chr() |>
    tibble(infile = _) |>
    mutate(
      moose_dir = str_extract(infile, ".*/[^/]+_\\d{8}"),
      outfile = path(moose_dir, str_replace(basename(infile), "clin_CT", "moose") |> str_remove("_segmentation_CT"))
    ) |>
    pwalk(function(infile, moose_dir, outfile) {
      file.copy(infile, "C:/Temp/todel/moosez", overwrite = TRUE)
      file.rename(infile, outfile)
      cat("\n", infile, "-->", outfile)
    })
  cat('\n')
}

# Rename LIFEx SUV outputs
rename_lx <- function(src_dir = gz_dir) {
  suv_fl <- list.files(src_dir, "SUVbw.*it[0-9]", full.names = TRUE, recursive = TRUE)
  for (i in suv_fl) {
    sub_long <- str_extract(i, "[A-Z0-9]+[BMPVX]\\d{10,13}_\\d{8}")
    sub_id   <- str_extract(sub_long, "[A-Z0-9]+[BMPVX]\\d{10,13}")
    prefix   <- str_remove(sub_id, "\\d{10,13}$")
    outfile  <- file.path(src_dir, prefix, sub_long, paste0("SUV_", sub_id, ".nii.gz"))
    mkdir_file(outfile)
    file.copy(i, outfile, overwrite = TRUE)
    print(outfile)
  }
}

# Move deprecated or raw files to cleanup root
move_deprecated <- function(src_dir = gz_dir, target_root = "C:/Temp/todel") {
  # src_dir <- gz_dir
  files_to_move <- list.files(src_dir, full.names = TRUE, recursive = TRUE) %>%
    keep(~ str_detect(.x, "ROI|Eq|SUVbw|json$|[a-z]\\.nii\\.gz")) %>%
    c(dir_ls(src_dir, regexp = "3D|moosez-", type = "directory", recurse = TRUE))
  
  walk(files_to_move, function(path0) {
    rel  <- path_rel(path0, start = src_dir)
    dest <- file.path(target_root, rel)
    dir_create(dirname(dest))
    if (dir_exists(path0)) { dir_copy(path0, dest, overwrite = TRUE); dir_delete(path0) }
    else { file_copy(path0, dest, overwrite = TRUE); file_delete(path0) }
  })
}

# Distribute processed files to permanent drive locations
distribute_processed_files <- function(processed_dir = gz_dir) {
  # processed_dir <- gz_dir
  cat('\n')
  
  df_to_distribute <- dir_parsing(processed_dir) |>
    select(dir, prefix, subjid_long) |>
    # (\(x) x[!duplicated(x),])() |>
    mutate(
      output_root = case_when(
        str_detect(prefix, "^(C1|C2|DC|EN|G5|G6|G7|GS|H1|H3|HC|HM|M1|M2|M3|M5|M7|M8|M9|OB|PSMA|PSMA2|PSMA3|TMJ|UR)[BMPVX]$") ~ "D:",
        str_detect(prefix, "^(AN|CA|CM|EY|G1|G2|G4|G8|G9|MU|N1|NM|NS|NU|O1|OB|OS|P1|PS|PY|RD|RH|SK|TR)[BMPVX]$") ~ "E:",
        TRUE ~ "C:/Temp/Unsorted"
      ),
      output_dir = file.path(output_root, prefix, subjid_long)
    ) |>
    arrange(output_root)
  
  for (i in 1:nrow(df_to_distribute)) {
    # i <- 1
    x <- df_to_distribute[i,] |>
      mutate(
        arx_dir = str_replace(output_dir, 'D:/|E:/', 'Z:/PETCTSEG/')
      )
    
    if (str_detect(x$dir, 'D:')|str_detect(x$dir, 'E:')) stop('!= Only files in gz_dir are allowed..')
    
    files_to_copy <- list.files(x$dir, full.names = TRUE)
    
    if (any(file.size(files_to_copy)==0)) stop('!= Zero-size file detected..') 
    if (!str_detect(x$prefix, 'HC[BMPVX]')) files_to_copy <- files_to_copy[!str_detect(basename(files_to_copy), 'PET_')] 
    
    mkdir_dir(x$output_dir)
    file.copy(files_to_copy, x$output_dir, overwrite = TRUE)
    cat('..', x$dir, "-->", x$output_dir)
    
    if (dir.exists('Z:/PETCTSEG/')) {
      mkdir_dir(x$arx_dir)
      file.copy(files_to_copy, x$arx_dir, overwrite = TRUE)
      cat(" -->", x$arx_dir)
    }
    
    cat("\n")
  }
}

# Simple vectorized approach to reorder rows from the center
reorder_from_center <- function(df) {
  n <- nrow(df)
  # Calculate the center index: floor((n + 1) / 2)
  center <- (n + 1) %/% 2
  
  # Sort indices based on absolute distance from the center, 
  # then use the index itself as a tie-breaker to ensure 'center-1' comes before 'center+1'
  df_reordered <- df[order(abs(1:n - center), 1:n), ]
  
  # Reset row names for a clean output
  rownames(df_reordered) <- NULL
  return(df_reordered)
}
