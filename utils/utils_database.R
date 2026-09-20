# Retrieve directories ending with department suffixes
get_filtered_dirs <- function(drive) {
  list.dirs(drive, full.names = TRUE, recursive = FALSE) %>%
    (\(x) x[str_detect(x, "[BMPVX]$")])() %>%
    normalizePath(winslash = "/", mustWork = FALSE)
}

# Import PACS Excel list
import_pacs_db <- function(src) {
  null_file <- list.files("C:/Temp/", full.names = TRUE, recursive = FALSE, pattern = src) %>%
    sort(decreasing = TRUE) %>% .[1]
  # Strictly preserved Korean column names
  choose.files(default = null_file) %>%
    read_xls() %>%
    rename(
      accn  = `ACCESSION NUM`, pname = `ENG NAME`, pid = 등록번호,
      pdate = 검사일시, exam = 검사명, dept = 처방과, room = 촬영실_상세
    ) %>%
    mutate(
      pdate = str_extract(pdate, "\\d{4}-\\d{2}-\\d{2}") %>% ymd()
      ) %>%
    distinct() %>% 
    filter(!is.na(pid))
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

update_seg_db <- function() {
  # 7. 언더스코어 포함해 패턴을 더 특정하게
  all_dirs <- dir_ls(gz_all, recurse = FALSE, type = 'directory', regexp = '_\\d{8}$')
  
  # 3. 변수명 구분: df -> df_subjects
  df_subjects <- tibble(subjdir = all_dirs) |>
    mutate(
      subjid_long = basename(subjdir),
      subjid = str_remove(subjid_long, "_\\d{8}")
    )
  
  
  df_hc    <- df_subjects |> filter( str_detect(subjid, '^HC[BMPVX]'))
  df_tmj   <- df_subjects |> filter( str_detect(subjid, '^TMJ[BMPVX]'))
  df_ent   <- df_subjects |> filter( str_detect(subjid, '^EN[BMPVX]'))
  df_nonhc <- df_subjects |> filter(!str_detect(subjid, '^(HC[BMPVX]|TMJ[BMPVX]|EN[BMPVX])'))
  
  # 6. 카테고리가 상호배타적으로 df2를 모두 커버하는지 sanity check
  stopifnot(
    nrow(df_hc) + nrow(df_tmj) + nrow(df_ent) + nrow(df_nonhc) == nrow(df_subjects)
  )
  
  # 4. 그룹별 module 매핑을 딕셔너리로 정리 후 map_dfr로 반복 처리
  module_map <- list(
    list(data = df_tmj,   modules = c('CT', 'SUV', tmj_modules)),
    list(data = df_ent,   modules = c('CT', 'SUV', ent_modules)),
    list(data = df_hc,    modules = c('CT', 'SUV', hc_modules)),
    list(data = df_nonhc, modules = c('CT', 'SUV', modules))
  )
  
  y <- map_dfr(module_map, function(x) {
    crossing(x$data, module = x$modules)
  }) |>
    mutate(
      segfile = file.path(subjdir, paste0(module, "_", subjid, ".nii.gz"))
    )
  
  # seg 중 하나라도 없으면 어떤 파일이 없는지 알려주고 에러로 중단
  cat('Checking missing segmentation files..\n')
  missing <- y |> filter(!file.exists(segfile))
  
  if (nrow(missing) > 0) {
    print(missing)
    stop(nrow(missing), " segmentation file(s) missing")
  }
  
  cat('Checking 0-byte files..\n')
  zero_files <- y |> filter(file_size(segfile)==0)
  
  if (nrow(zero_files) > 0) {
    print(zero_files)
    stop(nrow(zero_files), " segmentation file(s) are zero-files")
  }
  
  df_seg_ok <- y |> distinct(subjdir, subjid_long, subjid)
  
  # 2. Sys.Date()는 이미 Date 클래스이므로 ymd() 불필요
  df_final <- df_seg_ok |>
    mutate(
      accn = str_extract(subjid, '\\d{10,13}'),
      pdate = str_remove(subjid_long, '[A-Z0-9]+_') |> ymd()
    )|>
    select(accn, pdate, subjid, subjid_long, subjdir, pdate) |>
    arrange(desc(pdate)) 
  
  df_final$qdate <- Sys.Date()
  
  write_csv(df_final, path(share_dir, 'seg_db.csv'))
  write_csv(df_final, path(home_dir,  'seg_db.csv'))
  
  return(df_final)
}
