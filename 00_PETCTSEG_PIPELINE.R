## -- 00. initialize --##
rm(list = ls())
source("~/CONFIG.R")
source("~/PETCTSEG/CONFIG_PETCTSEG.R")

# set working directory
src_dir <- "C:/Temp/PETCTSRC/"
dir_create(src_dir)

# quick update segmentation database
seg_db <- read_csv(path(share_dir, "seg_db.csv"), show_col_types = FALSE)
seg_db$pdate[1]

#-- 01. preprocessing --#
dcm_df <- read_dicom_dir(dcm_dir = 'E:/PETCTDCM/')

py <- "C:/Python313/python.exe"

processx::run(
  py, 
  args = c("petct_batch_pipeline.py", "--workers", "12", "--overwrite"),
  wd = script_dir,         # 작업 폴더 지정
  echo = TRUE,             # 실시간 출력
  error_on_status = TRUE)  # 파이썬이 실패하면 R에서도 에러

delete_deprecated(src_dir)

# -- 02. segmentation --#
prepare_segmentation(src_dir, task = "all")

#-- 03. roi_measurement and summary --#
measure_roi(src_dir, overwrite = TRUE)
summarize_subject_measures(overwrite = TRUE)

#-- 04. distribute_processed_files --#
distribute_processed_files(gz_dir)

#-- 05. update segmentation database --#
update_seg_db()
quick_db_uptake()
